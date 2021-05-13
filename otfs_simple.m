function pkt_error = otfs_simple(snr_db, cheq)

% % set parameters 1
% cheq = 'mmse';
% snr_db = 8;

% set parameters 2
len_tb_bit = 6120;  % tb length
num_subc_bw = 132;  % number of subcarriers in bandwidth
num_fft = 256;      % fft size
num_delay = 132;     % number of subcarriers per user
num_doppler = 56;   % number of subcarriers per user
num_cp = 18;        % cp length
delay_spread_rms_us = 0.1;          % 100 ns
sample_rate = 15.36e6;
maximum_doppler_shift = 1.8519e3;   % 500 km/h
mcs = 8;

% set parameters 3
sim_option.override = true;     % override simulation options when true
sim_option.Qm = 4;              % number of bits per qam symbol
sim_option.coderate = 666;      % number of information bits per 1024 bits
num.num_qamsym_usr = num_delay*num_doppler;
num.t_usrfrm = (num_fft+num_cp)*num_doppler/sample_rate;
test_option.gpu_flag = false;
cc = nw_cc_prm(len_tb_bit);
rm = nw_rm_prm(len_tb_bit, mcs, num, cc, sim_option);

% create objects
tx_crc = comm.CRCGenerator(cc.gCRC24A);
rx_crc = comm.CRCDetector(cc.gCRC24A);
turbo_enc = comm.TurboEncoder('TrellisStructure', cc.tc_trellis, ...
    'InterleaverIndices', cc.PI+1);
turbo_dec = comm.TurboDecoder('TrellisStructure', cc.tc_trellis, ...
    'InterleaverIndices', cc.PI+1, 'NumIterations', cc.num_iter_max);
fading_ch = nrTDLChannel( ...
    'DelayProfile', 'TDL-C', ...
    'DelaySpread', delay_spread_rms_us*1e-6, ...
    'MaximumDopplerShift', maximum_doppler_shift, ...
    'NumTransmitAntennas', 1, ...
    'NumReceiveAntennas', 1, ...
    'SampleRate', sample_rate, ...
    'NormalizePathGains', true);

% calculate snr and noise variance
% snr_db_adj = snr_db+(10*log10(num_subc_bw/num_fft));
noise_var = 10^((-0.1)*snr_db);

%% tx

% generate bit stream
% rng(0)
tx_bit = randi([0 1], len_tb_bit, 1);

% pad bits
tx_bit_pad = [tx_bit; zeros(cc.F, 1)];

% buffer per code block
tx_bit_buff = reshape(tx_bit_pad, (len_tb_bit+cc.F), 1);

% generate crc
tx_bit_crc = tx_crc(tx_bit_buff);

% turbo encode data
tx_bit_enc = turbo_enc(tx_bit_crc);

% rate match
tx_bit_ratematch = tx_ratematch_r2(tx_bit_enc, 1, rm);

% modulate bit stream
tx_sym_serial = qammod(tx_bit_ratematch(:), 2^rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);

% simulate per user frame
tx_sym = reshape(tx_sym_serial, num.num_qamsym_usr, []);

% map qam symbols to user physical resource blocks
tx_sym_rbs = reshape(tx_sym, num_delay, num_doppler);

% 2d sfft
tx_sym_rbs_tf = sqrt(num_doppler/num_delay)*fft(ifft(tx_sym_rbs, [], 2), [], 1);

% map resource blocks (full)
tx_sym_bw = tx_sym_rbs_tf;

% map bandwidth to the center of fft range
tx_sym_nfft_shift = zeros(num_fft, num_doppler);
tx_sym_nfft_shift((num_fft/2)-(num_subc_bw/2)+1:(num_fft/2)+(num_subc_bw/2), :) = tx_sym_bw;
tx_sym_nfft = fftshift(tx_sym_nfft_shift, 1);

% ofdm modulate
tx_ofdmsym = sqrt(num_fft)*ifft(tx_sym_nfft, [], 1);

% add cp
tx_ofdmsym_cp = tx_ofdmsym([num_fft-num_cp+1:num_fft 1:num_fft], :);

% serialize
tx_ofdmsym_serial = tx_ofdmsym_cp(:);

%% channel

% pass signal through channel
% rng(0)
[tx_ofdmsym_faded, ch_path_gain] = fading_ch(squeeze(tx_ofdmsym_serial));

% add gaussian noise
% rng(0)
% rx_ofdmsym_serial = awgn(tx_ofdmsym_faded, snr_db_adj, 'measured');
% rx_ofdmsym_serial = awgn(tx_ofdmsym_faded, snr_db_adj, 1);
% rx_ofdmsym_serial = awgn(tx_ofdmsym_faded, snr_db, 10*log10(num_subc_bw/num_fft));
rx_ofdmsym_serial = awgn(tx_ofdmsym_faded, snr_db, 0);

% regenerate real channel
[~, ch_real_fullmatrix, ~, ~, ~] = gen_real_ch_r1(fading_ch, ch_path_gain, num_fft, num_cp, num_subc_bw, num_doppler, [], [], false, test_option);

% generate real one-tap time-frequency channel matrix
ch_real_rbs = zeros(num_delay, num_doppler);
for idx_ofdmsym = 1:num_doppler
    ch_real_rbs(:, idx_ofdmsym) = diag(ch_real_fullmatrix(:, :, idx_ofdmsym));  % diagonal term only
end

% generate effective time-frequency channel matrix
ch_real_eff = zeros(num_delay*num_doppler);
for idx_ofdmsym = 1:num_doppler
    ch_real_eff((idx_ofdmsym-1)*num_delay+1:idx_ofdmsym*num_delay, (idx_ofdmsym-1)*num_delay+1:idx_ofdmsym*num_delay) = ch_real_fullmatrix(:, :, idx_ofdmsym);
end

%% rx

% reshape
rx_ofdmsym_cp = reshape(rx_ofdmsym_serial, num_fft+num_cp, []);

% remove cp
rx_ofdmsym = rx_ofdmsym_cp(num_cp+1:end,:);

% ofdm demodulate
rx_sym_nfft = (1/sqrt(num_fft))*fft(rx_ofdmsym, [], 1);

% demap symbols in bandwidth from fft range
rx_sym_nfft_shift = fftshift(rx_sym_nfft, 1);
rx_sym_bw = rx_sym_nfft_shift((num_fft/2)-(num_subc_bw/2)+1:(num_fft/2)+(num_subc_bw/2), :);

% demap user resource block from whole bandwidth (full)
rx_sym_rbs_tf = rx_sym_bw;

% estimate channel
ch_est_rbs_tf = ch_real_rbs;

% x = tx_sym_rbs_tf(:);
% y = rx_sym_rbs_tf(:);
% h = ch_real_eff;
% hx = h*x;
% w = 1./ch_est_rbs_tf(:);
% w_mmse = conj(ch_est_rbs_tf(:))./(noise_var+abs(ch_est_rbs_tf(:)).^2);
% wy = w.*y;
% wy_mmse = w_mmse.*y;
% 
% wy_dd = sqrt(num_delay/num_doppler)*fft(ifft(reshape(wy, 132, 56), [], 1), [], 2);
% whx_dd = sqrt(num_delay/num_doppler)*fft(ifft(reshape(diag(w)*hx, 132, 56), [], 1), [], 2);
% wy_mmse_dd = sqrt(num_delay/num_doppler)*fft(ifft(reshape(wy_mmse, 132, 56), [], 1), [], 2);
% whx_mmse_dd = sqrt(num_delay/num_doppler)*fft(ifft(reshape(diag(w_mmse)*hx, 132, 56), [], 1), [], 2);
% 
% figure(1)
% subplot(2, 3, 1), scatter(real(wy_dd(:)), imag(wy_dd(:))), grid, legend('zf rx signal with noise'), axis([-5 5 -5 5])
% subplot(2, 3, 2), scatter(real(whx_dd(:)), imag(whx_dd(:))), grid, legend('zf rx signal without noise'), axis([-5 5 -5 5])
% subplot(2, 3, 3), scatter(real(wy_dd(:)-whx_dd(:)), imag(wy_dd(:)-whx_dd(:))), grid, legend('zf noise'), axis([-5 5 -5 5])
% subplot(2, 3, 4), scatter(real(wy_mmse_dd(:)), imag(wy_mmse_dd(:))), grid, legend('mmse rx signal with noise'), axis([-5 5 -5 5])
% subplot(2, 3, 5), scatter(real(whx_mmse_dd(:)), imag(whx_mmse_dd(:))), grid, legend('mmse rx signal without noise'), axis([-5 5 -5 5])
% subplot(2, 3, 6), scatter(real(wy_mmse_dd(:)-whx_mmse_dd(:)), imag(wy_mmse_dd(:)-whx_mmse_dd(:))), grid, legend('mmse noise'), axis([-5 5 -5 5])
% pkt_error = 0;
% fprintf('zf noise pwr: %8.4f, %8.4f (dB)\n', mean(abs(wy_dd-whx_dd).^2, 'all'), 10*log10(mean(abs(wy_dd-whx_dd).^2, 'all')))
% fprintf('mmse noise pwr: %8.4f, %8.4f (dB)\n', mean(abs(wy_mmse_dd-whx_mmse_dd).^2, 'all'), 10*log10(mean(abs(wy_mmse_dd-whx_mmse_dd).^2, 'all')))
% ch_est_rbs_dd = sqrt(num_delay/num_doppler)*fft(ifft(1./ch_est_rbs_tf, [], 1), [], 2);
% error_var = noise_var*mean(abs(sum(gen_eff_ch(ch_est_rbs_dd), 2)).^2, 'all');
% fprintf('est noise pwr: %8.4f, %8.4f (dB)\n', error_var, 10*log10(error_var))
% pause

% equalize channel
if strcmp(cheq, 'mmse')
    % calculate mmse channel
    ch_est_rbs_mmse_tf = conj(ch_est_rbs_tf)./(noise_var+abs(ch_est_rbs_tf).^2);
    
    % equalize channel (one-tap)
    rx_sym_rbs_eq_tf = rx_sym_rbs_tf.*ch_est_rbs_mmse_tf;
    
%     % calculate noise variance (original 1, performs well)
%     ch_est_rbs_mmse_dd = sqrt(num_delay/num_doppler)*fft(ifft(ch_est_rbs_mmse_tf, [], 1), [], 2);
%     ch_inv_dd = reshape(sum(gen_eff_ch(ch_est_rbs_mmse_dd), 2), num_delay, num_doppler);
%     noise_var_rbs = noise_var*(abs(ch_inv_dd).^2);
    
%     % calculate noise variance (original 2, performs well)
%     ch_est_rbs_mmse_dd = sqrt(num_delay/num_doppler)*fft(ifft(1./ch_est_rbs_tf, [], 1), [], 2);
%     ch_inv_dd = reshape(sum(gen_eff_ch(ch_est_rbs_mmse_dd), 2), num_delay, num_doppler);
%     noise_var_rbs = noise_var*(abs(ch_inv_dd).^2);
    
%     % calculate noise variance (original 3, performs well)
%     ch_est_rbs_mmse_dd = sqrt(num_delay/num_doppler)*fft(ifft(1./ch_est_rbs_tf, [], 1), [], 2);
%     ch_inv_dd = reshape(sum(abs(gen_eff_ch(ch_est_rbs_mmse_dd)).^2, 2), num_delay, num_doppler);
%     noise_var_rbs = noise_var*ch_inv_dd;
else
    % equalize channel (one-tap zf)
    rx_sym_rbs_eq_tf = rx_sym_rbs_tf./ch_est_rbs_tf;
end

% 2d invers sft
rx_sym_rbs = sqrt(num_delay/num_doppler)*fft(ifft(rx_sym_rbs_eq_tf, [], 1), [], 2);

% demap data symbols
rx_sym = rx_sym_rbs(:);

% calculate error variance
% error_var = noise_var_rbs(:);
% error_var = noise_var*mean(abs(1./ch_est_rbs_tf).^2, 'all');
% sfft_mtx = kron(dftmtx(num_doppler), conj(dftmtx(num_delay))/num_delay);
% isfft_mtx = kron(conj(dftmtx(num_doppler))/num_doppler, dftmtx(num_delay));
% error_var = noise_var*abs(sum(sfft_mtx*diag(1./ch_est_rbs_tf(:))*isfft_mtx, 2)).^2;
ch_est_rbs_dd = sqrt(num_delay/num_doppler)*fft(ifft(1./ch_est_rbs_tf, [], 1), [], 2);
error_var = noise_var*mean(abs(sum(gen_eff_ch(ch_est_rbs_dd), 2)).^2, 'all');
% error_var = 10*error_var;




% demap qam symbols
rx_bit_demod = (-1) * qamdemod(rx_sym.', 2^rm.Qm, 'UnitAveragePower', true, ...
    'OutputType', 'approxllr', 'NoiseVariance', error_var.');
rx_bit_demod = rx_bit_demod(:);

% rate match
rx_bit_ratematch = rx_ratematch_r2(rx_bit_demod, rm);

% turbo decode data
rx_bit_dec = turbo_dec(rx_bit_ratematch);

% detect crc
[rx_bit_crc_removed, rx_crc_error] = rx_crc(rx_bit_dec);

% buffer per code block
rx_bit_buff = reshape(rx_bit_crc_removed, [], 1);

% remove padded bits
rx_bit = rx_bit_buff(1:len_tb_bit);

% rearrange crc
rx_crc_error = sum(double(rx_crc_error))>0;

% calculate packet error
pkt_error = rx_crc_error | (symerr(tx_bit, rx_bit, 'column-wise') > 0);

end
