function pkt_error = ofdm_simple(snr_db, cheq)

% % set parameters 1
% cheq = 'mmse';
% snr_db = 6;

% set parameters 2
len_tb_bit = 6120;  % tb length
num_subc_bw = 132;  % number of subcarriers in bandwidth
num_fft = 256;      % fft size
num_subc = 132;     % number of subcarriers per user
num_ofdmsym = 56;   % number of subcarriers per user
num_cp = 18;        % cp length
delay_spread_rms_us = 0.1;          % 100 ns
sample_rate = 15.36e6;
maximum_doppler_shift = 1.8519e3;   % 500 km/h
mcs = 8;

% set parameters 3
sim_option.override = true;     % override simulation options when true
sim_option.Qm = 4;              % number of bits per qam symbol
sim_option.coderate = 666;      % number of information bits per 1024 bits
num.num_qamsym_usr = num_subc*num_ofdmsym;
num.t_usrfrm = (num_fft+num_cp)*num_ofdmsym/sample_rate;
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
snr_db_adj = snr_db+(10*log10(num_subc_bw/num_fft));
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
tx_sym_rbs = reshape(tx_sym, num_subc, num_ofdmsym);

% map resource blocks (full)
tx_sym_bw = tx_sym_rbs;

% map bandwidth to the center of fft range
tx_sym_nfft_shift = zeros(num_fft, num_ofdmsym);
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
rx_ofdmsym_serial = awgn(tx_ofdmsym_faded, snr_db_adj, 'measured');

% regenerate real channel
[~, ch_real_fullmatrix, ~, ~, ~] = gen_real_ch_r1(fading_ch, ch_path_gain, num_fft, num_cp, num_subc_bw, num_ofdmsym, [], [], false, test_option);

% generate real one-tap time-frequency channel
ch_real_rbs = zeros(num_subc, num_ofdmsym);
for idx_ofdmsym = 1:num_ofdmsym
    ch_real_rbs(:, idx_ofdmsym) = diag(ch_real_fullmatrix(:, :, idx_ofdmsym));  % diagonal term only
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
rx_sym_rbs = rx_sym_bw;

% estimate channel
ch_est_rbs = ch_real_rbs;

% equalize channel
if strcmp(cheq, 'mmse')
    % calculate mmse channel
    ch_est_rbs_mmse = conj(ch_est_rbs)./(noise_var+abs(ch_est_rbs).^2);
    
    % equalize channel (one-tap)
    rx_sym_rbs_eq = rx_sym_rbs.*ch_est_rbs_mmse;
    
    % calculate noise variance
%     noise_var_rbs = noise_var*(abs(ch_est_rbs_mmse).^2);
    noise_var_rbs = noise_var./(abs(ch_est_rbs).^2);
else
    % equalize channel (one-tap zf)
    rx_sym_rbs_eq = rx_sym_rbs./ch_est_rbs;
    
    % calculate noise variance
    noise_var_rbs = noise_var./(abs(ch_est_rbs).^2);
end

% demap data symbols
rx_sym = rx_sym_rbs_eq(:);

% calculate error variance
noise_var = noise_var_rbs(:);

% demap qam symbols
rx_bit_demod = (-1) * qamdemod(rx_sym.', 2^rm.Qm, 'UnitAveragePower', true, ...
    'OutputType', 'approxllr', 'NoiseVariance', noise_var.');
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
