% run a single ofdm symbol
% - sim, cc, rm, num: simulation parameters
% - snr_db: snr(db)
% - ch: channel parameter
% created: 2019.10.15
% modified:
% - channel coding added: 2019.10.21
% - rate matching added: 2019.10.22
% - structure updated: 2019.10.23
% - rate matching bug fixed: 2019.11.08
% - subframe buffer fixed: 2019.12.10
% - real channel estimation with single tone pilot added: 2020.02.09

function pkt_error = ofdm_single_run_r2(sim, cc, rm, num, snr_db, ch, chest_option, cheq_option)

% create turbo encoder/decoder
turbo_enc = comm.TurboEncoder('TrellisStructure', cc.tc_trellis, ...
    'InterleaverIndices', cc.PI+1);
turbo_dec = comm.TurboDecoder('TrellisStructure', cc.tc_trellis, ...
    'InterleaverIndices', cc.PI+1, 'NumIterations', cc.num_iter_max);

% create a rayleigh fading channel object
if ch.los
    fading_ch = comm.RicianChannel(...
        'SampleRate', num.sample_rate,...
        'PathDelays', ch.path_delays,...
        'AveragePathGains', ch.average_path_gains,...
        'KFactor', ch.k_factor,...
        'DirectPathDopplerShift', ch.maximum_doppler_shift,...
        'MaximumDopplerShift', ch.maximum_doppler_shift,...
        'DopplerSpectrum', ch.doppler_spectrum,...
        'PathGainsOutputPort', true, ...
        'NormalizePathGains', true);
    %     'DirectPathInitialPhase', 0.5,...
    %     'PathGainsOutputPort', true);
else
    fading_ch = comm.RayleighChannel(...
        'SampleRate', num.sample_rate, ...
        'PathDelays', ch.path_delays, ...
        'AveragePathGains', ch.average_path_gains, ...
        'NormalizePathGains', true, ...
        'MaximumDopplerShift', ch.maximum_doppler_shift, ...
        'DopplerSpectrum', ch.doppler_spectrum, ...
        'PathGainsOutputPort', true);
    %     'PathGainsOutputPort', true, ...
    %     'RandomStream', 'mt19937ar with seed', ...
    %     'Seed', 72, ...
    %     'Visualization', 'Impulse response');   % enables the impulse response channel visualization
end

% calculate noise variance
noise_var = (num.num_subc_usr/num.nfft)*(10^((-0.1)*snr_db));

% generate bit stream
tx_bit = randi([0 1], sim.len_tb_bit, 1);

% pad bits
tx_bit_pad = [tx_bit; zeros(cc.F, 1)];

% buffer per codeword
tx_bit_buff = reshape(tx_bit_pad, (sim.len_tb_bit+cc.F)/cc.C, cc.C);

% generate crc
tx_crc = crc.generator(cc.gCRC24A);
tx_bit_crc = generate(tx_crc, tx_bit_buff);

% turbo encode data per channel
tx_bit_enc = zeros(rm.D*3, cc.C);
for idx_cw = 1 : cc.C
    tx_bit_enc(:, idx_cw) = turbo_enc(tx_bit_crc(:, idx_cw));
end

% rate match
% tx_bit_ratematch = tx_ratematch(tx_bit_enc, rm);
tx_bit_ratematch = tx_ratematch_r2(tx_bit_enc, rm);

% modulate bit stream (only for 16qam)
tx_sym_serial = qammod(tx_bit_ratematch, 2^rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);

% simulate per subframe
tx_sym = reshape(tx_sym_serial, num.num_qamsym_usr, []);
rx_sym = zeros(size(tx_sym));
for idx_subfrm = 1 : size(tx_sym, 2)
    
    % extract subframe data per user
    tx_sym_data_subfrm = tx_sym(:, idx_subfrm);
    
    % map qam symbols to user physical resource blocks
    tx_sym_rbs = ofdm_sym_map(tx_sym_data_subfrm, num);
    
    % map user resource blocks to whole bandwidth
    tx_sym_bw = tx_sym_rbs;
    
    % map bandwidth to the center of fft range
    tx_sym_nfft_shift = zeros(num.nfft, num.num_ofdmsym_subfrm);
    tx_sym_nfft_shift((num.nfft/2)-(num.num_subc_usr/2)+1:(num.nfft/2)+(num.num_subc_usr/2), :) = tx_sym_bw;
    tx_sym_nfft = fftshift(tx_sym_nfft_shift, 1);
    
    % ofdm modulate
    tx_ofdmsym = sqrt(num.nfft) * ifft(tx_sym_nfft, [], 1);
    
    % add cp
    tx_ofdmsym_cp = tx_ofdmsym([num.nfft-num.num_cp+1:num.nfft 1:num.nfft], :);
    
    % serialize
    tx_ofdmsym_serial = tx_ofdmsym_cp(:);
    
    % pass signal through channel
    [tx_ofdmsym_faded, ch_path_gain] = fading_ch(tx_ofdmsym_serial);
    
    % add gaussian noise
    rx_ofdmsym_serial = awgn(tx_ofdmsym_faded, snr_db, 'measured');
    % rx_ofdm_sym_serial = awgn(tx_ofdm_sym_serial, snr_db, 'measured');
    
    % reshape
    rx_ofdmsym_cp = reshape(rx_ofdmsym_serial, num.nfft+num.num_cp, num.num_ofdmsym_subfrm);
    
    % remove cp
    rx_ofdmsym = rx_ofdmsym_cp(num.num_cp+1:end,:);
    
    % ofdm demodulate
    rx_sym_nfft = (1/sqrt(num.nfft)) * fft(rx_ofdmsym, [], 1);
    
    % demap symbols in bandwidth from fft range
    rx_sym_nfft_shift = fftshift(rx_sym_nfft, 1);
    rx_sym_bw = rx_sym_nfft_shift(num.nfft/2-num.num_subc_usr/2+1:num.nfft/2+num.num_subc_usr/2, :);
    
    % demap user resource block from whole bandwidth
    rx_sym_rbs = rx_sym_bw;
    
    % estimate channel (for uplink: channel estimation and equalization after demapping user physical resource block)
    ch_est_rbs = ofdm_ch_est(tx_sym_rbs, rx_sym_rbs, tx_ofdmsym_faded, num, chest_option, fading_ch, ch_path_gain);
    
    % equalize channel in tf-domain
    rx_sym_rbs_eq = otfs_ch_eq_tf(rx_sym_rbs, ch_est_rbs, noise_var, cheq_option);
    
    % demap data qam symbols
    rx_sym_data_subfrm = ofdm_sym_demap(rx_sym_rbs_eq, num);
    
    % buffer qam symbols
    rx_sym(:, idx_subfrm) = rx_sym_data_subfrm;
    
%     figure
%     subplot(2, 1, 1), plot(real(tx_sym_rbs(:)), '-b.'), hold on, plot(real(rx_sym_rbs_eq(:)), ':r.'), hold off
%     subplot(2, 1, 2), plot(imag(tx_sym_rbs(:)), '-b.'), hold on, plot(imag(rx_sym_rbs_eq(:)), ':r.'), hold off
%     figure
%     subplot(1, 2, 1), mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, real(rx_sym_rbs_eq-tx_sym_rbs))
%     subplot(1, 2, 2), mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, imag(rx_sym_rbs_eq-tx_sym_rbs))
%     fprintf('rmse:%10.4f\n', sqrt(mean((tx_sym_rbs(:)-rx_sym_rbs_eq(:)).^2)))
%     pause

end

% compensate channel estimation error variance
% qam_mse = mean(abs(tx_sym(:)-rx_sym(:)).^2);
if strcmp(chest_option, 'tf_ltedown') || strcmp(chest_option, 'tf_nr')
    error_var = noise_var + 0.12;
elseif strcmp(chest_option, 'tf_lteup')
    error_var = noise_var + 0.27;
elseif strcmp(chest_option, 'real')
    error_var = noise_var + 0.08;
else
    error_var = noise_var;
end

% demap the qam symbols
rx_bit_demod = (-1) * qamdemod(rx_sym(:), 2^rm.Qm, 'UnitAveragePower', true, ...
    'OutputType', 'approxllr', 'NoiseVariance', error_var);

% buffer per codeword
rx_bit_buff = reshape(rx_bit_demod, [], cc.C);  % rm.E*cc.C

% rate match
% rx_bit_ratematch = rx_ratematch(rx_bit_demod, rm);
rx_bit_ratematch = rx_ratematch_r2(rx_bit_buff, rm);

% turbo decode data
rx_bit_dec = zeros(cc.K, cc.C);
for idx_cw = 1 : cc.C
    rx_bit_dec(:, idx_cw) = turbo_dec(rx_bit_ratematch(:, idx_cw));
end

% detect crc
rx_crc = crc.detector(cc.gCRC24A);
[rx_bit_crc_removed, rx_crc_error] = detect(rx_crc, rx_bit_dec);

% assignin('base', 'rx_bit_dec', rx_bit_dec);
% assignin('base', 'rx_bit_pad', rx_bit_pad);
% assignin('base', 'rx_crc_error', rx_crc_error);
% pause

% serialize bits
rx_bit_pad = rx_bit_crc_removed(:);

% remove padded bits
rx_bit = rx_bit_pad(1 : sim.len_tb_bit);

if ~rx_crc_error
    pkt_error = symerr(tx_bit, rx_bit) > 0;
else
    pkt_error = true;
end

% fprintf('original noise variance: %10.4f    new noise variance: %10.4f\n', noise_var, qam_mse)
% figure, scatter(real(rx_sym_eq(:)), imag(rx_sym_eq(:)), 'r.'), hold on, scatter(real(tx_sym(:)), imag(tx_sym(:)), 'ko'), hold off, axis([-2 2 -2 2])
% figure
% subplot(2, 1, 1), plot(real(tx_sym(:)), '-b.'), hold on, plot(real(rx_sym(:)), ':r.'), hold off
% subplot(2, 1, 2), plot(imag(tx_sym(:)), '-b.'), hold on, plot(imag(rx_sym(:)), ':r.'), hold off
% fprintf('rmse:%10.4f\n', sqrt(mean((tx_sym(:)-rx_sym(:)).^2)))
% figure
% subplot(2, 1, 1), plot(tx_bit_ratematch(1:50), '-b.')
% subplot(2, 1, 2), plot(rx_bit_demod(1:50), '-b.')
% figure
% subplot(2, 1, 1), plot(real(tx_bit_enc(:)), '-b.'), hold on, plot(real(rx_bit_ratematch(:)), ':r.'), hold off
% subplot(2, 1, 2), plot(imag(tx_bit_enc(:)), '-b.'), hold on, plot(imag(rx_bit_ratematch(:)), ':r.'), hold off
% fprintf('rmse:%10.4f\n', sqrt(mean((tx_bit_enc(:)-rx_bit_ratematch(:)).^2)))
% pause

end
