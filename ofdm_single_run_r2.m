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
% - slot buffer fixed: 2019.12.10
% - real channel estimation with single tone pilot added: 2020.02.09
% - slot-based to user-frame-based simulation (user frame = multiple slots)

function [pkt_error, tx_papr, ch_mse, sym_err_var] = ofdm_single_run_r2(sim, cc, rm, num, snr_db, ch, chest_option, cheq_option, test_option)

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
noise_var = (num.num_subc_usr/num.num_fft)*(10^((-0.1)*snr_db));

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

% modulate bit stream
tx_sym_serial = qammod(tx_bit_ratematch(:), 2^rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);

% simulate per user frame (multiple slots per user frame)
tx_sym = reshape(tx_sym_serial, num.num_qamsym_usr, []);
rx_sym = zeros(size(tx_sym));
sym_pilot_error = zeros(num.num_subc_pilot_usr*num.num_ofdmsym_pilot_usr, size(tx_sym, 2));
tx_papr_usrfrm = zeros(1, size(tx_sym, 2));     % for test: papr
ch_mse_usrfrm = zeros(1, size(tx_sym, 2));      % for test: channel mse
for idx_usrfrm = 1 : size(tx_sym, 2)
    
    % extract user frame data per user
    tx_sym_data_usrfrm = tx_sym(:, idx_usrfrm);
    
    % map qam symbols to user physical resource blocks
    [tx_sym_rbs, tx_sym_pilot_usrfrm] = ofdm_sym_map(tx_sym_data_usrfrm, num);
    
    % map user resource blocks to whole bandwidth (temporary)
    tx_sym_bw = zeros(num.num_subc_bw, num.num_ofdmsym_usr);
    tx_sym_bw(1:size(tx_sym_rbs, 1), 1:size(tx_sym_rbs, 2)) = tx_sym_rbs;
    
    % map bandwidth to the center of fft range
    tx_sym_nfft_shift = zeros(num.num_fft, num.num_ofdmsym_usr);
    tx_sym_nfft_shift((num.num_fft/2)-(num.num_subc_bw/2)+1:(num.num_fft/2)+(num.num_subc_bw/2), :) = tx_sym_bw;
    tx_sym_nfft = fftshift(tx_sym_nfft_shift, 1);
    
    % ofdm modulate
    tx_ofdmsym = sqrt(num.num_fft) * ifft(tx_sym_nfft, [], 1);
    
    % add cp
    tx_ofdmsym_cp = tx_ofdmsym([num.num_fft-num.num_cp+1:num.num_fft 1:num.num_fft], :);
    
    % test: calculate papr
    if test_option.papr
        tx_papr_usrfrm(1, idx_usrfrm) = mean(max(abs(tx_ofdmsym_cp).^2, [], 1)./mean(abs(tx_ofdmsym_cp).^2, 1));
    end
    
    % serialize
    tx_ofdmsym_serial = tx_ofdmsym_cp(:);
    
    % pass signal through channel
    [tx_ofdmsym_faded, ch_path_gain] = fading_ch(tx_ofdmsym_serial);
    
    % add gaussian noise
    rx_ofdmsym_serial = awgn(tx_ofdmsym_faded, snr_db, 'measured');
%     rx_ofdmsym_serial = awgn(tx_ofdmsym_serial, snr_db, 'measured');
    
    % reshape
    rx_ofdmsym_cp = reshape(rx_ofdmsym_serial, num.num_fft+num.num_cp, []);
    
    % remove cp
    rx_ofdmsym = rx_ofdmsym_cp(num.num_cp+1:end, :);
    
    % ofdm demodulate
    rx_sym_nfft = (1/sqrt(num.num_fft)) * fft(rx_ofdmsym, [], 1);
    
    % demap symbols in bandwidth from fft range
    rx_sym_nfft_shift = fftshift(rx_sym_nfft, 1);
    rx_sym_bw = rx_sym_nfft_shift((num.num_fft/2)-(num.num_subc_bw/2)+1:(num.num_fft/2)+(num.num_subc_bw/2), :);
    
    % demap user resource block from whole bandwidth (temporary)
    rx_sym_rbs = rx_sym_bw(1:num.num_subc_usr, 1:num.num_ofdmsym_usr);
    
    % estimate channel (for uplink: channel estimation and equalization after demapping user physical resource block)
    ch_est_rbs = ofdm_ch_est(tx_sym_rbs, rx_sym_rbs, tx_ofdmsym_faded, num, chest_option, fading_ch, ch_path_gain);
    
    % test: calculate channel estimation error
    if test_option.ch_mse
        ch_real_rbs = ofdm_ch_est(tx_sym_rbs, rx_sym_rbs, tx_ofdmsym_faded, num, 'real', fading_ch, ch_path_gain);
        ch_mse_usrfrm(1, idx_usrfrm) = mean(abs(ch_real_rbs-ch_est_rbs).^2, 'all');
    end
    
    % equalize channel in tf-domain
    rx_sym_rbs_eq = ofdm_ch_eq(rx_sym_rbs, ch_est_rbs, noise_var, cheq_option);
    
    % demap data and pilot qam symbols
    [rx_sym_data_usrfrm, rx_sym_pilot_usrfrm] = ofdm_sym_demap(rx_sym_rbs_eq, num);
%     rx_sym_data_usrfrm = ofdm_sym_demap(rx_sym_rbs, num);
    
    % buffer qam symbols
    rx_sym(:, idx_usrfrm) = rx_sym_data_usrfrm;
    
    % buffer pilot errors
    sym_pilot_error(:, idx_usrfrm) = rx_sym_pilot_usrfrm(:)-tx_sym_pilot_usrfrm(:);   % for error variance calculation
    
%     fprintf('tx_data_pwr:%6.3f    rx_data_pwr:%6.3f    tx_pilot_pwr:%6.3f    rx_pilot_pwr:%6.3f\n', std(tx_sym(:)), std(rx_sym(:)), std(tx_sym_pilot_usrfrm(:)), std(rx_sym_pilot_usrfrm(:)))
%     fprintf('tx_data_pwr:%6.3f    rx_data_pwr:%6.3f    tx_pilot_pwr:%6.3f    rx_pilot_pwr:%6.3f\n', mean(abs(tx_sym(:)).^2), mean(abs(rx_sym(:)).^2), mean(abs(tx_sym_pilot_usrfrm(:)).^2), mean(abs(rx_sym_pilot_usrfrm(:)).^2))
%     fprintf('data_std:%6.3f    pilot_std:%6.3f\n', mean(abs(tx_sym(:)-rx_sym(:)).^2), mean(abs(tx_sym_pilot_usrfrm(:)-rx_sym_pilot_usrfrm(:)).^2))
%     fprintf('data_std:%6.3f    pilot_std:%6.3f\n', std(tx_sym(:)-rx_sym(:)), std(tx_sym_pilot_usrfrm(:)-rx_sym_pilot_usrfrm(:)))
%     pause
    
%     figure
%     subplot(2, 1, 1), plot(real(tx_sym_rbs(:)), '-b.'), hold on, plot(real(rx_sym_rbs_eq(:)), ':r.'), hold off
%     subplot(2, 1, 2), plot(imag(tx_sym_rbs(:)), '-b.'), hold on, plot(imag(rx_sym_rbs_eq(:)), ':r.'), hold off
%     figure
%     subplot(1, 2, 1), mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, real(rx_sym_rbs_eq-tx_sym_rbs))
%     subplot(1, 2, 2), mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, imag(rx_sym_rbs_eq-tx_sym_rbs))
%     fprintf('rmse:%10.4f\n', sqrt(mean((tx_sym_rbs(:)-rx_sym_rbs_eq(:)).^2)))
%     pause

end

% serialize qam symbols
rx_sym_serial = rx_sym(:);

% compensate channel estimation error variance
% qam_mse = mean(abs(tx_sym(:)-rx_sym(:)).^2);
% error_var = var(sym_pilot_error(:));    % cannot be used (pilot syms have less error than data syms)
% error_var = var(tx_sym_serial-rx_sym_serial);
if strcmp(chest_option, 'tf_ltedown') || strcmp(chest_option, 'tf_nr')
    error_var = noise_var + 0.12;
elseif strcmp(chest_option, 'tf_lteup')
    error_var = noise_var + 0.27;
elseif strcmp(chest_option, 'real')
    error_var = noise_var + 0.08;
else
    error_var = noise_var;
end

% test channel estimation error variance
if test_option.sym_err_var
    sym_err_var = var(tx_sym_serial-rx_sym_serial);
    fprintf('noise_var:%6.3f    sym_err_var:%6.3f    calc_err_var:%6.3f\n', noise_var, sym_err_var, error_var)
%     figure
%     subplot(2, 1, 1), plot(real(tx_sym_serial), '-b.'), hold on, plot(real(rx_sym_serial), '-r'), hold off, grid minor
%     subplot(2, 1, 2), plot(imag(tx_sym_serial), '-b.'), hold on, plot(imag(rx_sym_serial), '-r'), hold off, grid minor
    pause
else
    sym_err_var = 0;
end

% demap the qam symbols
rx_bit_demod = (-1) * qamdemod(rx_sym_serial, 2^rm.Qm, 'UnitAveragePower', true, ...
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

% serialize bits
rx_bit_pad = rx_bit_crc_removed(:);

% remove padded bits
rx_bit = rx_bit_pad(1 : sim.len_tb_bit);

% calculate packet error
if ~rx_crc_error
    pkt_error = symerr(tx_bit, rx_bit) > 0;
else
    pkt_error = true;
end

% calculate papr
if test_option.papr
    tx_papr = mean(tx_papr_usrfrm);
else
    tx_papr = [];
end

% calculate channel mse
if test_option.ch_mse
    ch_mse = mean(ch_mse_usrfrm);
else
    ch_mse = [];
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
