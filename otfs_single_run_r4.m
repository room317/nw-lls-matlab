% run a single otfs symbol
% - sim, cc, rm, num: simulation parameters
% - snr_db: snr(db)
% - ch: channel parameter
% - eq_dd: delay-doppler channel estimation (true/false)
% - test_dd_pilot: whole resource plane is allocated to pilot, for test only (true/false)
% created: 2019.10.15
% modified:
% - channel coding added: 2019.10.21
% - rate matching added: 2019.10.22
% - structure updated: 2019.10.23
% - rate matching bug fixed: 2019.11.08
% - slot buffer fixed: 2019.12.10
% - eq_dd(delay-doppler channel estimation) added: 2020.01.09
% - test_dd_pilot added: 2020.01.09
% - real channel estimation with single tone pilot added: 2020.02.19
% - real channel compensation in tf domain added: 2020.02.19
% - real channel compensation in dd domain added: 2020.02.19
% - variable name changed: 2020.06.09
% - slot-based to user-frame-based simulation (user frame = multiple slots)
% - full-tap real channel regeneration added: 2020.0908

function [pkt_error, tx_papr, ch_mse, sym_err_var] = otfs_single_run_r4(sim, cc, rm, num, snr_db, ch, chest_option, cheq_option, test_option)

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
        'NormalizePathGains', true, ...
        'KFactor', ch.k_factor,...
        'DirectPathDopplerShift', ch.maximum_doppler_shift,...
        'MaximumDopplerShift', ch.maximum_doppler_shift,...
        'DopplerSpectrum', ch.doppler_spectrum, ...
        'PathGainsOutputPort', true);
%         'DirectPathInitialPhase', 0.5,...
%         'PathGainsOutputPort', true);
else
    fading_ch = comm.RayleighChannel(...
        'SampleRate', num.sample_rate, ...
        'PathDelays', ch.path_delays, ...
        'AveragePathGains', ch.average_path_gains, ...
        'NormalizePathGains', true, ...
        'MaximumDopplerShift', ch.maximum_doppler_shift, ...
        'DopplerSpectrum', ch.doppler_spectrum,...
        'PathGainsOutputPort', true);
%         'RandomStream', 'mt19937ar with seed', ...
%         'Seed', ch_param.seed);
%         'PathGainsOutputPort', true);
%         'Visualization', 'Impulse response');   % enables the impulse response channel visualization

end

% calculate noise variance
noise_var = (num.num_subc_usr/num.num_fft)*(10 ^ ((-0.1)*snr_db));

% generate bit stream
tx_bit = randi([0 1], sim.len_tb_bit, 1);

% pad bits
tx_bit_pad = [tx_bit; zeros(cc.F, 1)];

% buffer per codeword
tx_bit_buff = reshape(tx_bit_pad, (sim.len_tb_bit+cc.F)/cc.C, cc.C);

% generate crc
tx_crc = crc.generator(cc.gCRC24A);
tx_bit_crc = generate(tx_crc, tx_bit_buff);

% turbo encode data
tx_bit_enc = zeros(rm.D*3, cc.C);
for idx_cw = 1 : cc.C
    tx_bit_enc(:, idx_cw) = turbo_enc(tx_bit_crc(:, idx_cw));
end

% rate match
tx_bit_ratematch = tx_ratematch_r2(tx_bit_enc, rm);

% modulate bit stream
tx_sym_serial = qammod(tx_bit_ratematch(:), 2^rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);

% simulate per user frame
tx_sym = reshape(tx_sym_serial, num.num_qamsym_usr, []);
rx_sym = zeros(size(tx_sym));
tx_papr_usrfrm = zeros(1, size(tx_sym, 2));     % for test: papr
ch_mse_usrfrm = zeros(1, size(tx_sym, 2));      % for test: channel mse
% sym_pilot_error = ...
%     zeros((num.num_delay_pilot_usr-2*num.num_delay_guard_usr)*(num.num_doppler_pilot_usr-2*num.num_doppler_guard_usr), ...
%     size(tx_sym, 2));                           % for pilot noise estimation (now working)
for idx_usrfrm = 1 : size(tx_sym, 2)
    
    % extract user frame data per user
    tx_sym_data_usrfrm = tx_sym(:, idx_usrfrm);
    
    % map data and pilot symbols
    if strcmp(chest_option, 'dd_tone')
%         [tx_sym_rbs_dd, tx_sym_pilot_usrfrm] = otfs_sym_map_r2(tx_sym_data_usrfrm, num, test_option);       % for pilot noise estimation (now working)
        [tx_sym_rbs_dd, ~] = otfs_sym_map_r2(tx_sym_data_usrfrm, num, test_option);
    else
        tx_sym_rbs_dd = reshape(tx_sym_data_usrfrm, num.num_delay_usr, []);
    end
    
    % 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
    tx_sym_rbs_tf = sqrt(num.num_ofdmsym_usr/num.num_subc_usr)*fft(ifft(tx_sym_rbs_dd, [], 2), [], 1);
    
    % map user resource blocks to whole bandwidth (temporary)
    tx_sym_bw = zeros(num.num_subc_bw, num.num_ofdmsym_usr);
    tx_sym_bw(1:size(tx_sym_rbs_tf, 1), 1:size(tx_sym_rbs_tf, 2)) = tx_sym_rbs_tf;
    
    % map bandwidth to the center of fft range
    tx_sym_nfft_shift = zeros(num.num_fft, num.num_ofdmsym_usr);
    tx_sym_nfft_shift((num.num_fft/2)-(num.num_subc_bw/2)+1:(num.num_fft/2)+(num.num_subc_bw/2), :) = tx_sym_bw;
    tx_sym_nfft = fftshift(tx_sym_nfft_shift, 1);
    
    % ofdm modulate
    tx_ofdmsym = sqrt(num.num_fft) * ifft(tx_sym_nfft, [], 1);
    
%     % test: walsh-hadamard
% %     % generate walsh-hadamard sequence
% %     spread_seq = 1;
% %     for idx = 1 : log2(16)
% %         spread_seq_update = kron([1 1; 1 -1], spread_seq);
% %         spread_seq = spread_seq_update;
% %     end
% 
%     tx_ofdmsym_wh = fwht(tx_ofdmsym);
%     
%     figure
%     subplot(2, 2, 1), mesh(1:size(tx_ofdmsym, 2), 1:size(tx_ofdmsym, 1), real(tx_ofdmsym)), axis_x = axis;
%     subplot(2, 2, 2), mesh(1:size(tx_ofdmsym, 2), 1:size(tx_ofdmsym, 1), imag(tx_ofdmsym)), axis_y = axis;
%     subplot(2, 2, 3), mesh(1:size(tx_ofdmsym_wh, 2), 1:size(tx_ofdmsym_wh, 1), real(tx_ofdmsym_wh)), axis(axis_x)
%     subplot(2, 2, 4), mesh(1:size(tx_ofdmsym_wh, 2), 1:size(tx_ofdmsym_wh, 1), imag(tx_ofdmsym_wh)), axis(axis_y)
%     pause
    
%     % Scrambler
%     scrambler_in = [header(8 : end) header_hcs data];
%     scrambler_out = lfsr(scrambler_in, SIM_CONFIG.HEADER_SCRAMBLER_INITIALIZATION);
%     scrambled_header = [header(1 : 7) scrambler_out(1 : 57)];
%     scrambled_data = scrambler_out(58 : end);
%     
%     % Descrambler
%     descrambler_in = [rx_header_zero_pad_removed(8 : end) rx_data_zero_pad_removed];
%     descrambler_out = lfsr(descrambler_in, SIM_CONFIG.HEADER_SCRAMBLER_INITIALIZATION);
%     descrambled_header = [rx_header_zero_pad_removed(1 : 7) descrambler_out(1 : 57)];
%     descrambled_data = descrambler_out(58 : end);
    
    % add cp
    tx_ofdmsym_cp = tx_ofdmsym([num.num_fft-num.num_cp+1:num.num_fft 1:num.num_fft], :);
    
    % test: calculate papr
    if test_option.papr
        tx_papr_usrfrm(1, idx_usrfrm) = mean(max(abs(tx_ofdmsym_cp).^2, [], 1)./mean(abs(tx_ofdmsym_cp).^2, 1));
    end
    
%     figure, mesh(1:size(tx_ofdmsym_cp, 2), 1:size(tx_ofdmsym_cp, 1), abs(tx_ofdmsym_cp))
%     pause
    
    % serialize
    tx_ofdmsym_serial = tx_ofdmsym_cp(:);
    
%     assignin('base', 'tx_sym_data_usrfrm', tx_sym_data_usrfrm)
%     assignin('base', 'tx_sym_rbs_dd', tx_sym_rbs_dd)
%     assignin ('base', 'tx_ofdmsym_serial', tx_ofdmsym_serial)
%     pause
    
    % pass signal through channel
    [tx_ofdmsym_faded, ch_path_gain] = fading_ch(tx_ofdmsym_serial);
    
    % add gaussian noise
    rx_ofdmsym_serial = awgn(tx_ofdmsym_faded, snr_db, 'measured');
    
%     fprintf('tx_sig: %6.3f    tx_sig_awgn: %6.3f\n', mean(abs(tx_ofdmsym_serial).^2, 'all'), mean(abs(rx_ofdmsym_serial).^2, 'all'))
%     pause
    
    % reshape
    rx_ofdmsym_cp = reshape(rx_ofdmsym_serial, num.num_fft+num.num_cp, num.num_ofdmsym_usr);
%     rx_ofdmsym_cp = reshape(tx_ofdmsym_serial, num.num_fft+num.num_cp, num.num_ofdmsym_usr);
    
    % remove cp
    rx_ofdmsym = rx_ofdmsym_cp(num.num_cp+1:end,:);
    
    % remove parts of received symbols
    if isscalar(test_option.partial_reception) && ismember(test_option.partial_reception, 1:num.num_ofdmsym_usr-1)
        rx_ofdmsym_part = [rx_ofdmsym(:, 1:test_option.partial_reception) zeros(size(rx_ofdmsym, 1), num.num_ofdmsym_usr-test_option.partial_reception)];
    else
        rx_ofdmsym_part = rx_ofdmsym;
    end
    
%     % test: walsh-hadamard
%     rx_ofdmsym_wh = ifwht(rx_ofdmsym);
%     
%     figure
%     subplot(2, 2, 1), mesh(1:size(rx_ofdmsym, 2), 1:size(rx_ofdmsym, 1), real(rx_ofdmsym))
%     subplot(2, 2, 2), mesh(1:size(rx_ofdmsym, 2), 1:size(rx_ofdmsym, 1), imag(rx_ofdmsym))
%     subplot(2, 2, 3), mesh(1:size(rx_ofdmsym_wh, 2), 1:size(rx_ofdmsym_wh, 1), real(rx_ofdmsym_wh))
%     subplot(2, 2, 4), mesh(1:size(rx_ofdmsym_wh, 2), 1:size(rx_ofdmsym_wh, 1), imag(rx_ofdmsym_wh))
%     pause
    
    % ofdm demodulate
    rx_sym_nfft = (1/sqrt(num.num_fft)) * fft(rx_ofdmsym_part, [], 1);
    
    % demap symbols in bandwidth from fft range
    rx_sym_nfft_shift = fftshift(rx_sym_nfft, 1);
    rx_sym_bw = rx_sym_nfft_shift((num.num_fft/2)-(num.num_subc_bw/2)+1:(num.num_fft/2)+(num.num_subc_bw/2), :);
    
    % demap user resource block from whole bandwidth (temporary)
    rx_sym_rbs_tf = rx_sym_bw(1:num.num_delay_usr, 1:num.num_doppler_usr);
    
    % estimate and equalize channel
    if strcmp(chest_option, 'dd_tone') && strcmp(cheq_option, 'ddeq')
        
        % 2d inverse sfft both for channel estimation and demodulation
        rx_sym_rbs_dd = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(rx_sym_rbs_tf, [], 1), [], 2);
        
        % estimate channel
        ch_est_rbs_dd = otfs_ch_est_dd_r1(rx_sym_rbs_dd, num, test_option);
        
        % equalize channel
        rx_sym_rbs_dd_eq = otfs_ch_eq_dd_r2(rx_sym_rbs_dd, ch_est_rbs_dd, num);
        
    elseif strcmp(chest_option, 'dd_tone') && ~strcmp(cheq_option, 'ddeq')
        
        % 2d inverse sfft for channel estimation only
        rx_sym_rbs_dd = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(rx_sym_rbs_tf, [], 1), [], 2);
        
        % estimate channel
        ch_est_rbs_dd = otfs_ch_est_dd_r1(rx_sym_rbs_dd, num, test_option);
%         figure, mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, abs(ch_est_rbs_dd))
%         pause
        
        % 2d inverse sfft for channel transformation
        ch_est_rbs_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(ch_est_rbs_dd, [], 2), [], 1);
        
        % equalize channel
        rx_sym_rbs_tf_eq = otfs_ch_eq_tf_r2(rx_sym_rbs_tf, ch_est_rbs_tf, noise_var, cheq_option);
        
        % 2d inverse sfft for demodulation
        rx_sym_rbs_dd_eq = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(rx_sym_rbs_tf_eq, [], 1), [], 2);
        
    elseif ~strcmp(chest_option, 'dd_tone') && strcmp(cheq_option, 'ddeq')
        
        % estimate channel
        ch_est_rbs_tf = otfs_ch_est_tf_r1(tx_sym_rbs_tf, tx_ofdmsym_faded, num, chest_option, fading_ch, ch_path_gain);
        
        % 2d inverse sfft for channel transformation
        ch_est_rbs_dd = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(ch_est_rbs_tf, [], 1), [], 2);
        
        % 2d inverse sfft for demodulation
        rx_sym_rbs_dd = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(rx_sym_rbs_tf, [], 1), [], 2);
        
        % equalize channel
        rx_sym_rbs_dd_eq = otfs_ch_eq_dd_r2(rx_sym_rbs_dd, ch_est_rbs_dd, num);
        
    else
        
        % estimate channel
        ch_est_rbs_tf = otfs_ch_est_tf_r1(tx_sym_rbs_tf, tx_ofdmsym_faded, num, chest_option, fading_ch, ch_path_gain);
        
        % equalize channel
        rx_sym_rbs_tf_eq = otfs_ch_eq_tf_r2(rx_sym_rbs_tf, ch_est_rbs_tf, noise_var, cheq_option);
        
        % 2d inverse sfft for demodulation
        rx_sym_rbs_dd_eq = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(rx_sym_rbs_tf_eq, [], 1), [], 2);
        
    end
    
    % test: calculate channel estimation error
    if test_option.ch_mse
        ch_real_rbs_tf = otfs_ch_est_tf_r1(tx_sym_rbs_tf, tx_ofdmsym_faded, num, 'real', fading_ch, ch_path_gain);
        
        % 2d inverse sfft for channel transformation
        if strcmp(chest_option, 'dd_tone') && strcmp(cheq_option, 'ddeq')
            ch_est_rbs_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(ch_est_rbs_dd, [], 2), [], 1);
        end
        
        ch_mse_usrfrm(1, idx_usrfrm) = mean(abs(ch_real_rbs_tf-ch_est_rbs_tf).^2, 'all');
        
        % test: dd-domain channel mse
        ch_real_rbs_dd = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(ch_real_rbs_tf, [], 1), [], 2);
        if ~strcmp(chest_option, 'dd_tone') && ~strcmp(cheq_option, 'ddeq')
            ch_est_rbs_dd = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(ch_est_rbs_tf, [], 1), [], 2);
        end
        fprintf('tf-domain: %6.3f    dd-domain: %6.3f\n', mean(abs(ch_real_rbs_tf-ch_est_rbs_tf).^2, 'all'), mean(abs(ch_real_rbs_dd-ch_est_rbs_dd).^2, 'all'))
%         figure
%         subplot(2, 2, 1), mesh(1:size(ch_real_rbs_tf, 2), 1:size(ch_real_rbs_tf, 1), real(ch_real_rbs_tf))
%         subplot(2, 2, 3), mesh(1:size(ch_real_rbs_tf, 2), 1:size(ch_real_rbs_tf, 1), imag(ch_real_rbs_tf))
%         subplot(2, 2, 2), mesh(1:size(ch_est_rbs_tf, 2), 1:size(ch_est_rbs_tf, 1), real(ch_est_rbs_tf))
%         subplot(2, 2, 4), mesh(1:size(ch_est_rbs_tf, 2), 1:size(ch_est_rbs_tf, 1), imag(ch_est_rbs_tf))
%         figure
%         subplot(2, 2, 1), mesh(1:size(ch_real_rbs_dd, 2), 1:size(ch_real_rbs_dd, 1), real(fftshift(fftshift(ch_real_rbs_dd, 1), 2)))
%         subplot(2, 2, 3), mesh(1:size(ch_real_rbs_dd, 2), 1:size(ch_real_rbs_dd, 1), imag(fftshift(fftshift(ch_real_rbs_dd, 1), 2)))
%         subplot(2, 2, 2), mesh(1:size(ch_est_rbs_dd, 2), 1:size(ch_est_rbs_dd, 1), real(fftshift(fftshift(ch_est_rbs_dd, 1), 2)))
%         subplot(2, 2, 4), mesh(1:size(ch_est_rbs_dd, 2), 1:size(ch_est_rbs_dd, 1), imag(fftshift(fftshift(ch_est_rbs_dd, 1), 2)))
%         figure
%         subplot(2, 3, 1), mesh(1:size(tx_sym_rbs_dd, 2), 1:size(tx_sym_rbs_dd, 1), real(tx_sym_rbs_dd))
%         subplot(2, 3, 4), mesh(1:size(tx_sym_rbs_dd, 2), 1:size(tx_sym_rbs_dd, 1), imag(tx_sym_rbs_dd))
%         subplot(2, 3, 2), mesh(1:size(rx_sym_rbs_dd, 2), 1:size(rx_sym_rbs_dd, 1), real(rx_sym_rbs_dd))
%         subplot(2, 3, 5), mesh(1:size(rx_sym_rbs_dd, 2), 1:size(rx_sym_rbs_dd, 1), imag(rx_sym_rbs_dd))
%         subplot(2, 3, 3), mesh(1:size(rx_sym_rbs_dd_eq, 2), 1:size(rx_sym_rbs_dd_eq, 1), real(rx_sym_rbs_dd_eq))
%         subplot(2, 3, 6), mesh(1:size(rx_sym_rbs_dd_eq, 2), 1:size(rx_sym_rbs_dd_eq, 1), imag(rx_sym_rbs_dd_eq))
%         
%         assignin('base', 'ch_real_rbs_tf', ch_real_rbs_tf)
%         assignin('base', 'ch_real_rbs_dd', ch_real_rbs_dd)
%         assignin('base', 'ch_est_rbs_tf', ch_est_rbs_tf)
%         assignin('base', 'ch_est_rbs_dd', ch_est_rbs_dd)
%         assignin('base', 'rx_sym_rbs_tf', rx_sym_rbs_tf)
%         assignin('base', 'rx_sym_rbs_dd', rx_sym_rbs_dd)
%         pause
        
    end
    
    % demap data qam symbols
    if strcmp(chest_option, 'dd_tone')
%         [rx_sym_data_usrfrm, rx_sym_pilot_usrfrm] = otfs_sym_demap_r2(rx_sym_rbs_dd_eq, num, test_option);      % for pilot noise estimation (now working)
        [rx_sym_data_usrfrm, ~] = otfs_sym_demap_r2(rx_sym_rbs_dd_eq, num, test_option);
    else
        rx_sym_data_usrfrm = reshape(rx_sym_rbs_dd_eq, num.num_delay_usr, []);
    end
    
    % buffer qam symbols
    rx_sym(:, idx_usrfrm) = rx_sym_data_usrfrm(:);
    
%     % buffer pilot errors (cannot be used since pilot syms have less error than data syms)
%     tx_sym_pilot_usrfrm_guard_removed = tx_sym_pilot_usrfrm(num.num_delay_guard_usr+1:end-num.num_delay_guard_usr, num.num_doppler_guard_usr+1:end-num.num_doppler_guard_usr);
%     rx_sym_pilot_usrfrm_guard_removed = rx_sym_pilot_usrfrm(num.num_delay_guard_usr+1:end-num.num_delay_guard_usr, num.num_doppler_guard_usr+1:end-num.num_doppler_guard_usr);
%     sym_pilot_error(:, idx_usrfrm) = rx_sym_pilot_usrfrm_guard_removed(:)-tx_sym_pilot_usrfrm_guard_removed(:);   % for error variance calculation
    
%     fprintf('tx_data_pwr:%6.3f    rx_data_pwr:%6.3f    tx_pilot_pwr:%6.3f    rx_pilot_pwr:%6.3f\n', std(tx_sym(:)), std(rx_sym(:)), std(tx_sym_pilot_usrfrm_guard_removed(:)), std(rx_sym_pilot_usrfrm_guard_removed(:)))
%     fprintf('tx_data_pwr:%6.3f    rx_data_pwr:%6.3f    tx_pilot_pwr:%6.3f    rx_pilot_pwr:%6.3f\n', mean(abs(tx_sym(:)).^2), mean(abs(rx_sym(:)).^2), mean(abs(tx_sym_pilot_usrfrm_guard_removed(:)).^2), mean(abs(rx_sym_pilot_usrfrm_guard_removed(:)).^2))
%     fprintf('data_std:%6.3f    pilot_std:%6.3f\n', mean(abs(tx_sym(:)-rx_sym(:)).^2), mean(abs(tx_sym_pilot_usrfrm_guard_removed(:)-rx_sym_pilot_usrfrm_guard_removed(:)).^2))
%     fprintf('data_std:%6.3f    pilot_std:%6.3f\n', std(tx_sym(:)-rx_sym(:)), std(tx_sym_pilot_usrfrm_guard_removed(:)-rx_sym_pilot_usrfrm_guard_removed(:)))
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
elseif strcmp(chest_option, 'dd_tone') && (test_option.otfs_map_plan == 1 || test_option.otfs_map_plan == 2)
    error_var = noise_var + 0.1;
elseif strcmp(chest_option, 'dd_tone') && (test_option.otfs_map_plan == 3)
    error_var = noise_var ^ 0.68; % 0.28
elseif strcmp(chest_option, 'dd_tone') && (test_option.otfs_map_plan == 4 || test_option.otfs_map_plan == 5)
    error_var = noise_var + 4.5;
elseif strcmp(chest_option, 'dd_tone') && (test_option.otfs_map_plan == 6)
    error_var = noise_var ^ 0.70; % 0.26
else
    error_var = noise_var;
end

% test channel estimation error variance
if test_option.sym_err_var
    sym_err_var = var(tx_sym_serial-rx_sym_serial);
%     fprintf('noise_var:%6.3f    sym_err_var:%6.3f    calc_err_var:%6.3f\n', noise_var, sym_err_var, error_var)
%     figure
%     subplot(2, 1, 1), plot(real(tx_sym_serial), '-b.'), hold on, plot(real(rx_sym_serial), '-r'), hold off, grid minor
%     subplot(2, 1, 2), plot(imag(tx_sym_serial), '-b.'), hold on, plot(imag(rx_sym_serial), '-r'), hold off, grid minor
%     pause
else
    sym_err_var = 0;
end

% demap the qam symbols
% noise_var_tmp = var(tx_sym_serial-rx_sym_serial);
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

% remove padded bits
rx_bit = rx_bit_crc_removed(1 : sim.len_tb_bit);

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

% assignin('base', 'rx_bit_demod', rx_bit_demod);
% assignin('base', 'rx_bit_ratematch', rx_bit_ratematch);
% assignin('base', 'rx_bit_dec', rx_bit_dec);
% assignin('base', 'rx_bit', rx_bit);
% assignin('base', 'tx_bit', tx_bit);
% figure(1), plot(rx_bit_demod, '-b.'), grid minor
% figure(2), plot(xor(tx_bit, rx_bit), '-b.'), grid minor
% var(tx_sym(:)-rx_sym(:))
% pause

% % dump
% assignin('base', 'tx_sym_rbs_dd', tx_sym_rbs_dd);
% assignin('base', 'tx_sym_rbs_tf', tx_sym_rbs_tf);
% assignin('base', 'tx_sym_bw', tx_sym_bw);
% assignin('base', 'tx_sym_nfft_shift', tx_sym_nfft_shift);
% assignin('base', 'tx_sym_nfft', tx_sym_nfft);
% assignin('base', 'tx_ofdmsym', tx_ofdmsym);
% assignin('base', 'tx_ofdmsym_cp', tx_ofdmsym_cp);
% assignin('base', 'tx_ofdmsym_serial', tx_ofdmsym_serial);
% assignin('base', 'tx_ofdmsym_faded', tx_ofdmsym_faded);
% assignin('base', 'rx_ofdmsym_serial', rx_ofdmsym_serial);
% assignin('base', 'rx_ofdmsym_cp', rx_ofdmsym_cp);
% assignin('base', 'rx_ofdmsym', rx_ofdmsym);
% assignin('base', 'rx_sym_nfft', rx_sym_nfft);
% assignin('base', 'rx_sym_nfft_shift', rx_sym_nfft_shift);
% assignin('base', 'rx_sym_bw', rx_sym_bw);
% assignin('base', 'rx_sym_rbs_dd', rx_sym_rbs_dd);
% assignin('base', 'rx_sym_rbs_tf', rx_sym_rbs_tf);
% assignin('base', 'rx_sym_rbs_tf_eq', rx_sym_rbs_tf_eq);
% assignin('base', 'rx_sym_rbs_dd_eq', rx_sym_rbs_dd_eq);
% pause

end
