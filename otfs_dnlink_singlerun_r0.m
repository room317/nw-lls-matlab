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

function [pkt_error, tx_papr, ch_mse, sym_err_var] = otfs_dnlink_singlerun_r0(sim, cc, rm, num, snr_db, ch, chest_option, cheq_option, test_option)

% create turbo encoder/decoder
turbo_enc = comm.TurboEncoder('TrellisStructure', cc.tc_trellis, ...
    'InterleaverIndices', cc.PI+1);
turbo_dec = comm.TurboDecoder('TrellisStructure', cc.tc_trellis, ...
    'InterleaverIndices', cc.PI+1, 'NumIterations', cc.num_iter_max);

%% common objects and parameters

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
else
    fading_ch = comm.RayleighChannel(...
        'SampleRate', num.sample_rate, ...
        'PathDelays', ch.path_delays, ...
        'AveragePathGains', ch.average_path_gains, ...
        'NormalizePathGains', true, ...
        'MaximumDopplerShift', ch.maximum_doppler_shift, ...
        'DopplerSpectrum', ch.doppler_spectrum,...
        'PathGainsOutputPort', true);
end

% calculate noise variance
noise_var = (num.num_subc_usr/num.num_fft)*(10^((-0.1)*snr_db));

%% base station (tx) operation per user (per user data generation)

% simulate per-user tx process
tx_bit = zeros(sim.len_tb_bit, num.num_usr);
tx_sym_serial = zeros(num.num_qamsym_usr*rm.num_usrfrm, num.num_usr);
tx_sym = zeros(num.num_qamsym_usr, rm.num_usrfrm, num.num_usr);
tx_sym_rbs_dd = zeros(num.num_delay_usr, num.num_doppler_usr, rm.num_usrfrm, num.num_usr);
tx_sym_rbs_tf = zeros(num.num_subc_usr, num.num_ofdmsym_usr, rm.num_usrfrm, num.num_usr);
tx_sym_bw = zeros(num.num_subc_bw, num.num_ofdmsym, rm.num_usrfrm);
for idx_usr = 1:num.num_usr
    
    % generate ofdm/otfs qam symbols per user
    [tx_bit_usr, tx_sym_serial_usr, tx_sym_usr] = gen_tx_qamsym(sim, cc, rm, num, turbo_enc);
    tx_bit(:, idx_usr) = tx_bit_usr;
    tx_sym_serial(:, idx_usr) = tx_sym_serial_usr;
    tx_sym(:, :, idx_usr) = tx_sym_usr;
    
    % simulate per-user-frame tx process
    for idx_usrfrm = 1:rm.num_usrfrm
        
        % extract user frame data per user
        tx_sym_data_usrfrm = squeeze(tx_sym(:, idx_usrfrm, idx_usr));
        
        % map data and pilot symbols
        if strcmp(chest_option, 'dd_tone')
            [tx_sym_rbs_usr_dd, ~] = otfs_sym_map_r2(tx_sym_data_usrfrm, num, test_option);
        else
            tx_sym_rbs_usr_dd = reshape(tx_sym_data_usrfrm, num.num_delay_usr, []);
        end
        tx_sym_rbs_dd(:, :, idx_usrfrm, idx_usr) = tx_sym_rbs_usr_dd;
        
        % 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
        tx_sym_rbs_usr_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(tx_sym_rbs_usr_dd, [], 2), [], 1);
        tx_sym_rbs_tf(:, :, idx_usrfrm, idx_usr) = tx_sym_rbs_usr_tf;
        
        % map user resource blocks to whole bandwidth
        usr_id = num.list_usr(idx_usr);                                         % user id
        idx_rb_usr = (ceil(usr_id/num.max_usr_slot)-1)*num.num_rb_usr+1;        % user rb index
        idx_slot_usr = mod(usr_id-1, num.max_usr_slot)*num.num_slot_usr+1;      % user slot index
        list_subc_usr = (idx_rb_usr-1)*num.num_subc_rb+1:(idx_rb_usr-1)*num.num_subc_rb+num.num_subc_usr;                       % user subcarrier list
        list_ofdmsym_usr = (idx_slot_usr-1)*num.num_ofdmsym_slot+1:(idx_slot_usr-1)*num.num_ofdmsym_slot+num.num_ofdmsym_usr;   % user ofdm symbol list
        tx_sym_bw(list_subc_usr, list_ofdmsym_usr, idx_usrfrm) = tx_sym_rbs_usr_tf;
    end
end

%% base station (tx) (serial downlink signal generation)

% proceed per-user-frame (common for all users)
tx_ofdmsym_serial = zeros((num.num_fft+num.num_cp)*num.num_ofdmsym, rm.num_usrfrm);
for idx_usrfrm = 1:rm.num_usrfrm
    
    % map bandwidth to the center of fft range
    tx_sym_nfft_shift = zeros(num.num_fft, num.num_ofdmsym);
    tx_sym_nfft_shift((num.num_fft/2)-(num.num_subc_bw/2)+1:(num.num_fft/2)+(num.num_subc_bw/2), :) = squeeze(tx_sym_bw(:, :, idx_usrfrm));
    tx_sym_nfft = fftshift(tx_sym_nfft_shift, 1);
    
    % ofdm modulate
    tx_ofdmsym = sqrt(num.num_fft)*ifft(tx_sym_nfft, [], 1);
    
    % add cp
    tx_ofdmsym_cp = tx_ofdmsym([num.num_fft-num.num_cp+1:num.num_fft 1:num.num_fft], :);
    
    % serialize
    tx_ofdmsym_serial(:, idx_usrfrm) = tx_ofdmsym_cp(:);
end

%% channel per user

% initialize
tx_ofdmsym_faded = zeros((num.num_fft+num.num_cp)*num.num_ofdmsym, rm.num_usrfrm, num.num_usr);
rx_ofdmsym_serial = zeros((num.num_fft+num.num_cp)*num.num_ofdmsym, rm.num_usrfrm, num.num_usr);
if test_option.ch_mse || strcmp(test_chest, 'real')
    ch_real_mat_tf = zeros(num.num_subc_bw, num.num_subc_bw, num.num_ofdmsym, rm.num_usrfrm, num.num_usr);
else
    ch_real_mat_tf = [];
end

% simulate per-user and per-user-frame
for idx_usr = 1:num.num_usr
    for idx_usrfrm = 1:rm.num_usrfrm
        
        % pass signal through channel
        [tx_ofdmsym_usr_faded, ch_path_gain_usr] = fading_ch(squeeze(tx_ofdmsym_serial(:, idx_usrfrm)));
        tx_ofdmsym_faded(:, idx_usrfrm, idx_usr) = tx_ofdmsym_usr_faded;
        
        % add gaussian noise
        rx_ofdmsym_serial(:, idx_usrfrm, idx_usr) = awgn(tx_ofdmsym_usr_faded, snr_db, 'measured');
        
        % regenerate real channel
        if test_option.ch_mse || strcmp(test_chest, 'real')
            [ch_real_mat_usr_tf, ~, ~, ~] = gen_real_ch_r1(fading_ch, ch_path_gain_usr, num.num_fft, num.num_cp, num.num_subc_bw, num.num_ofdmsym, [], [], false);
            ch_real_mat_tf(:, :, :, idx_usrfrm, idx_usr) = ch_real_mat_usr_tf;
        end
        
%         % test (demap test)
%         tx_ofdmsym_faded(:, idx_usrfrm, idx_usr) = tx_ofdmsym_serial(:, idx_usrfrm);
%         rx_ofdmsym_serial(:, idx_usrfrm, idx_usr) = tx_ofdmsym_serial(:, idx_usrfrm);
%         ch_real_mat_tf(:, :, :, idx_usrfrm, idx_usr) = repmat(eye(num.num_subc_bw), 1, 1, num.num_ofdmsym);
    end
    
    % copy channels when common channel is used for all users
    if test_option.common_usr_ch
        tx_ofdmsym_faded = repmat(tx_ofdmsym_faded(:, :, 1), 1, 1, num.num_usr);
        rx_ofdmsym_serial = repmat(rx_ofdmsym_serial(:, :, 1), 1, 1, num.num_usr);
        ch_real_mat_tf = repmat(ch_real_mat_tf(:, :, :, :, 1), 1, 1, 1, 1, num.num_usr);
        break
    end
end

%% mobile station (rx) operation per user

rx_bit = zeros(sim.len_tb_bit, num.num_usr);
rx_sym_serial = zeros(num.num_qamsym_usr*rm.num_usrfrm, num.num_usr);
rx_sym = zeros(num.num_qamsym_usr, rm.num_usrfrm, num.num_usr);
rx_sym_rbs_dd_eq = zeros(num.num_delay_usr, num.num_doppler_usr, rm.num_usrfrm, num.num_usr);
rx_sym_rbs_tf = zeros(num.num_subc_usr, num.num_ofdmsym_usr, rm.num_usrfrm, num.num_usr);
rx_crc_error = zeros(1, num.num_usr);
ch_real_eff_tf = zeros(num.num_subc_usr*num.num_ofdmsym_usr, num.num_subc_usr*num.num_ofdmsym_usr, rm.num_usrfrm, num.num_usr);
ch_est_rbs_tf = zeros(num.num_subc_usr, num.num_ofdmsym_usr, rm.num_usrfrm, num.num_usr);
ch_est_rbs_dd = zeros(num.num_delay_usr, num.num_doppler_usr, rm.num_usrfrm, num.num_usr);
error_var = zeros(1, num.num_usr);

for idx_usr = 1:num.num_usr
    for idx_usrfrm = 1:rm.num_usrfrm
        
        % reshape
        rx_ofdmsym_cp = reshape(squeeze(rx_ofdmsym_serial(:, idx_usrfrm, idx_usr)), num.num_fft+num.num_cp, num.num_ofdmsym);
        
        % remove cp
        rx_ofdmsym = rx_ofdmsym_cp(num.num_cp+1:end,:);
        
        % remove parts of received symbols
        if isscalar(test_option.partial_reception) && ismember(test_option.partial_reception, 1:num.num_ofdmsym-1)
            rx_ofdmsym_part = [rx_ofdmsym(:, 1:test_option.partial_reception) zeros(size(rx_ofdmsym, 1), num.num_ofdmsym-test_option.partial_reception)];
        else
            rx_ofdmsym_part = rx_ofdmsym;
        end
        
        % ofdm demodulate
        rx_sym_nfft = (1/sqrt(num.num_fft))*fft(rx_ofdmsym_part, [], 1);
        
        % demap symbols in bandwidth from fft range
        rx_sym_nfft_shift = fftshift(rx_sym_nfft, 1);
        rx_sym_bw = rx_sym_nfft_shift((num.num_fft/2)-(num.num_subc_bw/2)+1:(num.num_fft/2)+(num.num_subc_bw/2), :);
        
        % demap user resource block from whole bandwidth (temporary)
        usr_id = num.list_usr(idx_usr);                                         % user id
        idx_rb_usr = (ceil(usr_id/num.max_usr_slot)-1)*num.num_rb_usr+1;        % user rb index
        idx_slot_usr = mod(usr_id-1, num.max_usr_slot)*num.num_slot_usr+1;      % user slot index
        list_subc_usr = (idx_rb_usr-1)*num.num_subc_rb+1:(idx_rb_usr-1)*num.num_subc_rb+num.num_subc_usr;                       % user subcarrier list
        list_ofdmsym_usr = (idx_slot_usr-1)*num.num_ofdmsym_slot+1:(idx_slot_usr-1)*num.num_ofdmsym_slot+num.num_ofdmsym_usr;   % user ofdm symbol list
        rx_sym_rbs_usr_tf = rx_sym_bw(list_subc_usr, list_ofdmsym_usr);
        rx_sym_rbs_tf(:, :, idx_usrfrm, idx_usr) = rx_sym_rbs_usr_tf;
        
        % demap user real channel (for real channel estimation or channel estimation error check)
        [ch_real_rbs_usr_tf, ch_real_eff_usr_tf, ch_real_eff_usr_dd] = ...
            demap_real_ch(squeeze(ch_real_mat_tf(:, :, :, idx_usrfrm, idx_usr)), list_subc_usr, list_ofdmsym_usr, chest_option, cheq_option, test_option);
        
        % test: to check channel estimation error
        if test_option.ch_mse
            ch_real_eff_tf(:, :, idx_usrfrm, idx_usr) = ch_real_eff_usr_tf;
        end
        
        % estimate channel
        if strcmp(chest_option, 'real')
            % demap user channel
            ch_est_rbs_usr_tf = ch_real_rbs_usr_tf;
            ch_est_rbs_usr_dd = [];
        
            % test: to check channel estimation error
            if test_option.ch_mse && ~test_option.fulltap_eq
                ch_est_rbs_tf(:, :, idx_usrfrm, idx_usr) = ch_est_rbs_usr_tf;
            end
        elseif strcmp(chest_option, 'dd_tone')
            % 2d inverse sfft for channel estimation only
            rx_sym_rbs_usr_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(rx_sym_rbs_usr_tf, [], 1), [], 2);
            
            % estimate channel
            ch_est_rbs_usr_dd = otfs_ch_est_dd_r1(rx_sym_rbs_usr_dd, num, test_option);
            ch_est_rbs_usr_tf = [];
            
            % test: to check channel estimation error
            if test_option.ch_mse
                ch_est_rbs_dd(:, :, idx_usrfrm, idx_usr) = ch_est_rbs_usr_dd;
            end
        else
            % estimate channel
            ch_est_rbs_usr_tf = otfs_ch_est_tf_r2(squeeze(tx_sym_rbs_tf(:, :, idx_usrfrm, idx_usr)), squeeze(tx_ofdmsym_faded(:, idx_usrfrm, idx_usr)), list_subc_usr, list_ofdmsym_usr, num, chest_option);
            ch_est_rbs_usr_dd = [];
            
            % test: to check channel estimation error
            if test_option.ch_mse
                ch_est_rbs_tf(:, :, idx_usrfrm, idx_usr) = ch_est_rbs_usr_tf;
            end
        end
        
        % equalize channel
        if strcmp(cheq_option, 'ddeq_zf') || strcmp(cheq_option, 'ddeq_mmse')
            % check delay-doppler channel
            if isempty(ch_est_rbs_usr_dd)
                % 2d inverse sfft both for demodulation
                rx_sym_rbs_usr_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(rx_sym_rbs_usr_tf, [], 1), [], 2);
            end
            
            % equalize channel
            rx_sym_rbs_usr_dd_eq = otfs_ch_eq_dd_r3(rx_sym_rbs_usr_dd, ch_est_rbs_usr_tf, ch_est_rbs_usr_dd, ch_real_eff_usr_dd, num, noise_var, chest_option, cheq_option, test_option);
        elseif strcmp(cheq_option, 'tfeq_zf') || strcmp(cheq_option, 'tfeq_mmse')
            % equalize channel
            rx_sym_rbs_usr_tf_eq = otfs_ch_eq_tf_r3(rx_sym_rbs_usr_tf, ch_est_rbs_usr_tf, ch_est_rbs_usr_dd, ch_real_eff_usr_tf, num, noise_var, chest_option, cheq_option, test_option);
            
            % 2d inverse sfft for demodulation
            rx_sym_rbs_usr_dd_eq = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(rx_sym_rbs_usr_tf_eq, [], 1), [], 2);
        else
            error('cheq_option value must be one of these: {tfeq_zf, tfeq_mmse, ddeq_zf, ddeq_mmse}')
        end
        rx_sym_rbs_dd_eq(:, :, idx_usrfrm, idx_usr) = rx_sym_rbs_usr_dd_eq;
        
        % demap data qam symbols
        if strcmp(chest_option, 'dd_tone')
            [rx_sym_data_usrfrm, ~] = otfs_sym_demap_r2(rx_sym_rbs_usr_dd_eq, num, test_option);
        else
            rx_sym_data_usrfrm = rx_sym_rbs_usr_dd_eq(:);
        end
        
        % buffer qam symbols
        rx_sym(:, idx_usrfrm, idx_usr) = rx_sym_data_usrfrm;
        
    end
    
    % compensate channel estimation error variance
    % qam_mse = mean(abs(tx_sym(:)-rx_sym(:)).^2);
    % error_var = var(sym_pilot_error(:));    % cannot be used (pilot syms have less error than data syms)
    % error_var = var(tx_sym_serial-rx_sym_serial);
    if strcmp(chest_option, 'tf_ltedown') || strcmp(chest_option, 'tf_nr')
        error_var_usr = noise_var + 0.12;
    elseif strcmp(chest_option, 'tf_lteup')
        error_var_usr = noise_var + 0.27;
    elseif strcmp(chest_option, 'real')
        error_var_usr = noise_var + 0.08;
    elseif strcmp(chest_option, 'dd_tone') && (test_option.otfs_map_plan == 1 || test_option.otfs_map_plan == 2)
        error_var_usr = noise_var + 0.1;
    elseif strcmp(chest_option, 'dd_tone') && (test_option.otfs_map_plan == 3)
        error_var_usr = noise_var ^ 0.68; % 0.28
    elseif strcmp(chest_option, 'dd_tone') && (test_option.otfs_map_plan == 4 || test_option.otfs_map_plan == 5)
        error_var_usr = noise_var + 4.5;
    elseif strcmp(chest_option, 'dd_tone') && (test_option.otfs_map_plan == 6)
        error_var_usr = noise_var ^ 0.70; % 0.26
    else
        error_var_usr = noise_var;
    end
        
    % test: to check channel estimation error
    if test_option.ch_mse
        error_var(1, idx_usr) = error_var_usr;
    end

    % generate ofdm/otfs received bits per user
    [rx_bit_usr, rx_sym_serial_usr, rx_crc_error_usr] = gen_rx_bit(rx_sym(:, :, idx_usr), turbo_dec, error_var_usr, sim, cc, rm);
    rx_bit(:, idx_usr) = rx_bit_usr;
    rx_sym_serial(:, idx_usr) = rx_sym_serial_usr;
    rx_crc_error(:, idx_usr) = rx_crc_error_usr;
    
end

% calculate packet error
pkt_error = rx_crc_error | (symerr(tx_bit, rx_bit, 'column-wise') > 0);

% test: calculate papr
if test_option.papr
    tx_papr = mean(max(abs(tx_ofdmsym_serial).^2, [], 1)./mean(abs(tx_ofdmsym_serial).^2, 1));
else
    tx_papr = [];
end

% test: calculate channel mse
if test_option.ch_mse
    if strcmp(chest_option, 'real') && test_option.fulltap_eq
        ch_mse = zeros(1, num.num_usr);
    elseif strcmp(chest_option, 'dd_tone')
        ch_est_rbs_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(ch_est_rbs_dd, [], 2), [], 1);
        
        % generate effective tf channel
        ch_est_eff_tf = zeros(num.num_subc_usr*num.num_ofdmsym_usr, num.num_subc_usr*num.num_ofdmsym_usr, rm.num_usrfrm, num.num_usr);
        for idx_usr = 1:num.num_usr
            for idx_usrfrm = 1:rm.num_usrfrm
                ch_est_eff_tf(:, :, idx_usrfrm, idx_usr) = diag(reshape(ch_est_rbs_tf(:, :, idx_usrfrm, idx_usr), [], 1));
            end
        end
        
        % calculate channel estimation error (per time-frequency grid)
        ch_mse = mean(abs(ch_real_eff_tf-ch_est_eff_tf).^2, [1 2 3])*num.num_subc_usr*num.num_ofdmsym_usr;
        
%         fprintf('ch mse:')
%         fprintf(' %10.6f', ch_mse)
%         fprintf('\n')
%         assignin('base', 'ch_real_eff_tf', ch_real_eff_tf)
%         assignin('base', 'ch_est_eff_tf', ch_est_eff_tf)
    else
        % generate effective tf channel (per time-frequency grid)
        ch_est_eff_tf = zeros(num.num_subc_usr*num.num_ofdmsym_usr, num.num_subc_usr*num.num_ofdmsym_usr, rm.num_usrfrm, num.num_usr);
        for idx_usr = 1:num.num_usr
            for idx_usrfrm = 1:rm.num_usrfrm
                ch_est_eff_tf(:, :, idx_usrfrm, idx_usr) = diag(reshape(ch_est_rbs_tf(:, :, idx_usrfrm, idx_usr), [], 1));
            end
        end
        
        % calculate channel estimation error
        ch_mse = mean(abs(ch_real_eff_tf-ch_est_eff_tf).^2, [1 2 3])*num.num_subc_usr*num.num_ofdmsym_usr;
        
%         fprintf('ch mse:')
%         fprintf(' %10.6f', ch_mse)
%         fprintf('\n')
%         assignin('base', 'ch_real_eff_tf', ch_real_eff_tf)
%         assignin('base', 'ch_est_eff_tf', ch_est_eff_tf)
    end
%     for idx_usr=1:num.num_usr
%         fprintf('user %d:    tf channel est rmse: %6.3f\n', idx_usr, sqrt(ch_mse(1, idx_usr)))
%         for idx_usrfrm = 1:rm.num_usrfrm
%             test_ch_real_onetap_tf = reshape(diag(ch_real_eff_tf(:, :, idx_usrfrm, idx_usr)), num.num_subc_usr, num.num_ofdmsym_usr);
%             test_ch_real_onetap_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(test_ch_real_onetap_tf, [], 1), [], 2);
%             if strcmp(chest_option, 'real') && test_option.fulltap_eq
%                 test_ch_est_onetap_tf = test_ch_real_onetap_tf;
%                 test_ch_est_onetap_dd = test_ch_real_onetap_dd;
%             else
%                 test_ch_est_onetap_tf = ch_est_rbs_tf(:, :, idx_usrfrm, idx_usr);
%                 test_ch_est_onetap_dd = ch_est_rbs_dd(:, :, idx_usrfrm, idx_usr);
%             end
%             test_tx_sym_dd = tx_sym_rbs_dd(:, :, idx_usrfrm, idx_usr);
%             test_rx_sym_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(rx_sym_rbs_tf(:, :, idx_usrfrm, idx_usr), [], 1), [], 2);
%             test_rx_sym_eq_dd = rx_sym_rbs_dd_eq(:, :, idx_usrfrm, idx_usr);
%             
%             figure
%             subplot(2, 2, 1), mesh(1:size(test_ch_real_onetap_tf, 2), 1:size(test_ch_real_onetap_tf, 1), real(test_ch_real_onetap_tf))
%             subplot(2, 2, 3), mesh(1:size(test_ch_real_onetap_tf, 2), 1:size(test_ch_real_onetap_tf, 1), imag(test_ch_real_onetap_tf))
%             subplot(2, 2, 2), mesh(1:size(test_ch_est_onetap_tf, 2), 1:size(test_ch_est_onetap_tf, 1), real(test_ch_est_onetap_tf))
%             subplot(2, 2, 4), mesh(1:size(test_ch_est_onetap_tf, 2), 1:size(test_ch_est_onetap_tf, 1), imag(test_ch_est_onetap_tf))
%             
%             figure
%             subplot(2, 2, 1), mesh(1:size(test_ch_real_onetap_dd, 2), 1:size(test_ch_real_onetap_dd, 1), real(fftshift(fftshift(test_ch_real_onetap_dd, 1), 2)))
%             subplot(2, 2, 3), mesh(1:size(test_ch_real_onetap_dd, 2), 1:size(test_ch_real_onetap_dd, 1), imag(fftshift(fftshift(test_ch_real_onetap_dd, 1), 2)))
%             subplot(2, 2, 2), mesh(1:size(test_ch_est_onetap_dd, 2), 1:size(test_ch_est_onetap_dd, 1), real(fftshift(fftshift(test_ch_est_onetap_dd, 1), 2)))
%             subplot(2, 2, 4), mesh(1:size(test_ch_est_onetap_dd, 2), 1:size(test_ch_est_onetap_dd, 1), imag(fftshift(fftshift(test_ch_est_onetap_dd, 1), 2)))
%             
%             figure
%             subplot(2, 3, 1), mesh(1:size(test_tx_sym_dd, 2), 1:size(test_tx_sym_dd, 1), real(test_tx_sym_dd))
%             subplot(2, 3, 4), mesh(1:size(test_tx_sym_dd, 2), 1:size(test_tx_sym_dd, 1), imag(test_tx_sym_dd))
%             subplot(2, 3, 2), mesh(1:size(test_rx_sym_dd, 2), 1:size(test_rx_sym_dd, 1), real(test_rx_sym_dd))
%             subplot(2, 3, 5), mesh(1:size(test_rx_sym_dd, 2), 1:size(test_rx_sym_dd, 1), imag(test_rx_sym_dd))
%             subplot(2, 3, 3), mesh(1:size(test_rx_sym_eq_dd, 2), 1:size(test_rx_sym_eq_dd, 1), real(test_rx_sym_eq_dd))
%             subplot(2, 3, 6), mesh(1:size(test_rx_sym_eq_dd, 2), 1:size(test_rx_sym_eq_dd, 1), imag(test_rx_sym_eq_dd))
%             pause
%         end
%     end
else
    ch_mse = [];
end

% test: channel estimation error variance
if test_option.sym_err_var
    sym_err_var = var(tx_sym_serial-rx_sym_serial, 0, 1);
    
%     for idx_usr = 1:num.num_usr
%         fprintf('user %d:    noise_var:%6.3f    sym_err_var:%6.3f    calc_err_var:%6.3f\n', idx_usr, noise_var, sym_err_var(1, idx_usr), error_var(1, idx_usr))
%         
%         figure
%         subplot(2, 1, 1), plot(real(tx_sym_serial(:, idx_usr)), '-b.'), hold on, plot(real(rx_sym_serial(:, idx_usr)), '-r'), hold off, grid minor
%         subplot(2, 1, 2), plot(imag(tx_sym_serial(:, idx_usr)), '-b.'), hold on, plot(imag(rx_sym_serial(:, idx_usr)), '-r'), hold off, grid minor
%         pause
%     end
else
    sym_err_var = [];
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
