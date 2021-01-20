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

function [pkt_error, tx_papr, ch_mse, sym_err_var] = ofdm_dnlink_singlerun_r1(sim, cc, rm, num, snr_db, ch, chest_option, cheq_option, test_option)

%% common objects and parameters

% create crc generator and detector objects
tx_crc = comm.CRCGenerator(cc.gCRC24A);
rx_crc = comm.CRCDetector(cc.gCRC24A);

% create turbo encoder/decoder objects
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
else
    fading_ch = comm.RayleighChannel(...
        'SampleRate', num.sample_rate, ...
        'PathDelays', ch.path_delays, ...
        'AveragePathGains', ch.average_path_gains, ...
        'NormalizePathGains', true, ...
        'MaximumDopplerShift', ch.maximum_doppler_shift, ...
        'DopplerSpectrum', ch.doppler_spectrum, ...
        'PathGainsOutputPort', true);
end

% calculate snr and noise variance
%   - snr_db(input) is consigered to be the snr(db) for user data in delay-doppler domain
%   - snr_db_adj is consigered to be the snr(db) in time domain wave in the air
%   - noise_var is noise variance for user data in delay-doppler domain
snr_db_adj = snr_db+(10*log10((num.num_usr*num.num_rb_usr*num.num_slot_usr)/(num.num_rb*num.num_slot)))+(10*log10(num.num_subc_bw/num.num_fft));
noise_var = 10^((-0.1)*snr_db);

% % calculate snr and noise variance (primary)
% snr_db_adj = snr_db-(10*log10((num.num_rb*num.num_slot)/(num.num_usr*num.num_rb_usr*num.num_slot_usr)));    % adjust time-domain snr wrt resource size
% noise_var = 10^((-0.1)*snr_db)*(num.num_subc_bw/num.num_fft);                                               % calculate noise variance at dd resource block

% % calculate snr and noise variance (temp)
% snr_db_adj = snr_db-(10*log10((num.num_fft+num.num_cp)/num.num_fft));
% % snr_db_adj = snr_db+(10*log10((num.num_usr*num.num_rb_usr*num.num_slot_usr)/(num.num_rb*num.num_slot)))+(10*log10((num.num_subc_bw/num.num_fft)))-(10*log10((num.num_fft+num.num_cp)/num.num_fft));
% % noise_var = 10^((-0.1)*snr_db_adj)*(num.num_subc_bw/num.num_fft)*((num.num_usr*num.num_rb_usr*num.num_slot_usr)/(num.num_rb*num.num_slot));
% noise_var = 10^((-0.1)*snr_db);

%% base station (tx) operation per user (per user data generation)

% initialize
tx_bit = zeros(sim.len_tb_bit, num.num_usr);
tx_sym_serial = zeros(num.num_qamsym_usr*sum(rm.num_usrfrm_cb), num.num_usr);
tx_sym = zeros(num.num_qamsym_usr, sum(rm.num_usrfrm_cb), num.num_usr);
tx_sym_rbs = zeros(num.num_subc_usr, num.num_ofdmsym_usr, sum(rm.num_usrfrm_cb), num.num_usr);
tx_sym_bw = zeros(num.num_subc_bw, num.num_ofdmsym, sum(rm.num_usrfrm_cb));

% simulate per-user tx process
for idx_usr = 1:num.num_usr
    
    % generate bit stream
    tx_bit(:, idx_usr) = randi([0 1], sim.len_tb_bit, 1);
    
    % pad bits
    tx_bit_pad_usr = [tx_bit(:, idx_usr); zeros(cc.F, 1)];
    
    % buffer per code block
    tx_bit_buff_usr = reshape(tx_bit_pad_usr, (sim.len_tb_bit+cc.F)/cc.C, cc.C);
    
    % simulate per-cb tx process
    for idx_cb = 1:cc.C
        
        % generate crc
        tx_bit_crc_cb = tx_crc(tx_bit_buff_usr(:, idx_cb));
        
        % turbo encode data
        tx_bit_enc_cb = turbo_enc(tx_bit_crc_cb);
        
        % rate match
        tx_bit_ratematch_cb = tx_ratematch_r2(tx_bit_enc_cb, idx_cb, rm);
        
        % modulate bit stream
        % accumulate output per code block serially
        tx_sym_serial_cb = qammod(tx_bit_ratematch_cb(:), 2^rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);
        tx_sym_serial(sum(rm.num_usrfrm_cb(1:idx_cb-1))*num.num_qamsym_usr+1:sum(rm.num_usrfrm_cb(1:idx_cb))*num.num_qamsym_usr, idx_usr) = tx_sym_serial_cb;
        
        % simulate per user frame
        tx_sym_cb = reshape(tx_sym_serial_cb, num.num_qamsym_usr, []);
        tx_sym(:, sum(rm.num_usrfrm_cb(1:idx_cb-1))+1:sum(rm.num_usrfrm_cb(1:idx_cb)), idx_usr) = tx_sym_cb;
    end
    
    % simulate per-user-frame tx process
    for idx_usrfrm = 1:sum(rm.num_usrfrm_cb)
        
        % extract user frame data per user
        tx_sym_data_usrfrm = squeeze(tx_sym(:, idx_usrfrm, idx_usr));
        
        % map qam symbols to user physical resource blocks
        [tx_sym_rbs_usr, ~] = ofdm_sym_map(tx_sym_data_usrfrm, num);
        tx_sym_rbs(:, :, idx_usrfrm, idx_usr) = tx_sym_rbs_usr;
        
        % map user resource blocks to whole bandwidth
        usr_id = num.list_usr(idx_usr);                                         % user id
        idx_rb_usr = (ceil(usr_id/num.max_usr_slot)-1)*num.num_rb_usr+1;        % user rb index
        idx_slot_usr = mod(usr_id-1, num.max_usr_slot)*num.num_slot_usr+1;      % user slot index
        list_subc_usr = (idx_rb_usr-1)*num.num_subc_rb+1:(idx_rb_usr-1)*num.num_subc_rb+num.num_subc_usr;                       % user subcarrier list
        list_ofdmsym_usr = (idx_slot_usr-1)*num.num_ofdmsym_slot+1:(idx_slot_usr-1)*num.num_ofdmsym_slot+num.num_ofdmsym_usr;   % user ofdm symbol list
        tx_sym_bw(list_subc_usr, list_ofdmsym_usr, idx_usrfrm) = tx_sym_rbs_usr;
    end
end

%% base station (tx) (serial downlink signal generation)

% proceed per-user-frame (common for all users)
tx_ofdmsym_serial = zeros((num.num_fft+num.num_cp)*num.num_ofdmsym, sum(rm.num_usrfrm_cb));
for idx_usrfrm = 1:sum(rm.num_usrfrm_cb)
    
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
tx_ofdmsym_faded = zeros((num.num_fft+num.num_cp)*num.num_ofdmsym, sum(rm.num_usrfrm_cb), num.num_usr);
rx_ofdmsym_serial = zeros((num.num_fft+num.num_cp)*num.num_ofdmsym, sum(rm.num_usrfrm_cb), num.num_usr);
if ~test_option.awgn && (test_option.ch_mse || strcmp(chest_option, 'real'))
    ch_real_mat_tf = zeros(num.num_subc_bw, num.num_subc_bw, num.num_ofdmsym, sum(rm.num_usrfrm_cb), num.num_usr);
else
    ch_real_mat_tf = [];
end

% simulate per-user and per-user-frame
for idx_usr = 1:num.num_usr
    for idx_usrfrm = 1:sum(rm.num_usrfrm_cb)
        
        % pass signal through channel
        if test_option.awgn
            tx_ofdmsym_usr_faded = tx_ofdmsym_serial(:, idx_usrfrm);
        else
            [tx_ofdmsym_usr_faded, ch_path_gain_usr] = fading_ch(squeeze(tx_ofdmsym_serial(:, idx_usrfrm)));
        end
        tx_ofdmsym_faded(:, idx_usrfrm, idx_usr) = tx_ofdmsym_usr_faded;
        
        % add gaussian noise
        if test_option.fading_only
            rx_ofdmsym_serial(:, idx_usrfrm, idx_usr) = tx_ofdmsym_usr_faded;
        else
            rx_ofdmsym_serial(:, idx_usrfrm, idx_usr) = awgn(tx_ofdmsym_usr_faded, snr_db_adj, 'measured');
%             rx_ofdmsym_serial(:, idx_usrfrm, idx_usr) = tx_ofdmsym_usr_faded+(sqrt(1/2)*10^(-0.1*snr_db_adj)*complex(randn(size(tx_ofdmsym_usr_faded)), randn(size(tx_ofdmsym_usr_faded))));
        end
        
        % regenerate real channel
        if ~test_option.awgn && (test_option.ch_mse || strcmp(chest_option, 'real'))
            [~, ch_real_mat_usr_tf, ~, ~, ~] = gen_real_ch_r1(fading_ch, ch_path_gain_usr, num.num_fft, num.num_cp, num.num_subc_bw, num.num_ofdmsym, [], [], false, test_option);
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
        if ~test_option.awgn
            ch_real_mat_tf = repmat(ch_real_mat_tf(:, :, :, :, 1), 1, 1, 1, 1, num.num_usr);
        end
        break
    end
end

%% mobile station (rx) operation per user

% initialize
rx_bit = zeros(sim.len_tb_bit, num.num_usr);
rx_crc_error = zeros(1, num.num_usr);
rx_sym_serial = zeros(num.num_qamsym_usr*sum(rm.num_usrfrm_cb), num.num_usr);
rx_sym = zeros(num.num_qamsym_usr, sum(rm.num_usrfrm_cb), num.num_usr);
rx_sym_rbs = zeros(num.num_subc_usr, num.num_ofdmsym_usr, sum(rm.num_usrfrm_cb), num.num_usr);
rx_sym_rbs_eq = zeros(num.num_subc_usr, num.num_ofdmsym_usr, sum(rm.num_usrfrm_cb), num.num_usr);
error_var = zeros(num.num_qamsym_usr, sum(rm.num_usrfrm_cb), num.num_usr);
if ~test_option.awgn
    ch_real_eff_tf = zeros(num.num_subc_usr*num.num_ofdmsym_usr, num.num_subc_usr*num.num_ofdmsym_usr, sum(rm.num_usrfrm_cb), num.num_usr);
    ch_est_rbs = zeros(num.num_subc_usr, num.num_ofdmsym_usr, sum(rm.num_usrfrm_cb), num.num_usr);
end

for idx_usr = 1:num.num_usr
    for idx_usrfrm = 1:sum(rm.num_usrfrm_cb)
        
        % reshape
        rx_ofdmsym_cp = reshape(squeeze(rx_ofdmsym_serial(:, idx_usrfrm, idx_usr)), num.num_fft+num.num_cp, []);
        
        % remove cp
        rx_ofdmsym = rx_ofdmsym_cp(num.num_cp+1:end,:);
        
        % ofdm demodulate
        rx_sym_nfft = (1/sqrt(num.num_fft))*fft(rx_ofdmsym, [], 1);
        
        % demap symbols in bandwidth from fft range
        rx_sym_nfft_shift = fftshift(rx_sym_nfft, 1);
        rx_sym_bw = rx_sym_nfft_shift((num.num_fft/2)-(num.num_subc_bw/2)+1:(num.num_fft/2)+(num.num_subc_bw/2), :);
        
        % demap user resource block from whole bandwidth (temporary)
        usr_id = num.list_usr(idx_usr);                                         % user id
        idx_rb_usr = (ceil(usr_id/num.max_usr_slot)-1)*num.num_rb_usr+1;        % user rb index
        idx_slot_usr = mod(usr_id-1, num.max_usr_slot)*num.num_slot_usr+1;      % user slot index
        list_subc_usr = (idx_rb_usr-1)*num.num_subc_rb+1:(idx_rb_usr-1)*num.num_subc_rb+num.num_subc_usr;                       % user subcarrier list
        list_ofdmsym_usr = (idx_slot_usr-1)*num.num_ofdmsym_slot+1:(idx_slot_usr-1)*num.num_ofdmsym_slot+num.num_ofdmsym_usr;   % user ofdm symbol list
        rx_sym_rbs_usr = rx_sym_bw(list_subc_usr, list_ofdmsym_usr);
        rx_sym_rbs(:, :, idx_usrfrm, idx_usr) = rx_sym_rbs_usr;
        
        if test_option.awgn
            rx_sym_rbs_usr_eq = rx_sym_rbs_usr;
            rx_sym_rbs_eq(:, :, idx_usrfrm, idx_usr) = rx_sym_rbs_usr_eq;
        else
            % demap user real channel (for real channel estimation or channel estimation error check)
            if test_option.ch_mse || strcmp(chest_option, 'real')
                [ch_real_rbs_usr_tf, ch_real_eff_usr_tf, ~] = ...
                    demap_real_ch(squeeze(ch_real_mat_tf(:, :, :, idx_usrfrm, idx_usr)), list_subc_usr, list_ofdmsym_usr, chest_option, cheq_option, test_option);
                
%                 % temp
%                 ch_real_rbs_usr_tf = ones(size(ch_real_rbs_usr_tf));
%                 ch_real_eff_usr_tf = eye(size(ch_real_eff_usr_tf, 1));
                
            else
                ch_real_rbs_usr_tf = [];
                ch_real_eff_usr_tf = [];
            end
            
            % estimate channel (for uplink: channel estimation and equalization after demapping user physical resource block)
            if strcmp(chest_option, 'real')
                ch_est_rbs_usr = ch_real_rbs_usr_tf;
            else
                ch_est_rbs_usr = ofdm_ch_est_r1(squeeze(tx_sym_rbs(:, :, idx_usrfrm, idx_usr)), rx_sym_rbs_usr, squeeze(tx_ofdmsym_faded(:, idx_usrfrm, idx_usr)), list_subc_usr, list_ofdmsym_usr, num, chest_option);
            end
            
%             % recalculate noise variance (considering interference)
%             if strcmp(chest_option, 'real')
% %                 noise_var_intf = noise_var+mean(abs(squeeze(tx_sym_rbs(:, :, idx_usrfrm, idx_usr))).^2, 'all')*mean(abs(sum(ch_real_eff_usr_tf-diag(diag(ch_real_eff_usr_tf)), 2)).^2, 'all');
%                 noise_var_intf = mean(abs(squeeze(tx_sym_rbs(:, :, idx_usrfrm, idx_usr))).^2, 'all')*mean(abs(sum(ch_real_eff_usr_tf-diag(diag(ch_real_eff_usr_tf)), 2)).^2, 'all');
%             else
%                 noise_var_intf = noise_var;
%             end
            
            % equalize channel in tf-domain
%             [rx_sym_rbs_usr_eq, noise_var_adj] = ofdm_ch_eq_r1(rx_sym_rbs_usr, ch_est_rbs_usr, ch_real_eff_usr_tf, num, noise_var_intf, chest_option, cheq_option, test_option);
            [rx_sym_rbs_usr_eq, noise_var_adj] = ofdm_ch_eq_r1(rx_sym_rbs_usr, ch_est_rbs_usr, ch_real_eff_usr_tf, num, noise_var, chest_option, cheq_option, test_option);
            rx_sym_rbs_eq(:, :, idx_usrfrm, idx_usr) = rx_sym_rbs_usr_eq;
            
            % test: to check and calculate channel estimation error
            if test_option.ch_mse
                ch_real_eff_tf(:, :, idx_usrfrm, idx_usr) = ch_real_eff_usr_tf;
                ch_est_rbs(:, :, idx_usrfrm, idx_usr) = ch_est_rbs_usr;
            end
            
%             assignin('base', 'tx_sym_rbs_usr', squeeze(tx_sym_rbs(:, :, idx_usrfrm, idx_usr)))
%             assignin('base', 'rx_sym_rbs_usr', rx_sym_rbs_usr)
%             assignin('base', 'ch_est_rbs_usr', ch_est_rbs_usr)
%             assignin('base', 'noise_var', noise_var)
%             assignin('base', 'noise_var_intf', noise_var_intf)
%             assignin('base', 'ch_real_eff_usr_tf', ch_real_eff_usr_tf)
            
%             assignin('base', 'rx_sym_rbs_usr', rx_sym_rbs_usr)
%             assignin('base', 'ch_est_rbs_usr', ch_est_rbs_usr)
%             assignin('base', 'tx_sym_rbs_usr', tx_sym_rbs_usr)
%             pause
            
        end
        
        % demap data and pilot qam symbols
        [rx_sym_data_usrfrm, ~] = ofdm_sym_demap(rx_sym_rbs_usr_eq, num);
        
        % buffer qam symbols
        rx_sym(:, idx_usrfrm, idx_usr) = rx_sym_data_usrfrm;
        
        % calculate error variance
        if test_option.awgn
            noise_var_usrfrm = noise_var;
        else
            noise_var_usrfrm = ofdm_sym_demap(noise_var_adj, num);
        end
%         error_var_usr = noise_var./(abs(ch_est_rbs_usr).^2);
%         [error_var_usrfrm, ~] = ofdm_sym_demap(error_var_usr, num);
%         error_var(:, idx_usrfrm, idx_usr) = error_var_usrfrm;
        error_var(:, idx_usrfrm, idx_usr) = noise_var_usrfrm;
        
    end
    
%     % compensate channel estimation error variance
%     % qam_mse = mean(abs(tx_sym(:)-rx_sym(:)).^2);
%     % error_var = var(sym_pilot_error(:));    % cannot be used (pilot syms have less error than data syms)
%     % error_var = var(tx_sym_serial-rx_sym_serial);
%     if strcmp(chest_option, 'tf_ltedown') || strcmp(chest_option, 'tf_nr')
%         error_var_usr = noise_var + 0.12;
%     elseif strcmp(chest_option, 'tf_lteup')
%         error_var_usr = noise_var + 0.27;
%     elseif strcmp(chest_option, 'real')
%         error_var_usr = noise_var + 0.08;
%     else
%         error_var_usr = noise_var;
%     end
%     
%     % test: to check channel estimation error
%     if test_option.ch_mse
%         error_var(1, idx_usr) = error_var_usr;
%     end
    
%     % generate ofdm/otfs received bits per user
%     [rx_bit_usr, rx_sym_serial_usr, rx_crc_error_usr] = gen_rx_bit(rx_sym(:, :, idx_usr), rx_crc, turbo_dec, error_var(:, :, idx_usr), sim, cc, rm);
    
    % simulate per-cb rx process
    rx_bit_crc_removed_cb = zeros(cc.K-cc.L, cc.C);
    rx_crc_error_cb = zeros(1, cc.C);
    for idx_cb = 1:cc.C
        
        % serialize qam symbols
        rx_sym_cb = rx_sym(:, sum(rm.num_usrfrm_cb(1:idx_cb-1))+1:sum(rm.num_usrfrm_cb(1:idx_cb)), idx_usr);
        rx_sym_serial_cb = reshape(rx_sym_cb, [], 1);
        rx_sym_serial(sum(rm.num_usrfrm_cb(1:idx_cb-1))*num.num_qamsym_usr+1:sum(rm.num_usrfrm_cb(1:idx_cb))*num.num_qamsym_usr, idx_usr) = rx_sym_serial_cb;
        
        % serialize noise variance
        error_var_cb = reshape(error_var(:, sum(rm.num_usrfrm_cb(1:idx_cb-1))+1:sum(rm.num_usrfrm_cb(1:idx_cb)), idx_usr), [], 1);
%         error_var_cb = abs(rx_sym_serial_cb-tx_sym_serial).^2;
        
        % demap the qam symbols
        % noise_var_tmp = var(tx_sym_serial-rx_sym_serial);
%         rx_bit_demod_cb = (-1) * qamdemod(rx_sym_serial_cb, 2^rm.Qm, 'UnitAveragePower', true, ...
%             'OutputType', 'approxllr', 'NoiseVariance', error_var_cb);
        
        rx_bit_demod_cb = zeros(4, length(rx_sym_serial_cb));
        for i = 1:length(rx_sym_serial_cb)
            rx_bit_demod_cb(:, i) = (-1) * qamdemod(rx_sym_serial_cb(i), 2^rm.Qm, 'UnitAveragePower', true, ...
                'OutputType', 'approxllr', 'NoiseVariance', error_var_cb(i));
        end
        rx_bit_demod_cb = rx_bit_demod_cb(:);
        
        % rate match
        % output fixed (rm.Kw)
        rx_bit_ratematch_cb = rx_ratematch_r2(rx_bit_demod_cb, rm);
        
        % turbo decode data
        rx_bit_dec_cb = turbo_dec(rx_bit_ratematch_cb);
        
        % detect crc
        [rx_bit_crc_removed_cb(:, idx_cb), rx_crc_error_cb(1, idx_cb)] = rx_crc(rx_bit_dec_cb);
    end
    
    % buffer per code block
    rx_bit_buff_usr = reshape(rx_bit_crc_removed_cb, [], 1);
    
    % remove padded bits
    rx_bit(:, idx_usr) = rx_bit_buff_usr(1 : sim.len_tb_bit);
    
    % rearrange crc
    rx_crc_error(:, idx_usr) = sum(double(rx_crc_error_cb)) > 0;
end

% calculate packet error
pkt_error = rx_crc_error | (symerr(tx_bit, rx_bit, 'column-wise') > 0);

% calculate papr
if test_option.papr
    tx_papr = mean(max(abs(tx_ofdmsym_serial).^2, [], 1)./mean(abs(tx_ofdmsym_serial).^2, 1));
else
    tx_papr = [];
end

% test: calculate channel estimation error
if ~test_option.awgn && test_option.ch_mse
    if strcmp(chest_option, 'real') && test_option.fulltap_eq
        ch_mse = zeros(1, num.num_usr);
    else
        % generate effective tf channel (per time-frequency grid)
        ch_est_eff = zeros(num.num_subc_usr*num.num_ofdmsym_usr, num.num_subc_usr*num.num_ofdmsym_usr, sum(rm.num_usrfrm_cb), num.num_usr);
        for idx_usr = 1:num.num_usr
            for idx_usrfrm = 1:sum(rm.num_usrfrm_cb)
                ch_est_eff(:, :, idx_usrfrm, idx_usr) = diag(reshape(ch_est_rbs(:, :, idx_usrfrm, idx_usr), [], 1));
            end
        end
        
        % calculate channel estimation error
        ch_mse = mean(abs(ch_real_eff_tf-ch_est_eff).^2, [1 2 3])*num.num_subc_usr*num.num_ofdmsym_usr;
        
%         fprintf('ch mse:')
%         fprintf(' %10.6f', ch_mse)
%         fprintf('\n')
%         assignin('base', 'ch_real_eff_tf', ch_real_eff_tf)
%         assignin('base', 'ch_est_eff', ch_est_eff)
    end
    for idx_usr=1:num.num_usr
        fprintf('user %d:    tf channel est rmse: %6.3f\n', idx_usr, sqrt(ch_mse(1, idx_usr)))
%         for idx_usrfrm = 1:rm.num_usrfrm
%             test_ch_real_onetap_tf = reshape(diag(ch_real_eff_tf(:, :, idx_usrfrm, idx_usr)), num.num_subc_usr, num.num_ofdmsym_usr);
%             if strcmp(chest_option, 'real') && test_option.fulltap_eq
%                 test_ch_est_onetap_tf = test_ch_real_onetap_tf;
%             else
%                 test_ch_est_onetap_tf = ch_est_rbs(:, :, idx_usrfrm, idx_usr);
%             end
%             test_tx_sym = tx_sym_rbs(:, :, idx_usrfrm, idx_usr);
%             test_rx_sym = rx_sym_rbs(:, :, idx_usrfrm, idx_usr);
%             test_rx_sym_eq = rx_sym_rbs_eq(:, :, idx_usrfrm, idx_usr);
%             
%             figure
%             subplot(2, 2, 1), mesh(1:size(test_ch_real_onetap_tf, 2), 1:size(test_ch_real_onetap_tf, 1), real(test_ch_real_onetap_tf))
%             subplot(2, 2, 3), mesh(1:size(test_ch_real_onetap_tf, 2), 1:size(test_ch_real_onetap_tf, 1), imag(test_ch_real_onetap_tf))
%             subplot(2, 2, 2), mesh(1:size(test_ch_est_onetap_tf, 2), 1:size(test_ch_est_onetap_tf, 1), real(test_ch_est_onetap_tf))
%             subplot(2, 2, 4), mesh(1:size(test_ch_est_onetap_tf, 2), 1:size(test_ch_est_onetap_tf, 1), imag(test_ch_est_onetap_tf))
%             
%             figure
%             subplot(2, 3, 1), mesh(1:size(test_tx_sym, 2), 1:size(test_tx_sym, 1), real(test_tx_sym))
%             subplot(2, 3, 4), mesh(1:size(test_tx_sym, 2), 1:size(test_tx_sym, 1), imag(test_tx_sym))
%             subplot(2, 3, 2), mesh(1:size(test_rx_sym, 2), 1:size(test_rx_sym, 1), real(test_rx_sym))
%             subplot(2, 3, 5), mesh(1:size(test_rx_sym, 2), 1:size(test_rx_sym, 1), imag(test_rx_sym))
%             subplot(2, 3, 3), mesh(1:size(test_rx_sym_eq, 2), 1:size(test_rx_sym_eq, 1), real(test_rx_sym_eq))
%             subplot(2, 3, 6), mesh(1:size(test_rx_sym_eq, 2), 1:size(test_rx_sym_eq, 1), imag(test_rx_sym_eq))
%             pause
%         end
    end
else
    ch_mse = [];
end

% test channel estimation error variance
if test_option.sym_err_var
    sym_err_var = var(tx_sym_serial-rx_sym_serial, 0, 1);
    
    for idx_usr = 1:num.num_usr
        fprintf('user %d:    noise_var:%6.3f    sym_err_var:%6.3f    calc_err_var:%6.3f\n', idx_usr, noise_var, sym_err_var(1, idx_usr), error_var(1, idx_usr))
        
%         figure
%         subplot(2, 1, 1), plot(real(tx_sym_serial(:, idx_usr)), '-b.'), hold on, plot(real(rx_sym_serial(:, idx_usr)), '-r'), hold off, grid minor
%         subplot(2, 1, 2), plot(imag(tx_sym_serial(:, idx_usr)), '-b.'), hold on, plot(imag(rx_sym_serial(:, idx_usr)), '-r'), hold off, grid minor
%         pause
    end
else
    sym_err_var = [];
end

% dump
assignin('base', 'tx_bit', tx_bit);
assignin('base', 'tx_sym', tx_sym);
assignin('base', 'tx_sym_rbs', tx_sym_rbs);
assignin('base', 'tx_sym_bw', tx_sym_bw);
assignin('base', 'tx_sym_nfft_shift', tx_sym_nfft_shift);
assignin('base', 'tx_sym_nfft', tx_sym_nfft);
assignin('base', 'tx_ofdmsym', tx_ofdmsym);
assignin('base', 'tx_ofdmsym_cp', tx_ofdmsym_cp);
assignin('base', 'tx_ofdmsym_serial', tx_ofdmsym_serial);
assignin('base', 'tx_ofdmsym_faded', tx_ofdmsym_faded);

assignin('base', 'rx_ofdmsym_serial', rx_ofdmsym_serial);
assignin('base', 'rx_ofdmsym_cp', rx_ofdmsym_cp);
assignin('base', 'rx_ofdmsym', rx_ofdmsym);
assignin('base', 'rx_sym_nfft', rx_sym_nfft);
assignin('base', 'rx_sym_nfft_shift', rx_sym_nfft_shift);
assignin('base', 'rx_sym_bw', rx_sym_bw);
assignin('base', 'rx_sym_rbs', rx_sym_rbs);
assignin('base', 'rx_sym_rbs_eq', rx_sym_rbs_eq);
assignin('base', 'rx_sym', rx_sym);
assignin('base', 'rx_bit', rx_bit);

assignin('base', 'fading_ch', fading_ch)
assignin('base', 'ch_real_rbs_usr_tf', ch_real_rbs_usr_tf);
assignin('base', 'ch_real_eff_usr_tf', ch_real_eff_usr_tf);
assignin('base', 'ch_real_eff_tf', ch_real_eff_tf);
assignin('base', 'ch_est_rbs_usr', ch_est_rbs_usr);
assignin('base', 'ch_est_rbs', ch_est_rbs);
assignin('base', 'ch_path_gain_usr', ch_path_gain_usr);
assignin('base', 'snr_db_adj', snr_db_adj);
assignin('base', 'error_var', error_var)
assignin('base', 'noise_var', noise_var)
% assignin('base', 'noise_var_intf', noise_var_intf)
assignin('base', 'noise_var_adj', noise_var_adj)
assignin('base', 'noise_var_usrfrm', noise_var_usrfrm)

end
