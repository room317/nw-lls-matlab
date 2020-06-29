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
% - subframe buffer fixed: 2019.12.10
% - eq_dd(delay-doppler channel estimation) added: 2020.01.09
% - test_dd_pilot added: 2020.01.09
% - real channel estimation with single tone pilot added: 2020.02.19
% - real channel compensation in tf domain added: 2020.02.19
% - real channel compensation in dd domain added: 2020.02.19
% - variable name changed: 2020.06.09

function pkt_error = otfs_single_run_r4(sim, cc, rm, num, snr_db, ch, chest_option, cheq_option, map_plan)

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
noise_var = (num.num_subc_usr/num.nfft)*(10 ^ ((-0.1)*snr_db));

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
% tx_bit_ratematch = tx_ratematch(tx_bit_enc, rm);
tx_bit_ratematch = tx_ratematch_r2(tx_bit_enc, rm);

% modulate bit stream
tx_sym_serial = qammod(tx_bit_ratematch(:), 2^rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);

% simulate per subframe
tx_sym = reshape(tx_sym_serial, num.num_qamsym_usr, []);
rx_sym = zeros(size(tx_sym));
for idx_subfrm = 1 : size(tx_sym, 2)
    
    % extract subframe data per user
    tx_sym_data_subfrm = tx_sym(:, idx_subfrm);
    
    % map data and pilot symbols
    tx_sym_rbs_dd = otfs_sym_map_r2(tx_sym_data_subfrm, num, map_plan);
    
    % 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
    tx_sym_rbs_tf = sqrt(num.num_ofdmsym_usr/num.num_subc_usr)*fft(ifft(tx_sym_rbs_dd, [], 2), [], 1);
    
    % map user resource blocks to whole bandwidth (temporary)
    tx_sym_bw = zeros(num.num_subc_bw, num.num_ofdmsym_subfrm);
    tx_sym_bw(1:size(tx_sym_rbs_tf, 1), 1:size(tx_sym_rbs_tf, 2)) = tx_sym_rbs_tf;
    
    % map bandwidth to the center of fft range
    tx_sym_nfft_shift = zeros(num.nfft, num.num_ofdmsym_subfrm);
    tx_sym_nfft_shift((num.nfft/2)-(num.num_subc_bw/2)+1:(num.nfft/2)+(num.num_subc_bw/2), :) = tx_sym_bw;
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
    
    % reshape
    rx_ofdmsym_cp = reshape(rx_ofdmsym_serial, num.nfft+num.num_cp, num.num_ofdmsym_subfrm);
    
    % remove cp
    rx_ofdmsym = rx_ofdmsym_cp(num.num_cp+1:end,:);
    
    % ofdm demodulate
    rx_sym_nfft = (1/sqrt(num.nfft)) * fft(rx_ofdmsym, [], 1);
    
    % demap symbols in bandwidth from fft range
    rx_sym_nfft_shift = fftshift(rx_sym_nfft, 1);
    rx_sym_bw = rx_sym_nfft_shift((num.nfft/2)-(num.num_subc_bw/2)+1:(num.nfft/2)+(num.num_subc_bw/2), :);
    
    % demap user resource block from whole bandwidth
    rx_sym_rbs_tf = rx_sym_bw(1:num.num_delay_usr, 1:num.num_doppler_usr);
    
    % estimate and equalize channel
    if strcmp(chest_option, 'dd_tone') && strcmp(cheq_option, 'ddeq')
        
        % 2d inverse sfft both for channel estimation and demodulation
        rx_sym_rbs_dd = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(rx_sym_rbs_tf, [], 1), [], 2);
        
        % estimate channel
        ch_est_rbs_dd = otfs_ch_est_dd_r1(rx_sym_rbs_dd, num, map_plan);
        
        % equalize channel
        rx_sym_rbs_dd_eq = otfs_ch_eq_dd_r2(rx_sym_rbs_dd, ch_est_rbs_dd, num);
        
    elseif strcmp(chest_option, 'dd_tone') && ~strcmp(cheq_option, 'ddeq')
        
        % 2d inverse sfft for channel estimation only
        rx_sym_rbs_dd = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(rx_sym_rbs_tf, [], 1), [], 2);
        
        % estimate channel
        ch_est_rbs_dd = otfs_ch_est_dd_r1(rx_sym_rbs_dd, num, map_plan);
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
    
    % demap data qam symbols
    [rx_sym_data_subfrm, ~] = otfs_sym_demap_r2(rx_sym_rbs_dd_eq, num, map_plan);
    
    % buffer qam symbols
    rx_sym(:, idx_subfrm) = rx_sym_data_subfrm(:);
end

% serialize qam symbols
rx_sym_serial = rx_sym(:);

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
% assignin('base', 'tx_sym_dd_ndft', tx_sym_rbs_dd);
% assignin('base', 'tx_sym_tf_ndft', tx_sym_rbs_tf);
% assignin('base', 'tx_sym_tf_subc', tx_sym_bw);
% assignin('base', 'tx_sym_tf_nfft', tx_sym_nfft_shift);
% assignin('base', 'tx_sym_tf_nfft_shift', tx_sym_nfft);
% assignin('base', 'tx_ofdm_sym', tx_ofdmsym);
% assignin('base', 'tx_ofdm_sym_cp', tx_ofdmsym_cp);
% assignin('base', 'tx_ofdm_sym_serial', tx_ofdmsym_serial);
% assignin('base', 'tx_ofdm_sym_faded', tx_ofdmsym_faded);
% assignin('base', 'rx_ofdm_sym_serial', rx_ofdmsym_serial);
% assignin('base', 'rx_ofdm_sym_cp', rx_ofdmsym_cp);
% assignin('base', 'rx_ofdm_sym', rx_ofdmsym);
% assignin('base', 'rx_sym_tf_nfft', rx_sym_nfft);
% assignin('base', 'rx_sym_tf_nfft_shift', rx_sym_nfft_shift);
% assignin('base', 'rx_sym_tf_subc', rx_sym_bw);
% assignin('base', 'rx_sym_tf_ndft', rx_sym_rbs);
% assignin('base', 'rx_sym_dd_ndft_eq', rx_sym_dd_ndft_eq);

end
