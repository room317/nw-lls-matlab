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

function pkt_error = otfs_single_run_r3(sim, cc, rm, num, snr_db, ch, chest_option, cheq_option)

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
        'NormalizePathGains', true);
%         'DirectPathInitialPhase', 0.5,...
%         'PathGainsOutputPort', true);
else
    fading_ch = comm.RayleighChannel(...
        'SampleRate', num.sample_rate, ...
        'PathDelays', ch.path_delays, ...
        'AveragePathGains', ch.average_path_gains, ...
        'NormalizePathGains', true, ...
        'MaximumDopplerShift', ch.maximum_doppler_shift, ...
        'DopplerSpectrum', ch.doppler_spectrum);
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

% generate crc
tx_crc = crc.generator(cc.gCRC24A);
tx_bit_crc = generate(tx_crc, tx_bit_pad);

% turbo encode data
tx_bit_enc = turbo_enc(tx_bit_crc);

% rate match
tx_bit_ratematch = tx_ratematch(tx_bit_enc, rm);

% modulate bit stream
tx_sym = qammod(tx_bit_ratematch, 2^rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);

% simulate per subframe
tx_sym_subfrm = reshape(tx_sym, num.num_qamsym_usr, []);
rx_sym_subfrm = zeros(size(tx_sym_subfrm));
for subfrm_idx = 1 : size(tx_sym_subfrm, 2)
    
    if strcmp(chest_option, 'real')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% pilot frame to get real channel %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % generate single tone in dd domain
        tx_sym_dd_ndft = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
        tx_sym_dd_ndft(1, 1) = sqrt(num.num_subc_usr*num.num_ofdmsym_usr);
        
        % 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
        tx_sym_tf_ndft = sqrt(num.num_ofdmsym_usr/num.num_subc_usr)*fft(ifft(tx_sym_dd_ndft, [], 2), [], 1);
        
        % map otfs symbols to resource block (otfs)
        tx_sym_tf_subc = tx_sym_tf_ndft;
        
        % map to the resource block (nfft x symbols)
        tx_sym_tf_nfft = zeros(num.nfft, num.num_ofdmsym_usr);
        tx_sym_tf_nfft((num.nfft/2)-(num.num_subc_bw/2)+1:(num.nfft/2)+(num.num_subc_bw/2), :) = tx_sym_tf_subc;
        tx_sym_tf_nfft_shift = fftshift(tx_sym_tf_nfft, 1);
        
        % ofdm modulate
        tx_ofdm_sym = sqrt(num.nfft) * ifft(tx_sym_tf_nfft_shift, [], 1);
        
        % add cp
        tx_ofdm_sym_cp = tx_ofdm_sym([num.nfft-num.num_cp+1:num.nfft 1:num.nfft], :);
        
        % serialize
        tx_ofdm_sym_serial = tx_ofdm_sym_cp(:);
        
        % pass signal through channel
        random_seed = rng;
        tx_ofdm_sym_faded = fading_ch(tx_ofdm_sym_serial);
        
        % add gaussian noise
        rx_ofdm_sym_serial = tx_ofdm_sym_faded;
        
        % estimate real channel in tf-domain
        ch_tf_ndft_real = otfs_ch_est_tf(tx_ofdm_sym_serial, rx_ofdm_sym_serial, num, chest_option);
        
        if strcmp(cheq_option, 'ddeq')
            % reshape
            rx_ofdm_sym_cp = reshape(rx_ofdm_sym_serial, num.nfft+num.num_cp, num.num_ofdmsym_usr);
            
            % remove cp
            rx_ofdm_sym = rx_ofdm_sym_cp(num.num_cp+1:end,:);
            
            % ofdm demodulate the symbol
            rx_sym_tf_nfft = (1/sqrt(num.nfft)) * fft(rx_ofdm_sym, [], 1);
            
            % demap subcarriers (subcarrier x symbols)
            rx_sym_tf_nfft_shift = fftshift(rx_sym_tf_nfft, 1);
            rx_sym_tf_subc = rx_sym_tf_nfft_shift(num.nfft/2-num.num_subc_bw/2+1:num.nfft/2+num.num_subc_bw/2, :);
            
            % demap symbols from subcarriers
            rx_sym_tf_ndft = rx_sym_tf_subc;
            
            % 2d inverse sfft
            rx_sym_dd_ndft_subfrm = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(rx_sym_tf_ndft, [], 1), [], 2);
            
            % estimate real channel in dd-domain
            ch_dd_ndft_real = rx_sym_dd_ndft_subfrm;
        else
            ch_dd_ndft_real = [];
        end
    else
        ch_tf_ndft_real = [];
        ch_dd_ndft_real = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%
    %%% normal frame %%%
    %%%%%%%%%%%%%%%%%%%%
    
    % map data and pilot symbols
    if strcmp(chest_option, 'dd_singletone')
%         [tx_sym_dd_ndft, idx_pilot_sym] = otfs_sym_map(tx_sym_subfrm(:, subfrm_idx), num);
        tx_sym_dd_ndft = otfs_sym_map_r1(tx_sym_subfrm(:, subfrm_idx), num, 2);
    else
        % generate pilot symbols
        tx_sym_dd_ndft_pilot = randi([0 1], num.num_subc_usr, num.num_ofdmsym_pilot_usr)*2-1;
        
        % reshape data symbols
        tx_sym_dd_ndft_data = reshape(tx_sym_subfrm(:, subfrm_idx), num.num_subc_usr, num.num_ofdmsym_data_usr);
    
        % map data and pilot symbols
        tx_sym_dd_ndft = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
        tx_sym_dd_ndft(:, num.idx_ofdmsym_data_usr) = tx_sym_dd_ndft_data;
        tx_sym_dd_ndft(:, num.idx_ofdmsym_pilot_usr) = tx_sym_dd_ndft_pilot;
    end
    
    % 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
    tx_sym_tf_ndft = sqrt(num.num_ofdmsym_usr/num.num_subc_usr)*fft(ifft(tx_sym_dd_ndft, [], 2), [], 1);
    
    % map otfs symbols to resource block (otfs)
    tx_sym_tf_subc = tx_sym_tf_ndft;
    
    % map to the resource block (nfft x symbols)
    tx_sym_tf_nfft = zeros(num.nfft, num.num_ofdmsym_usr);
    tx_sym_tf_nfft((num.nfft/2)-(num.num_subc_bw/2)+1:(num.nfft/2)+(num.num_subc_bw/2), :) = tx_sym_tf_subc;
    tx_sym_tf_nfft_shift = fftshift(tx_sym_tf_nfft, 1);
    
    % ofdm modulate
    tx_ofdm_sym = sqrt(num.nfft) * ifft(tx_sym_tf_nfft_shift, [], 1);
    
    % add cp
    tx_ofdm_sym_cp = tx_ofdm_sym([num.nfft-num.num_cp+1:num.nfft 1:num.nfft], :);
    
    % serialize
    tx_ofdm_sym_serial = tx_ofdm_sym_cp(:);
    
    % pass signal through channel
    if strcmp(chest_option, 'real')
        rng(random_seed)
    end
    tx_ofdm_sym_faded = fading_ch(tx_ofdm_sym_serial);
    
    % add gaussian noise
    rx_ofdm_sym_serial = awgn(tx_ofdm_sym_faded, snr_db, 'measured');
    
    % reshape
    rx_ofdm_sym_cp = reshape(rx_ofdm_sym_serial, num.nfft+num.num_cp, num.num_ofdmsym_usr);
    
    % remove cp
    rx_ofdm_sym = rx_ofdm_sym_cp(num.num_cp+1:end,:);
    
    % ofdm demodulate the symbol
    rx_sym_tf_nfft = (1/sqrt(num.nfft)) * fft(rx_ofdm_sym, [], 1);
    
    % demap subcarriers (subcarrier x symbols)
    rx_sym_tf_nfft_shift = fftshift(rx_sym_tf_nfft, 1);
    rx_sym_tf_subc = rx_sym_tf_nfft_shift(num.nfft/2-num.num_subc_bw/2+1:num.nfft/2+num.num_subc_bw/2, :);
    
    % demap symbols from subcarriers
    rx_sym_tf_ndft = rx_sym_tf_subc;
    
    % estimate and equalize channel
    if strcmp(cheq_option, 'tfeq_zf') || strcmp(cheq_option, 'tfeq_mmse')
        % estimate channel
        if strcmp(chest_option, 'real')
            % get real channel
            ch_est_tf_ndft = ch_tf_ndft_real;
        elseif  strcmp(chest_option, 'perfect')
            % estimate channel in tf-domain
            ch_est_tf_ndft = otfs_ch_est_tf(tx_ofdm_sym_serial, tx_ofdm_sym_faded, num, chest_option);
        elseif  strcmp(chest_option, 'tf_lteup')
            % estimate channel in tf-domain
            ch_est_tf_ndft = otfs_ch_est_tf(tx_ofdm_sym_serial, rx_ofdm_sym_serial, num, chest_option);
        else
            % 2d inverse sfft
            rx_sym_dd_ndft = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(rx_sym_tf_ndft, [], 1), [], 2);
            
            % estimate channel in dd-domain
            ch_est_dd_ndft = otfs_ch_est_dd(rx_sym_dd_ndft, num);
            
            % transform dd channel to tf channel
            ch_est_tf_ndft = sqrt(num.num_ofdmsym_usr/num.num_subc_usr)*fft(ifft(ch_est_dd_ndft, [], 2), [], 1);
        end
        
        % equalize channel in tf-domain
        rx_sym_tf_ndft_eq = otfs_ch_eq_tf(rx_sym_tf_ndft, ch_est_tf_ndft, noise_var, cheq_option);
        
        % 2d inverse sfft
        rx_sym_dd_ndft_eq = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(rx_sym_tf_ndft_eq, [], 1), [], 2);
%         rx_sym_dd_ndft_eq = circshift(rx_sym_dd_ndft_eq, num.num_delay_guard);
        
%         figure, mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, abs(tx_sym_dd_ndft))
%         figure, mesh(1:num.num_ofdmsym_usr, 1:num.nfft+num.num_cp, abs(rx_ofdm_sym_cp))
%         figure, mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, abs(rx_sym_dd_ndft))
%         figure, mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, abs(rx_sym_dd_ndft_eq))
%         figure, plot(1:num.num_subc_usr, abs(tx_sym_dd_ndft))
%         figure, plot(1:num.nfft+num.num_cp, abs(tx_ofdm_sym_cp))
%         figure, plot(1:num.nfft+num.num_cp, abs(rx_ofdm_sym_cp))
%         figure, plot(1:num.num_subc_usr, abs(rx_sym_dd_ndft))
%         figure, plot(1:num.num_subc_usr, abs(rx_sym_dd_ndft_eq))
%         figure, plot(1:num.num_subc_usr, abs(circshift(rx_sym_dd_ndft_eq, num.num_delay_guard)))
%         pause
    else
        % 2d inverse sfft
        rx_sym_dd_ndft = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(rx_sym_tf_ndft, [], 1), [], 2);
        
        % estimate channel
        if strcmp(chest_option, 'real')
            % get real channel
            ch_est_dd_ndft = ch_dd_ndft_real;
        elseif  strcmp(chest_option, 'dd_singletone')
            % estimate channel in dd-domain
            ch_est_dd_ndft = otfs_ch_est_dd(rx_sym_dd_ndft, num);
        else
            ch_est_dd_ndft = [];
        end
        
        % equalize channel in dd-domain
        rx_sym_dd_ndft_eq = otfs_ch_eq_dd(rx_sym_dd_ndft, ch_est_dd_ndft, num);
    end
    
    % demap data qam symbols
    if strcmp(chest_option, 'dd_singletone')
        [rx_sym_dd_ndft_data, ~] = otfs_sym_demap_r1(rx_sym_dd_ndft_eq, num, 2);
    else
        rx_sym_dd_ndft_data = rx_sym_dd_ndft_eq(:, num.idx_ofdmsym_data_usr);
    end
    
    % buffer qam symbols
    rx_sym_subfrm(:, subfrm_idx) = rx_sym_dd_ndft_data(:);
end

% serialize qam symbols
rx_sym = rx_sym_subfrm(:);

% demap the qam symbols
noise_var_tmp = var(tx_sym(:)-rx_sym(:));
rx_bit_demod = (-1) * qamdemod(rx_sym(:), 2^rm.Qm, 'UnitAveragePower', true, ...
    'OutputType', 'approxllr', 'NoiseVariance', noise_var_tmp);

% rate match
rx_bit_ratematch = rx_ratematch(rx_bit_demod, rm);

% turbo decode data
rx_bit_dec = turbo_dec(rx_bit_ratematch);

% detect crc
rx_crc = crc.detector(cc.gCRC24A);
[rx_bit_pad, rx_crc_error] = detect(rx_crc, rx_bit_dec);

% remove padded bits
rx_bit = rx_bit_pad(1 : sim.len_tb_bit);

if rx_crc_error
    pkt_error = true;
else
    pkt_error = symerr(tx_bit, rx_bit) > 0;
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

% dump
assignin('base', 'tx_sym_dd_ndft', tx_sym_dd_ndft);
assignin('base', 'tx_sym_tf_ndft', tx_sym_tf_ndft);
assignin('base', 'tx_sym_tf_subc', tx_sym_tf_subc);
assignin('base', 'tx_sym_tf_nfft', tx_sym_tf_nfft);
assignin('base', 'tx_sym_tf_nfft_shift', tx_sym_tf_nfft_shift);
assignin('base', 'tx_ofdm_sym', tx_ofdm_sym);
assignin('base', 'tx_ofdm_sym_cp', tx_ofdm_sym_cp);
assignin('base', 'tx_ofdm_sym_serial', tx_ofdm_sym_serial);
assignin('base', 'tx_ofdm_sym_faded', tx_ofdm_sym_faded);
assignin('base', 'rx_ofdm_sym_serial', rx_ofdm_sym_serial);
assignin('base', 'rx_ofdm_sym_cp', rx_ofdm_sym_cp);
assignin('base', 'rx_ofdm_sym', rx_ofdm_sym);
assignin('base', 'rx_sym_tf_nfft', rx_sym_tf_nfft);
assignin('base', 'rx_sym_tf_nfft_shift', rx_sym_tf_nfft_shift);
assignin('base', 'rx_sym_tf_subc', rx_sym_tf_subc);
assignin('base', 'rx_sym_tf_ndft', rx_sym_tf_ndft);
assignin('base', 'rx_sym_dd_ndft_eq', rx_sym_dd_ndft_eq);

end
