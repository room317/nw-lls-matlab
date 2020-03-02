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

function pkt_error = otfs_single_run_r1(sim, cc, rm, num, snr_db, ch, eq_dd, test_dd_pilot, test_real_ch)

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
noise_var = 10 ^ ((-0.1)*snr_db);

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
tx_sym_subfrm = reshape(tx_sym, num.len_rb_sym_user, []);
rx_sym_dd_ndft = zeros(size(tx_sym_subfrm));
for subfrm_idx = 1 : size(tx_sym_subfrm, 2)
    
    if test_real_ch
        %%%%%%%%%%%%%%%%%%%%
        %%%%%%%% pilot frame
        %%%%%%%%%%%%%%%%%%%%
        
        % generate single tone in dd domain
        tx_sym_dd_ndft = zeros(num.ndft, num.num_ofdmsym_per_subframe);
        tx_sym_dd_ndft((num.ndft/2)+1, (num.num_ofdmsym_per_subframe/2)+1) = sqrt(num.ndft*num.num_ofdmsym_per_subframe);
        
        % 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
        tx_sym_tf_ndft = sqrt(num.num_ofdmsym_per_subframe/num.ndft)*fft(ifft(tx_sym_dd_ndft, [], 2), [], 1);
        
        % map otfs symbols to resource block (otfs)
        tx_sym_tf_subc = tx_sym_tf_ndft;
        
        % map to the resource block (nfft x symbols)
        tx_sym_tf_nfft = zeros(num.nfft, num.num_ofdmsym_per_subframe);
        tx_sym_tf_nfft((num.nfft/2)-(num.num_subcarrier/2)+1:(num.nfft/2)+(num.num_subcarrier/2), :) = tx_sym_tf_subc;
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
        
        % reshape
        rx_ofdm_sym_cp = reshape(rx_ofdm_sym_serial, num.nfft+num.num_cp, num.num_ofdmsym_per_subframe);
        
        % remove cp
        rx_ofdm_sym = rx_ofdm_sym_cp(num.num_cp+1:end,:);
%         assignin('base', 'rx_ofdm_sym', rx_ofdm_sym);
%         figure, mesh(1:num.num_ofdmsym_per_subframe, 1:num.nfft, abs(rx_ofdm_sym))
%         pause
        
        % ofdm demodulate the symbol
        rx_sym_tf_nfft = (1/sqrt(num.nfft)) * fft(rx_ofdm_sym, [], 1);
        
        % demap subcarriers (subcarrier x symbols)
        rx_sym_tf_nfft_shift = fftshift(rx_sym_tf_nfft, 1);
        rx_sym_tf_subc = rx_sym_tf_nfft_shift(num.nfft/2-num.num_subcarrier/2+1:num.nfft/2+num.num_subcarrier/2, :);
        
        % demap symbols from subcarriers
        rx_sym_tf_ndft = rx_sym_tf_subc;
        
        % estimate real channel
        [~, ch_tf_real] = otfs_ch_comp_tf_mmse_perfect(rx_sym_tf_ndft, tx_ofdm_sym_serial, tx_ofdm_sym_faded, false, [], num, noise_var);
        
        % 2d inverse sfft
        ch_dd_real_shift = sqrt(num.ndft/num.num_ofdmsym_per_subframe)*fft(ifft(ch_tf_real, [], 1), [], 2);
        ch_dd_real = fftshift(fftshift(ch_dd_real_shift, 1), 2);
    else
        ch_tf_real = [];
        ch_dd_real = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% normal frame
    %%%%%%%%%%%%%%%%%%%%%
    
    % map data and pilot symbols
    if eq_dd
        [tx_sym_dd_ndft, idx_pilot_sym] = otfs_sym_map(tx_sym_subfrm(:, subfrm_idx), num, test_dd_pilot);
    else
        % generate pilot symbols
        tx_sym_dd_ndft_pilot = randi([0 1], num.ndft, num.num_pilot_ofdmsym_per_subframe)*2-1;
        
        % reshape data symbols
        tx_sym_dd_ndft_data = reshape(tx_sym_subfrm(:, subfrm_idx), num.ndft, num.num_data_ofdmsym_per_subframe);
    
        % map data and pilot symbols
        tx_sym_dd_ndft = zeros(num.ndft, num.num_ofdmsym_per_subframe);
        tx_sym_dd_ndft(:, num.idx_data_ofdmsym) = tx_sym_dd_ndft_data;
        tx_sym_dd_ndft(:, num.idx_pilot_ofdmsym) = tx_sym_dd_ndft_pilot;
    end
    
    % 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
    tx_sym_tf_ndft = sqrt(num.num_ofdmsym_per_subframe/num.ndft)*fft(ifft(tx_sym_dd_ndft, [], 2), [], 1);
    
    % map otfs symbols to resource block (otfs)
    tx_sym_tf_subc = tx_sym_tf_ndft;
    
    % map to the resource block (nfft x symbols)
    tx_sym_tf_nfft = zeros(num.nfft, num.num_ofdmsym_per_subframe);
    tx_sym_tf_nfft((num.nfft/2)-(num.num_subcarrier/2)+1:(num.nfft/2)+(num.num_subcarrier/2), :) = tx_sym_tf_subc;
    tx_sym_tf_nfft_shift = fftshift(tx_sym_tf_nfft, 1);
    
    % ofdm modulate
    tx_ofdm_sym = sqrt(num.nfft) * ifft(tx_sym_tf_nfft_shift, [], 1);
    
    % add cp
    tx_ofdm_sym_cp = tx_ofdm_sym([num.nfft-num.num_cp+1:num.nfft 1:num.nfft], :);
    
    % serialize
    tx_ofdm_sym_serial = tx_ofdm_sym_cp(:);
    
    % pass signal through channel
    if test_real_ch
        rng(random_seed)
    end
    tx_ofdm_sym_faded = fading_ch(tx_ofdm_sym_serial);
    
    % add gaussian noise
    rx_ofdm_sym_serial = awgn(tx_ofdm_sym_faded, snr_db, 'measured');
%     rx_ofdm_sym_serial = awgn(tx_ofdm_sym_serial, snr_db, 'measured');
%     rx_ofdm_sym_serial = tx_ofdm_sym_faded;
    
    % reshape
    rx_ofdm_sym_cp = reshape(rx_ofdm_sym_serial, num.nfft+num.num_cp, num.num_ofdmsym_per_subframe);
    
    % remove cp
    rx_ofdm_sym = rx_ofdm_sym_cp(num.num_cp+1:end,:);
    
    % ofdm demodulate the symbol
    rx_sym_tf_nfft = (1/sqrt(num.nfft)) * fft(rx_ofdm_sym, [], 1);
    
    % demap subcarriers (subcarrier x symbols)
    rx_sym_tf_nfft_shift = fftshift(rx_sym_tf_nfft, 1);
    rx_sym_tf_subc = rx_sym_tf_nfft_shift(num.nfft/2-num.num_subcarrier/2+1:num.nfft/2+num.num_subcarrier/2, :);
    
    % demap symbols from subcarriers
    rx_sym_tf_ndft = rx_sym_tf_subc;
    
    % estimate and compensate channel
    if eq_dd
        % 2d inverse sfft
        rx_sym_dd_ndft_subfrm = sqrt(num.ndft/num.num_ofdmsym_per_subframe)*fft(ifft(rx_sym_tf_ndft, [], 1), [], 2);
        
        % estimate and compensate channel in dd domain
        [rx_sym_dd_ndft_eq_subfrm, ~, ~] = otfs_ch_comp_dd_mmse(rx_sym_dd_ndft_subfrm, rx_sym_tf_ndft, idx_pilot_sym, num, noise_var, test_dd_pilot, test_real_ch, ch_dd_real);
%         [rx_sym_dd_ndft_eq_subfrm, aaaaaa, ~] = otfs_ch_comp_dd_mmse(rx_sym_dd_ndft_subfrm, rx_sym_tf_ndft, idx_pilot_sym, num, noise_var, test_dd_pilot, false, ch_dd_real);
        
%         % for test (comment this block before simulation!)
%         [~, ch_est_dd_ndft_shift_subfrm, ch_est_tf_ndft_subfrm] = otfs_ch_comp_dd_mmse(rx_sym_dd_ndft_subfrm, rx_sym_tf_ndft, idx_pilot_sym, num, noise_var, test_dd_pilot);
%         [ch_est_mse_tf_subfrm, ch_est_mse_dd_subfrm] = test_otfs_real_ch(tx_ofdm_sym_serial, tx_ofdm_sym_faded, ch_est_dd_ndft_shift_subfrm, ch_est_tf_ndft_subfrm, num);
%         fprintf('channel estimation error in tf: %5.3f    channel estimation error in dd: %5.3f\n', ch_est_mse_tf_subfrm, ch_est_mse_dd_subfrm)
%         pause
        
        % demap data qam symbols
        rx_sym_dd_ndft_data_subfrm = otfs_sym_demap(rx_sym_dd_ndft_eq_subfrm, idx_pilot_sym, num);
    else
        % estimate and compensate channel in tf domain
        [rx_sym_tf_ndft_eq, ~] = otfs_ch_comp_tf_mmse_perfect(rx_sym_tf_ndft, tx_ofdm_sym_serial, tx_ofdm_sym_faded, test_real_ch, ch_tf_real, num, noise_var);
%         [rx_sym_tf_ndft_eq, aaaaaa] = otfs_ch_comp_tf_mmse_perfect(rx_sym_tf_ndft, tx_ofdm_sym_serial, tx_ofdm_sym_faded, false, ch_tf_real, num, noise_var);
%         [rx_sym_tf_ndft_eq, ~] = otfs_ch_comp_tf_mmse(rx_sym_tf_ndft, tx_ofdm_sym_serial, tx_ofdm_sym_faded, num, noise_var);
%         figure, subplot(1, 2, 1), mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft, abs(ch_tf_real)), subplot(1, 2, 2), mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft, abs(aaaaaa))
%         pause
        
        % 2d inverse sfft
        rx_sym_dd_ndft_eq_subfrm = sqrt(num.ndft/num.num_ofdmsym_per_subframe)*fft(ifft(rx_sym_tf_ndft_eq, [], 1), [], 2);
        
        % demap data qam symbols
        rx_sym_dd_ndft_data_subfrm = rx_sym_dd_ndft_eq_subfrm(:, num.idx_data_ofdmsym);
    end
    
    % buffer qam symbols
    rx_sym_dd_ndft(:, subfrm_idx) = rx_sym_dd_ndft_data_subfrm(:);
end

% demap the qam symbols
rx_bit_demod = (-1) * qamdemod(rx_sym_dd_ndft(:), 2^rm.Qm, 'UnitAveragePower', true, ...
    'OutputType', 'llr', 'NoiseVariance', noise_var);

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

end
