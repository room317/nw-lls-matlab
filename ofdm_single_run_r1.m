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

function pkt_error = ofdm_single_run_r1(sim, cc, rm, num, snr_db, ch, chest_option)

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
    %     'DirectPathInitialPhase', 0.5,...
    %     'PathGainsOutputPort', true);
else
    fading_ch = comm.RayleighChannel(...
        'SampleRate', num.sample_rate, ...
        'PathDelays', ch.path_delays, ...
        'AveragePathGains', ch.average_path_gains, ...
        'NormalizePathGains', true, ...
        'MaximumDopplerShift', ch.maximum_doppler_shift, ...
        'DopplerSpectrum', ch.doppler_spectrum);
    %     'PathGainsOutputPort', true, ...
    %     'RandomStream', 'mt19937ar with seed', ...
    %     'Seed', 72, ...
    %     'Visualization', 'Impulse response');   % enables the impulse response channel visualization
end

% calculate noise variance
noise_var = (num.ndft/num.nfft)*(10^((-0.1)*snr_db));

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

% modulate bit stream (only for 16qam)
tx_sym = qammod(tx_bit_ratematch, 2^rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);

% simulate per subframe
tx_sym_subfrm = reshape(tx_sym, num.len_rb_sym_user, []);
rx_sym_eq = zeros(num.ndft*num.num_data_ofdmsym_per_subframe, size(tx_sym_subfrm, 2));
for subfrm_idx = 1 : size(tx_sym_subfrm, 2)
    
    if strcmp(chest_option, 'real')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% pilot frame to get real channel %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % generate impulse train
        tx_sym_subfrm_reshape = ones(num.ndft, num.num_ofdmsym_per_subframe);
        
        % map qam symbols to resource block 
        tx_sym_map = zeros(num.nfft, num.num_ofdmsym_per_subframe);
        tx_sym_map((num.nfft/2)-(num.ndft/2)+1:(num.nfft/2)+(num.ndft/2), :) = tx_sym_subfrm_reshape;
        tx_sym_map_shift = fftshift(tx_sym_map, 1);
        
        % ofdm modulate
        tx_ofdm_sym = sqrt(num.nfft) * ifft(tx_sym_map_shift, [], 1);
        
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
        
        % ofdm demodulate the symbol
        rx_sym_map = (1/sqrt(num.nfft)) * fft(rx_ofdm_sym, [], 1);
        
        % demap qam symbols
        rx_sym_map_shift = fftshift(rx_sym_map, 1);
        rx_sym_demap = rx_sym_map_shift(num.nfft/2-num.ndft/2+1:num.nfft/2+num.ndft/2, :);
        
        % estimate and compensate channel
        [~, ch_real] = ofdm_ch_comp_tf_mmse_perfect(rx_sym_demap, tx_ofdm_sym_serial, tx_ofdm_sym_faded, num, noise_var, false, []);
    else
        ch_real = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%
    %%% normal frame %%%
    %%%%%%%%%%%%%%%%%%%%
    
    % generate pilot symbols
    tx_sym_subfrm_pilot = randi([0 1], num.ndft, num.num_pilot_ofdmsym_per_subframe)*2-1;
    
    % map qam symbols to resource block 
    tx_sym_subfrm_reshape = reshape(tx_sym_subfrm(:, subfrm_idx), num.ndft, num.num_data_ofdmsym_per_subframe);
    tx_sym_map = zeros(num.nfft, num.num_ofdmsym_per_subframe);
    tx_sym_map((num.nfft/2)-(num.ndft/2)+1:(num.nfft/2)+(num.ndft/2), num.idx_data_ofdmsym) = tx_sym_subfrm_reshape;   % data mapping
    tx_sym_map((num.nfft/2)-(num.ndft/2)+1:(num.nfft/2)+(num.ndft/2), num.idx_pilot_ofdmsym) = tx_sym_subfrm_pilot;   % pilot mapping
    tx_sym_map_shift = fftshift(tx_sym_map, 1);
    
    % ofdm modulate
    tx_ofdm_sym = sqrt(num.nfft) * ifft(tx_sym_map_shift, [], 1);
    
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
    % rx_ofdm_sym_serial = awgn(tx_ofdm_sym_serial, snr_db, 'measured');
    
    % reshape
    rx_ofdm_sym_cp = reshape(rx_ofdm_sym_serial, num.nfft+num.num_cp, num.num_ofdmsym_per_subframe);
    
    % remove cp
    rx_ofdm_sym = rx_ofdm_sym_cp(num.num_cp+1:end,:);
    
    % ofdm demodulate the symbol
    rx_sym_map = (1/sqrt(num.nfft)) * fft(rx_ofdm_sym, [], 1);
    
    % demap qam symbols
    rx_sym_map_shift = fftshift(rx_sym_map, 1);
    rx_sym_demap = rx_sym_map_shift(num.nfft/2-num.ndft/2+1:num.nfft/2+num.ndft/2, :);
    
    % estimate and compensate channel
    if strcmp(chest_option, 'tf_pilot')
        [rx_sym_eq_subfrm, ~] = ofdm_ch_comp_tf_mmse(rx_sym_demap, tx_sym_subfrm_pilot, num, noise_var);
    else
        [rx_sym_eq_subfrm, ~] = ofdm_ch_comp_tf_mmse_perfect(rx_sym_demap, tx_ofdm_sym_serial, tx_ofdm_sym_faded, num, noise_var, chest_option, ch_real);
    end
    
    % buffer qam symbols
    rx_sym_eq(:, subfrm_idx) = rx_sym_eq_subfrm(:);
end

% demap the qam symbols
rx_bit_demod = (-1) * qamdemod(rx_sym_eq(:), 2^rm.Qm, 'UnitAveragePower', true, ...
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
