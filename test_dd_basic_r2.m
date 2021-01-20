% test dd basic_r1 (real channel testing)
% to observe real channel reproduction in dd domain
%   - test_synch: synch position for test
%   - test_singletone: to check channel change
%   - test_scope: print and plot results
%   - test_seed: synch position
%   - test_timesaving: save time to regenerate channels
%   - test_latency: latency test with partial reception (number of symbols received)
% created: 2020.03.23
% modified:
%   - 2020.03.28:
%   - 2020.04.29:
%   - 2020.06.30: walsh_hadamard spreading (not working as intended)
%   - 2020.06.30: papr calculation
%   - 2020.06.30: symbol scrambling (not working as intended)
%   - 2020.07.19: latency test with partial reception (test_latency)
% memo
%   - https://kr.mathworks.com/help/comm/ref/comm.rayleighchannel-system-object.html#d120e191883

function [ch_rmse, papr_db] = test_dd_basic_r2(test_synch, test_tone, test_scope, test_seed, test_timesaving, test_latency)

% set parameter
num_ofdm_subc = 1024;
num_ofdm_sym = 14;
num_subc = 600; % 32;
num_sym = 14; % 32;
qam_size = 16;
len_cp = 72; % num_fft/4; % 0;
f_subc = 15e3;
f_s = f_subc*num_ofdm_subc;
% t_sym = (num_ofdm_subc+len_cp)/f_s;
cfo_norm = 0; % normalized cfo. cfo = cfo_norm/t_sym;
carrier_freq_mhz = 4000;
velocity_kmh = 3; % 120;
idx_fading = 9; % 9: TDL-A, 12: TDL-D
delay_spread_rms_us = 0.1; % 0.1; %0.1e-6;
test_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, delay_spread_rms_us);

% % check error
% if test_synch > len_cp
%     error('test_synch should not be greater than %d.\n', len_cp)
% end

if test_seed >= 0
    rng(test_seed)
else
    test_seed = rng;
end

% create a rayleigh fading channel object
rng(test_seed)
switch idx_fading
    case {12, 13}
        fading_ch = comm.RicianChannel(...
            'SampleRate', f_s,...
            'PathDelays', test_ch.path_delays,...
            'AveragePathGains', test_ch.average_path_gains,...
            'KFactor', test_ch.k_factor,...
            'DirectPathDopplerShift', test_ch.maximum_doppler_shift,...
            'MaximumDopplerShift', test_ch.maximum_doppler_shift,...
            'PathGainsOutputPort', true, ...
            'DopplerSpectrum', test_ch.doppler_spectrum,...
            'NormalizePathGains', true);
    otherwise
        fading_ch = comm.RayleighChannel(...
            'SampleRate', f_s, ...
            'PathDelays', test_ch.path_delays, ...
            'AveragePathGains', test_ch.average_path_gains, ...
            'NormalizePathGains', true, ...
            'MaximumDopplerShift', test_ch.maximum_doppler_shift, ...
            'PathGainsOutputPort', true, ...
            'DopplerSpectrum', test_ch.doppler_spectrum);
end

if isempty(test_tone)
    % generate bit stream
    tx_bit = randi([0 qam_size-1], num_subc*num_sym, 1);
    
    % modulate bit stream
    tx_sym = qammod(tx_bit, qam_size, 'UnitAveragePower', true);
    
    % reshape
    tx_sym_dd = reshape(tx_sym, num_subc, []);
else
    % dd tone generation
    tx_sym_dd = zeros(num_subc, num_sym);
    tx_sym_dd(test_tone(1)+1, test_tone(2)+1) = sqrt(num_subc*num_sym);
end

% % scramble data
% tx_sym_dd(:, test_tone(2)+1) = lfsr_soft(tx_sym_dd(:, test_tone(2)+1), [1 1 1 1 1 1 1]);

% tx_sym_dd = zeros(num_subc, num_sym);
% tx_sym_dd(tx_pilot_position(1)+1, tx_pilot_position(2)+1) = sqrt(num_subc*num_sym);

% 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
tx_sym_tf = sqrt(num_sym/num_subc) * fft(ifft(tx_sym_dd, [], 2), [], 1);

% map to tf domain
idx_ofdm_map = [num_ofdm_subc/2-num_subc/2, num_ofdm_sym/2-num_sym/2];
tx_sym_map_shift_tf = zeros(num_ofdm_subc, num_ofdm_sym);
tx_sym_map_shift_tf(idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(2)+1:idx_ofdm_map(2)+num_sym) = tx_sym_tf;
tx_sym_map_tf = fftshift(tx_sym_map_shift_tf, 1);

% ofdm modulate
tx_ofdm_sym = sqrt(num_ofdm_subc) * ifft(tx_sym_map_tf, [], 1);

% spread with walsh-hadamard
% tx_ofdm_sym_wh = fwht(tx_ofdm_sym);
tx_ofdm_sym_wh = tx_ofdm_sym;

% add cp (cyclic prefix)
tx_ofdm_sym_cp = [tx_ofdm_sym_wh(end-len_cp+1:end, :); tx_ofdm_sym_wh];

% serialize
tx_sig = tx_ofdm_sym_cp(:);

% pass signal through channel
rng(test_seed)
[tx_sig_faded, ch_path_gain] = fading_ch(tx_sig);       % ch_path_gain: normally constant per path, vary when doppler exists
ch_info = info(fading_ch);
ch_filter_coeff = ch_info.ChannelFilterCoefficients;    % ch_filter_coeff: channel filter coefficient per path
num_path = length(fading_ch.PathDelays);                % num_path: number of path

% reproduce real channel
if test_timesaving == 1
    % time saving (matrix operation + for loop: a little faster than full matrix operation)
    len_ch_filter = size(ch_filter_coeff, 2);
    ch_real_mat_t = zeros(num_ofdm_subc, num_ofdm_subc, num_sym);
    ch_real_halfmap_mat_tf = zeros(num_ofdm_subc, num_ofdm_subc, num_sym);
    coeff_mat = zeros(num_ofdm_subc, num_ofdm_subc, num_path);
    ch_path_gain_reshape = reshape(ch_path_gain, num_ofdm_subc+len_cp, num_ofdm_sym, num_path);
    ch_path_gain_reshape(1:len_cp, :, :) = [];
    ch_path_gain_reshape = ch_path_gain_reshape(:, idx_ofdm_map(2)+1:idx_ofdm_map(2)+num_sym, :);
    for idx_ch_path = 1:num_path
        coeff_row = circshift([ch_filter_coeff(idx_ch_path, end:-1:1) zeros(1, num_ofdm_subc - len_ch_filter)], -len_ch_filter+1);
        coeff_col = [ch_filter_coeff(idx_ch_path, :) zeros(1, num_ofdm_subc - len_ch_filter)];
        coeff_mat(:, :, idx_ch_path) = toeplitz(coeff_col, coeff_row);
    end
    for idx_ch_sym = 1:num_sym
%         path_gain_mat = repmat(ch_path_gain_reshape(:, idx_ch_sym, :), 1, num_ofdm_subc, 1);
        path_gain_mat = ch_path_gain_reshape(:, idx_ch_sym, :);
        ch_mat_per_path = path_gain_mat .* coeff_mat;
        ch_real_mat_t(:, :, idx_ch_sym) = sum(ch_mat_per_path, 3);
        ch_real_halfmap_mat_tf(:, :, idx_ch_sym) = ifft(fft(ch_real_mat_t(:, :, idx_ch_sym), [], 1), [], 2);
    end
    ch_real_mat_shift_tf = fftshift(fftshift(ch_real_halfmap_mat_tf, 2), 1);
    ch_real_mat_tf = ch_real_mat_shift_tf(idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, :);
elseif test_timesaving == 2
    % time saving (full matrix operation)
    len_ch_filter = size(ch_filter_coeff, 2);
    coeff_mat = zeros(num_ofdm_subc, num_ofdm_subc, 1, num_path);
    ch_path_gain_reshape = reshape(ch_path_gain, num_ofdm_subc+len_cp, 1, num_ofdm_sym, num_path);
    ch_path_gain_reshape(1:len_cp, :, :, :) = [];
    ch_path_gain_reshape = ch_path_gain_reshape(:, :, idx_ofdm_map(2)+1:idx_ofdm_map(2)+num_sym, :);
    for idx_ch_path = 1:num_path
        coeff_row = circshift([ch_filter_coeff(idx_ch_path, end:-1:1) zeros(1, num_ofdm_subc - len_ch_filter)], -len_ch_filter+1);
        coeff_col = [ch_filter_coeff(idx_ch_path, :) zeros(1, num_ofdm_subc - len_ch_filter)];
        coeff_mat(:, :, 1, idx_ch_path) = toeplitz(coeff_col, coeff_row);
    end
%     path_gain_mat = repmat(ch_path_gain_reshape, 1, num_ofdm_subc, 1, 1);
%     coeff_total_mat = repmat(coeff_mat, 1, 1, num_sym, 1);
    path_gain_mat = ch_path_gain_reshape;
    coeff_total_mat = coeff_mat;
    ch_mat_per_path = path_gain_mat .* coeff_total_mat;
    ch_real_mat_t = squeeze(sum(ch_mat_per_path, 4));
    ch_real_halfmap_mat_tf = ifft(fft(ch_real_mat_t, [], 1), [], 2);
    ch_real_mat_shift_tf = fftshift(fftshift(ch_real_halfmap_mat_tf, 2), 1);
    ch_real_mat_tf = ch_real_mat_shift_tf(idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, :);
else
    % original
    len_ch_filter = size(ch_filter_coeff, 2);
    ch_real_map_mat_t = zeros(num_ofdm_subc, num_ofdm_subc, num_ofdm_sym);
    ch_real_map_mat_tf = zeros(num_ofdm_subc, num_ofdm_subc, num_ofdm_sym);
    for idx_ch_sym = 1:num_ofdm_sym
        ch_mat_per_path = zeros(num_ofdm_subc, num_ofdm_subc, num_path);
        for idx_ch_path = 1:num_path
            coeff_row = circshift([ch_filter_coeff(idx_ch_path, end:-1:1) zeros(1, num_ofdm_subc - len_ch_filter)], -len_ch_filter+1);
            coeff_col = [ch_filter_coeff(idx_ch_path, :) zeros(1, num_ofdm_subc - len_ch_filter)];
            coeff_mat = toeplitz(coeff_col, coeff_row);
            path_gain_mat = diag(ch_path_gain((num_ofdm_subc+len_cp)*(idx_ch_sym-1)+len_cp+1:(num_ofdm_subc+len_cp)*idx_ch_sym, idx_ch_path));
            ch_mat_per_path(:, :, idx_ch_path) = path_gain_mat * coeff_mat;
        end
        ch_real_map_mat_t(:, :, idx_ch_sym) = sum(ch_mat_per_path, 3);
        ch_real_map_mat_tf(:, :, idx_ch_sym) = dftmtx(num_ofdm_subc)*ch_real_map_mat_t(:, :, idx_ch_sym)*conj(dftmtx(num_ofdm_subc))/num_ofdm_subc;
    end
    ch_real_mat_shift_tf = fftshift(fftshift(ch_real_map_mat_tf, 2), 1);
    ch_real_mat_tf = ch_real_mat_shift_tf(idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(2)+1:idx_ofdm_map(2)+num_sym);
end

% add gaussian noise
rx_sig = tx_sig_faded;

% cfo
cfo_sample = 0:length(rx_sig)-1;
cfo_vec = cfo_sample.'*cfo_norm/(num_ofdm_subc+len_cp);
rx_sig_cfo = rx_sig(:) .* exp(-1i*2*pi*cfo_vec);

% synchronization
rx_sig_synch = circshift(rx_sig_cfo, -len_cp-test_synch);

% reshape
rx_ofdm_sym_cp = reshape(rx_sig_synch, num_ofdm_subc+len_cp, []);

% remove cp (with synch: to observe impact of time shift)
rx_ofdm_sym = rx_ofdm_sym_cp(1:num_ofdm_subc, :);

% remove parts of received symbols
if isscalar(test_latency) && ismember(test_latency, 1:num_ofdm_sym-1)
    rx_ofdm_sym_part = [rx_ofdm_sym(:, 1:test_latency) zeros(size(rx_ofdm_sym, 1), num_ofdm_sym-test_latency)];
else
    rx_ofdm_sym_part = rx_ofdm_sym;
end

% despread with walsh-hadamard
% rx_ofdm_sym_wh = ifwht(rx_ofdm_sym);
rx_ofdm_sym_wh = rx_ofdm_sym_part;

% ofdm demodulate the symbol
rx_sym_map_tf = (1/sqrt(num_ofdm_subc)) * fft(rx_ofdm_sym_wh, [], 1);

% demap
rx_sym_map_shift_tf = fftshift(rx_sym_map_tf, 1);
rx_sym_tf = rx_sym_map_shift_tf(idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(2)+1:idx_ofdm_map(2)+num_sym);

% 2d inverse sfft
rx_sym_dd = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf, [], 1), [], 2);

% % descramble data
% for idx = 1 : num_sym
%     rx_sym_dd(:, idx) = lfsr_soft(rx_sym_dd(:, idx), [1 1 1 1 1 1 1]);
% end

% output: channel estimation error
if isempty(test_tone)
    ch_est_tf_perfect = rx_sym_tf ./ tx_sym_tf;     % tf domain
    ch_est_tf = ch_est_tf_perfect;
else
    ch_est_dd = circshift(rx_sym_dd, (-1)*test_tone);
    ch_est_tf = sqrt(num_sym/num_subc) * fft(ifft(ch_est_dd, [], 2), [], 1);
end
ch_mse = zeros(1, num_sym);
for idx_ce_sym = 1:num_sym
%     ce_est_mse(idx_ce_sym) = mean(abs(ch_real_mat_tf(:, :, idx_ce_sym)-diag(ch_est_tf(:, idx_ce_sym))).^2, 'all');  % off-diagonal term included
    ch_mse(idx_ce_sym) = mean(abs(diag(ch_real_mat_tf(:, :, idx_ce_sym))-ch_est_tf(:, idx_ce_sym)).^2, 'all');  % diagonal term only
end
ch_rmse = sqrt(mean(ch_mse));

% output: papr
papr_db = 10*log10(mean(max(abs(tx_ofdm_sym_cp).^2, [], 1)./mean(abs(tx_ofdm_sym_cp).^2, 1)));

if test_scope
    % reproduce channel output
    ch_in = tx_sig;
    ch_delayed = zeros(size(ch_in, 1), num_path);
    for idx_path = 1:num_path
        ch_delayed(:, idx_path) = filter(ch_filter_coeff(idx_path, :), 1, ch_in);
    end
    ch_out = sum(ch_path_gain .* ch_delayed, 2);
    tx_sig_faded_rep = ch_out;
    
    figure
    subplot(2, 1, 1), plot(real(tx_sig_faded), '-b.'), hold on, plot(real(tx_sig_faded_rep), ':r.'), hold off, grid minor
    subplot(2, 1, 2), plot(imag(tx_sig_faded), '-b.'), hold on, plot(imag(tx_sig_faded_rep), ':r.'), hold off, grid minor
    pause
    
    % test real channel
    if ~test_timesaving
        % reproduce time domain rx signals
        test_rx_ofdm_sym_est = zeros(size(rx_ofdm_sym));
        for idx_sym = 1 : num_ofdm_sym
            test_rx_ofdm_sym_est(:, idx_sym) = ch_real_map_mat_t(:, :, idx_sym) * tx_ofdm_sym(:, idx_sym);
        end
        figure
        subplot(2, 1, 1), plot(real(rx_ofdm_sym(:)), '-b.'), hold on, plot(real(test_rx_ofdm_sym_est(:)), ':r.'), hold off, grid minor
        subplot(2, 1, 2), plot(imag(rx_ofdm_sym(:)), '-b.'), hold on, plot(imag(test_rx_ofdm_sym_est(:)), ':r.'), hold off, grid minor
        pause
        
        % reproduce frequency domain rx signals
        test_rx_sym_map_tf_est = zeros(size(rx_sym_map_tf));
        for idx_sym = 1 : num_ofdm_sym
            test_rx_sym_map_tf_est(:, idx_sym) = ch_real_map_mat_tf(:, :, idx_sym) * tx_sym_map_tf(:, idx_sym);
        end
        figure
        subplot(2, 1, 1), plot(real(rx_sym_map_tf(:)), '-b.'), hold on, plot(real(test_rx_sym_map_tf_est(:)), ':r.'), hold off, grid minor
        subplot(2, 1, 2), plot(imag(rx_sym_map_tf(:)), '-b.'), hold on, plot(imag(test_rx_sym_map_tf_est(:)), ':r.'), hold off, grid minor
        pause
        
        % reproduce frequency domain rx signals
        test_rx_sym_tf_est = zeros(size(rx_sym_tf));
        for idx_sym = 1 : num_sym
            test_rx_sym_tf_est(:, idx_sym) = ch_real_mat_tf(:, :, idx_sym) * tx_sym_tf(:, idx_sym);
        end
        figure
        subplot(2, 1, 1), plot(real(rx_sym_tf(:)), '-b.'), hold on, plot(real(test_rx_sym_tf_est(:)), ':r.'), hold off, grid minor
        subplot(2, 1, 2), plot(imag(rx_sym_tf(:)), '-b.'), hold on, plot(imag(test_rx_sym_tf_est(:)), ':r.'), hold off, grid minor
        pause
    end
    
    % reproduce frequency domain rx signals
    test_rx_sym_tf_est = zeros(size(rx_sym_tf));
    for idx_sym = 1 : num_sym
        test_rx_sym_tf_est(:, idx_sym) = ch_real_mat_tf(:, :, idx_sym) * tx_sym_tf(:, idx_sym);
    end
    figure
    subplot(2, 1, 1), plot(real(rx_sym_tf(:)), '-b.'), hold on, plot(real(test_rx_sym_tf_est(:)), ':r.'), hold off, grid minor
    subplot(2, 1, 2), plot(imag(rx_sym_tf(:)), '-b.'), hold on, plot(imag(test_rx_sym_tf_est(:)), ':r.'), hold off, grid minor
    pause
    
    % compare tx and rx signals
    ch_real_tf = zeros(num_subc, num_sym);
    for idx_ce_sym = 1:num_sym
        ch_real_tf(:, idx_ce_sym) = diag(ch_real_mat_tf(:, :, idx_ce_sym));  % diagonal term only
    end
    ch_real_dd = circshift(sqrt(num_subc/num_sym) * fft(ifft(ch_real_tf, [], 1), [], 2), test_tone);
    figure
    subplot(2, 3, 1), mesh(1:size(tx_sym_dd, 2), 1:size(tx_sym_dd, 1), real(tx_sym_dd))
    subplot(2, 3, 4), mesh(1:size(tx_sym_dd, 2), 1:size(tx_sym_dd, 1), imag(tx_sym_dd))
    subplot(2, 3, 2), mesh(1:size(rx_sym_dd, 2), 1:size(rx_sym_dd, 1), real(rx_sym_dd))
    subplot(2, 3, 5), mesh(1:size(rx_sym_dd, 2), 1:size(rx_sym_dd, 1), imag(rx_sym_dd))
    subplot(2, 3, 3), mesh(1:size(ch_real_dd, 2), 1:size(ch_real_dd, 1), real(ch_real_dd))
    subplot(2, 3, 6), mesh(1:size(ch_real_dd, 2), 1:size(ch_real_dd, 1), imag(ch_real_dd))
    pause
    
    % check channel estimation
    for idx_sym = 1 : num_sym
        figure(10)
        subplot(1, 2, 1), mesh(1:num_subc, 1:num_subc, abs(diag(ch_est_tf(:, idx_sym)))), test_axis = axis;
        subplot(1, 2, 2), mesh(1:num_subc, 1:num_subc, abs(ch_real_mat_tf(:, :, idx_sym))), axis(test_axis)
        pause
    end
end

% dump
if test_scope
    assignin('base', 'tx_sym_dd', tx_sym_dd);
    assignin('base', 'tx_sym_tf', tx_sym_tf);
    assignin('base', 'tx_sym_map_tf', tx_sym_map_tf);
    assignin('base', 'tx_ofdm_sym', tx_ofdm_sym);
    assignin('base', 'tx_ofdm_sym_cp', tx_ofdm_sym_cp);
    assignin('base', 'tx_sig', tx_sig);
    assignin('base', 'tx_sig_faded', tx_sig_faded);
    assignin('base', 'rx_ofdm_sym_cp', rx_ofdm_sym_cp);
    assignin('base', 'rx_ofdm_sym', rx_ofdm_sym);
    assignin('base', 'rx_sym_map_tf', rx_sym_map_tf);
    assignin('base', 'rx_sym_tf', rx_sym_tf);
    assignin('base', 'rx_sym_dd', rx_sym_dd);
end

end
