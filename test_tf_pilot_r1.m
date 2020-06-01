% test dd pilot_r1 (real channel testing)
%   - test_metric: snr in db scale
%   - test_synch: synch position
%   - test_scope: print and plot results (no awgn)
%   - test_seed: random seed for channel
%   - test_chest: channel estimation option {'tf', 'perfect', 'real'}
%   - test_cheq: channel eq. option {'tf', 'tf_mmse'}
%   - test_dmrs: dmrs type {'lte_down', 'lte_up', 'nr'}
% created: 2020.03.06
% modified:
%   - 2020.03.02: prefix/postfix option 
%   - 2020.03.03: rician channel added
%   - 2020.03.06: pilot estimation added
%   - 2020.04.29: removed unused options
%   - 2020.05.19: tf pilot estimation added

function [qam_error, num_qam_per_pkt, ch_est_rmse, papr_db] = test_tf_pilot_r1(test_metric, test_synch, test_scope, test_seed, test_chest, test_cheq, test_dmrs)

% test
snr_db = test_metric;

% set parameter
num_ofdm_subc = 1024; % 64;
num_ofdm_sym = 14; % 64;
num_subc = 600; % 32;
num_sym = 14; % 32;
f_subc = 15e3;
f_s = f_subc*num_ofdm_subc;
qam_size = 16;
len_cp = 72; % floor(num_pilot_subc/2);
sym_avg_win_size = 1;   % symbol average window size for chest

% set channel
cfo_norm = 0; % normalized cfo. cfo = cfo_norm/t_sym;
carrier_freq_mhz = 4000;
velocity_kmh = 120;
idx_fading = 9; % 9: TDL-A, 12: TDL-D
delay_spread_rms_us = 0.01; %0.1e-6;
test_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, snr_db, delay_spread_rms_us);

% configure dmrs position
if (mod(num_subc, 12) == 0) && (num_sym == 14) && strcmp(test_dmrs, 'lte_down')
    idx_pilot_sym = [1 5 8 12];
    idx_pilot_subc = {(6:6:num_subc), (3:6:num_subc), (6:6:num_subc), (3:6:num_subc)};
elseif (mod(num_subc, 12) == 0) && (num_sym == 14) && strcmp(test_dmrs, 'nr')
    idx_pilot_sym = [3 6 9 12];
    idx_pilot_subc = {(2:2:num_subc), (2:2:num_subc), (2:2:num_subc), (2:2:num_subc)};
else
    idx_pilot_sym = 4:7:num_sym;
    idx_pilot_subc = repmat({(1:num_subc)}, 1, length(idx_pilot_sym));
end

% set pilot parameter
idx_sym = 1:num_sym;
idx_subc = 1:num_subc;
num_pilot_sym = length(idx_pilot_sym);
num_pilot_subc = length(idx_pilot_subc{1});

% create a rayleigh fading channel object
if test_seed >= 0
    rng(test_seed)
end
switch idx_fading
    case {12, 13}
        fading_ch = comm.RicianChannel(...
            'SampleRate', f_s,...
            'PathDelays', test_ch.path_delays,...
            'AveragePathGains', test_ch.average_path_gains,...
            'KFactor', test_ch.k_factor,...
            'DirectPathDopplerShift', test_ch.maximum_doppler_shift,...
            'MaximumDopplerShift', test_ch.maximum_doppler_shift,...
            'DopplerSpectrum', test_ch.doppler_spectrum,...
            'PathGainsOutputPort', true, ...
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

% calculate noise variance
noise_var = 10 ^ ((-0.1)*snr_db);

% calculate num. qam symbols per packet
num_qam_per_pkt = (num_subc*num_sym)-(num_pilot_subc*num_pilot_sym);

% generate bit stream
tx_bit = randi([0 qam_size-1], num_qam_per_pkt, 1);

% modulate bit stream
tx_sym_data = qammod(tx_bit, qam_size, 'UnitAveragePower', true);
fprintf('tx_sym_data  : %10.4f\n', mean(abs(tx_sym_data(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% generate pilot symbols (random sequence for now)
tx_sym_pilot_tf = 1/sqrt(2)*complex(randi([0 1], num_pilot_subc, num_pilot_sym)*2-1, randi([0 1], num_pilot_subc, num_pilot_sym)*2-1);

% reshape data
tx_sym_data1 = tx_sym_data(1:num_subc*(num_sym-num_pilot_sym));
tx_sym_data2 = tx_sym_data(num_subc*(num_sym-num_pilot_sym)+1:end);
tx_sym_data1_tf = reshape(tx_sym_data1, num_subc, num_sym-num_pilot_sym);
tx_sym_data2_tf = reshape(tx_sym_data2, num_subc-num_pilot_subc, num_pilot_sym);

% map data and pilot
idx_data_sym = idx_sym;
idx_data_sym(idx_pilot_sym) = [];
tx_sym_tf = zeros(num_subc, num_sym);
tx_sym_tf(:, idx_data_sym) = tx_sym_data1_tf;
tx_sym_pilot_ce_tf = zeros(num_subc, num_sym);  % for channel estimation
for i = 1:num_pilot_sym
    idx_data_subc = idx_subc;
    idx_data_subc(idx_pilot_subc{i}) = [];
    tx_sym_tf(idx_data_subc, idx_pilot_sym(i)) = tx_sym_data2_tf(:, i);
    tx_sym_tf(idx_pilot_subc{i}, idx_pilot_sym(i)) = tx_sym_pilot_tf(:, i);
    tx_sym_pilot_ce_tf(idx_pilot_subc{i}, idx_pilot_sym(i)) = tx_sym_pilot_tf(:, i);
end
fprintf('tx_sym_tf    : %10.4f\n', mean(abs(tx_sym_tf(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% map to tf domain
idx_ofdm_map = [floor((num_ofdm_subc-num_subc)/2), floor((num_ofdm_sym-num_sym)/2)];
tx_sym_map_shift_tf = zeros(num_ofdm_subc, num_ofdm_sym);
tx_sym_map_shift_tf(idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(2)+1:idx_ofdm_map(2)+num_sym) = tx_sym_tf;
tx_sym_map_tf = fftshift(tx_sym_map_shift_tf, 1);
fprintf('tx_sym_map_tf: %10.4f\n', mean(abs(tx_sym_map_tf(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% ofdm modulate
tx_ofdm_sym = sqrt(num_ofdm_subc) * ifft(tx_sym_map_tf, [], 1);
fprintf('tx_ofdm_sym  : %10.4f\n', mean(abs(tx_ofdm_sym(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% add cp (cyclic prefix)
tx_ofdm_sym_cp = [tx_ofdm_sym(end-len_cp+1:end, :); tx_ofdm_sym];       % cyclic-prefix

% serialize
tx_sig = tx_ofdm_sym_cp(:);
fprintf('tx_sig       : %10.4f\n', mean(abs(tx_sig(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% calculate papr
papr_db = 10*log10(max(abs(tx_sig).^2)/mean(abs(tx_sig).^2));

% pass signal through channel
[tx_sig_faded, ch_path_gain] = fading_ch(tx_sig);       % ch_path_gain: normally constant per path, vary when doppler exists
ch_info = info(fading_ch);
ch_filter_coeff = ch_info.ChannelFilterCoefficients;    % ch_filter_coeff: channel filter coefficient per path
num_path = length(fading_ch.PathDelays);                % num_path: number of path

assignin('base', 'ch_path_gain', ch_path_gain);
fprintf('tx_sig_faded : %10.4f\n', mean(abs(tx_sig_faded(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% reproduce real channel (time saving)
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
    path_gain_mat = ch_path_gain_reshape(:, idx_ch_sym, :);
    ch_mat_per_path = path_gain_mat .* coeff_mat;
    ch_real_mat_t(:, :, idx_ch_sym) = sum(ch_mat_per_path, 3);
    ch_real_halfmap_mat_tf(:, :, idx_ch_sym) = ifft(fft(ch_real_mat_t(:, :, idx_ch_sym), [], 1), [], 2);
end
ch_real_mat_shift_tf = fftshift(fftshift(ch_real_halfmap_mat_tf, 2), 1);
ch_real_mat_tf = ch_real_mat_shift_tf(idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, :);

% add gaussian noise
if test_scope
    rx_sig = tx_sig_faded;
else
    rx_sig = awgn(tx_sig_faded, snr_db, 'measured');
end
fprintf('rx_sig       : %10.4f\n', mean(abs(rx_sig(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% compensate cfo (to observe impact of frequency shift)
cfo_sample = 0:length(rx_sig)-1;
cfo_vec = cfo_sample.'*cfo_norm/(num_ofdm_subc+len_cp);
rx_sig_cfo = rx_sig(:) .* exp(-1i*2*pi*cfo_vec);

% synchronize
rx_sig_synch = circshift(rx_sig_cfo, -len_cp-test_synch);

% reshape
rx_ofdm_sym_cp = reshape(rx_sig_synch, num_ofdm_subc+len_cp, []);

% remove cp (with synch: to observe impact of time shift)
rx_ofdm_sym = rx_ofdm_sym_cp(1:num_ofdm_subc, :);
fprintf('rx_ofdm_sym  : %10.4f\n', mean(abs(rx_ofdm_sym(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% ofdm demodulate the symbol
rx_sym_map_tf = (1/sqrt(num_ofdm_subc)) * fft(rx_ofdm_sym, [], 1);
fprintf('rx_sym_map_tf: %10.4f\n', mean(abs(rx_sym_map_tf(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% demap
rx_sym_map_shift_tf = fftshift(rx_sym_map_tf, 1);
rx_sym_tf = rx_sym_map_shift_tf(idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(2)+1:idx_ofdm_map(2)+num_sym);
fprintf('rx_sym_tf    : %10.4f\n', mean(abs(rx_sym_tf(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% demap pilots and average pilots
rx_sym_pilot_tf = zeros(num_pilot_subc, num_pilot_sym);
rx_sym_pilot_map_tf = zeros(num_subc, num_sym);     % for channel estimation
rx_sym_pilot_avg_tf = zeros(num_subc, num_sym);     % for channel estimation
rx_sym_pilot_interp_sym_tf = zeros(num_subc, num_sym);  % for channel estimation
rx_sym_pilot_interp_subc_tf = zeros(num_subc, num_sym);  % for channel estimation
for i = 1:num_pilot_sym
    % demap
    rx_sym_pilot_tf(:, i) = rx_sym_tf(idx_pilot_subc{i}, idx_pilot_sym(i));
    
    % compensate pilots
    rx_sym_pilot_ch_tf = rx_sym_pilot_tf(:, i) .* conj(tx_sym_pilot_tf(:, i));
    
    % map pilot channels
    rx_sym_pilot_map_tf(idx_pilot_subc{i}, idx_pilot_sym(i)) = rx_sym_pilot_ch_tf;
    
    % find symbol index within averaging window
    idx_pilot_avg = idx_pilot_sym(idx_pilot_sym>=max(1, idx_pilot_sym(i)-sym_avg_win_size) & idx_pilot_sym<=idx_pilot_sym(i));
    len_pilot_avg = length(idx_pilot_avg);
    
    % initialize
    list_pilot_avg = zeros(num_subc, num_pilot_sym);
    cnt_pilot_avg = zeros(num_pilot_subc, 1);
    
    % average pilots across symbols
    for j = i-len_pilot_avg+1:i
        list_pilot_avg(idx_pilot_subc{j}, j) = rx_sym_pilot_map_tf(idx_pilot_subc{j}, idx_pilot_sym(j));
        cnt_pilot_avg = cnt_pilot_avg+double(ismember(idx_pilot_subc{j}, idx_pilot_subc{i})).';
    end
    sum_pilot_avg = sum(list_pilot_avg, 2);
    rx_sym_pilot_avg_tf(idx_pilot_subc{i}, idx_pilot_sym(i)) = sum_pilot_avg(idx_pilot_subc{i}, 1)./cnt_pilot_avg;
    
    % generate virtual pilots (ignored)
    % interpolate 1st dimension
    rx_sym_pilot_interp_sym_tf(:, idx_pilot_sym(i)) = ...
        interp1(idx_pilot_subc{i}, rx_sym_pilot_avg_tf(idx_pilot_subc{i}, idx_pilot_sym(i)), idx_subc, 'linear', 'extrap');
end

% interpolate 2nd dimension
for i = 1:num_subc
    rx_sym_pilot_interp_subc_tf(i, :) = ...
        interp1(idx_pilot_sym, rx_sym_pilot_interp_sym_tf(i, idx_pilot_sym), idx_sym, 'linear', 'extrap');
end

% extract diagonal elements of real channel
ch_est_tf_real = zeros(num_subc, num_sym);
for idx_sym = 1:num_sym
    ch_est_tf_real(:, idx_sym) = diag(ch_real_mat_tf(:, :, idx_sym));  % diagonal term only
end

% estimate channel
if strcmp(test_chest, 'tf')
    % estimate tf channel
    ch_est_tf = rx_sym_pilot_interp_subc_tf;
elseif strcmp(test_chest, 'perfect')
    % estimate tf channel
    ch_est_tf_perfect = rx_sym_tf ./ tx_sym_tf;
    ch_est_tf = ch_est_tf_perfect;
elseif strcmp(test_chest, 'real')
    % get real channel
    ch_est_tf = ch_est_tf_real;
else
    error('test_chest value must be one of these: {''tf'', ''perfect'', ''real''}')
end

% calculate channel rmse
ch_est_rmse = sqrt(mean(abs(ch_est_tf_real-ch_est_tf).^2, 'all'));  % diagonal term only

% check channel estimation
if test_scope
    ch_real_mat_diag_tf = zeros(num_subc, num_sym);
    for idx_sym = 1 : num_sym
        ch_real_mat_diag_tf(:, idx_sym) = diag(ch_real_mat_tf(:, :, idx_sym));
    end
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(ch_est_tf)), test_axis = axis;
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(ch_real_mat_diag_tf)), axis(test_axis)
    pause
%     for idx_sym = 1 : num_sym
%         figure(10)
%         subplot(1, 2, 1), mesh(1:num_subc, 1:num_subc, abs(diag(ch_est_tf(:, idx_sym)))), test_axis = axis;
%         subplot(1, 2, 2), mesh(1:num_subc, 1:num_subc, abs(ch_real_mat_tf(:, :, idx_sym))), axis(test_axis)
%         pause
%     end
end

% equalize channel
if strcmp(test_cheq, 'tf')
    % equalize channel in tf domain
    rx_sym_tf_eq = rx_sym_tf ./ ch_est_tf;
elseif strcmp(test_cheq, 'tf_mmse')
    % equalize channel in tf domain
    ch_est_tf_mmse = conj(ch_est_tf) ./ (noise_var+abs(ch_est_tf).^2);
    rx_sym_tf_eq = rx_sym_tf .* ch_est_tf_mmse;
else
    error('test_chest value must be one of these: {''tf'', ''tf_mmse''}')
end
fprintf('rx_sym_tf_eq : %10.4f\n', mean(abs(rx_sym_tf_eq(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% demap data
rx_sym_data1_tf = rx_sym_tf_eq(:, idx_data_sym);
rx_sym_data2_tf = zeros(num_subc-num_pilot_subc, num_pilot_sym);
for i = 1:num_pilot_sym
    % demap
    idx_data_subc = idx_subc;
    idx_data_subc(idx_pilot_subc{i}) = [];
    rx_sym_data2_tf(:, i) = rx_sym_tf_eq(idx_data_subc, idx_pilot_sym(i));
end

% reshape data
rx_sym_data = [rx_sym_data1_tf(:); rx_sym_data2_tf(:)];
fprintf('rx_sym_data  : %10.4f\n', mean(abs(rx_sym_data(:)).^2))    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% demodulate qam symbols
rx_bit = qamdemod(rx_sym_data, qam_size, 'UnitAveragePower', true);

% calculate bit errors
if isempty(tx_bit)
    qam_error = 0;
else
    qam_error = symerr(tx_bit, rx_bit);
end

% plot
if test_scope
    figure, plot(tx_bit, '-b.'), hold on, plot(rx_bit, ':r.'), hold off, grid minor
    figure
    subplot(2, 1, 1), plot(real(tx_sym_data(:)), '-b.'), hold on, plot(real(rx_sym_data(:)), ':r.'), hold off, grid minor
    subplot(2, 1, 2), plot(imag(tx_sym_data(:)), '-b.'), hold on, plot(imag(rx_sym_data(:)), ':r.'), hold off, grid minor
    figure, scatter(real(rx_sym_data(:)), imag(rx_sym_data(:)), 'r.'), hold on, scatter(real(tx_sym_data(:)), imag(tx_sym_data(:)), 'ko'), hold off, axis([-2 2 -2 2])
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_tf)), test_axis = axis; subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_tf_eq)), axis(test_axis)
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_tf)), test_axis = axis; subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_tf)), axis(test_axis)
    figure, subplot(1, 2, 1), mesh(1:num_ofdm_sym, 1:num_ofdm_subc, abs(tx_sym_map_tf)), subplot(1, 2, 2), mesh(1:num_ofdm_sym, 1:num_ofdm_subc, abs(rx_sym_map_tf))
    figure, subplot(1, 2, 1), mesh(1:num_ofdm_sym, 1:num_ofdm_subc, abs(tx_ofdm_sym)), subplot(1, 2, 2), mesh(1:num_ofdm_sym, 1:num_ofdm_subc, abs(rx_ofdm_sym))
    figure, subplot(1, 2, 1), mesh(1:num_ofdm_sym, 1:num_ofdm_subc+len_cp, abs(tx_ofdm_sym_cp)), subplot(1, 2, 2), mesh(1:num_ofdm_sym, 1:num_ofdm_subc+len_cp, abs(rx_ofdm_sym_cp))
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(ch_est_tf_real)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(ch_est_tf))
    
    % channel test
    test_rx_sym_tf_est = tx_sym_tf .* ch_est_tf;
    figure
    subplot(2, 1, 1), plot(real(rx_sym_tf(:)), '-b.'), hold on, plot(real(test_rx_sym_tf_est(:)), ':r.'), hold off, title('real(rx sig) and real(rx sig est) in tf-domain'), grid minor
    subplot(2, 1, 2), plot(imag(rx_sym_tf(:)), '-b.'), hold on, plot(imag(test_rx_sym_tf_est(:)), ':r.'), hold off, title('imag(rx sig) and imag(rx sig est) in tf-domain'), grid minor
end

% dump
if test_scope
    assignin('base', 'tx_sym_data', tx_sym_data);
    assignin('base', 'tx_sym_tf', tx_sym_tf);
    assignin('base', 'tx_sym_map_tf', tx_sym_map_tf);
    assignin('base', 'tx_ofdm_sym', tx_ofdm_sym);
    assignin('base', 'tx_ofdm_sym_cp', tx_ofdm_sym_cp);
    assignin('base', 'tx_sig', tx_sig);
    assignin('base', 'tx_sig_faded', tx_sig_faded);
    assignin('base', 'rx_sig', rx_sig);
    assignin('base', 'rx_ofdm_sym_cp', rx_ofdm_sym_cp);
    assignin('base', 'rx_ofdm_sym', rx_ofdm_sym);
    assignin('base', 'rx_sym_map_tf', rx_sym_map_tf);
    assignin('base', 'rx_sym_map_shift_tf', rx_sym_map_shift_tf);
    assignin('base', 'rx_sym_tf', rx_sym_tf);
    assignin('base', 'rx_sym_tf_eq', rx_sym_tf_eq);
    assignin('base', 'rx_sym_data', rx_sym_data);
    assignin('base', 'ch_est_tf', ch_est_tf);
    assignin('base', 'test_ch', test_ch);
end

end
