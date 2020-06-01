% test dd pilot_r3 (real channel testing)
%   - snr_db: snr in db scale
%   - test_synch: synch position
%   - test_scope: print and plot results
%   - test_seed: random seed for channel
%   - test_chest: channel estimation option {'dd', 'perfect', 'real'}
%   - test_cheq: channel eq. option {'tf', 'tf_mmse', 'dd'}
% created: 2020.03.06
% modified:
%   - 2020.03.02: prefix/postfix option
%   - 2020.03.03: rician channel added
%   - 2020.03.06: pilot estimation added
%   - 2020.04.29: removed unused options

function [qam_error, num_qam_per_pkt, ch_est_rmse, papr_db] = test_dd_pilot_r3(test_metric, test_synch, test_scope, test_seed, test_chest, test_cheq)

% test
snr_db = test_metric;

% set parameter
num_ofdm_subc = 1024; % 64;
num_ofdm_sym = 14; % 64;
num_subc = 600; % 32;
num_sym = 14; % 32;
num_pilot_subc = 84; % delay spread axis
num_pilot_sym = 14;

f_subc = 15e3;
f_s = f_subc*num_ofdm_subc;
qam_size = 16;
len_cp = 72; % floor(num_pilot_subc/2);
cfo_norm = 0; % normalized cfo. cfo = cfo_norm/t_sym;
carrier_freq_mhz = 4000;
velocity_kmh = 120;
idx_fading = 9; % 9: TDL-A, 12: TDL-D
delay_spread_rms_us = 0.01; %0.1e-6;
test_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, snr_db, delay_spread_rms_us);

% guard area
num_pilot_guard_subc = 14; % 20;
num_pilot_guard_sym = 0;

% null area
num_data_null_subc_head = 5; % 1;
num_data_null_subc_tail = 5; % 5;

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

% calculate noise variance
noise_var = 10 ^ ((-0.1)*snr_db);

if strcmp(test_chest, 'dd')
    % calculate num. qam symbols per packet
    num_qam_per_pkt = ...
        (num_subc*num_sym)- ...
        (num_pilot_subc*num_pilot_sym)- ...
        ((num_data_null_subc_head+num_data_null_subc_tail)*num_pilot_sym);
    
    % generate bit stream
    tx_bit = randi([0 qam_size-1], num_qam_per_pkt, 1);
    
    % modulate bit stream
    tx_sym_data = qammod(tx_bit, qam_size, 'UnitAveragePower', true);
    
%     % help map
%     nr=16;
%     np=5; p=zeros(np, 1); p(ceil(np/2)+1)=1;
%     nd1=floor((nr-np)/2); d1=ones(nd1,1);
%     nd2=ceil((nr-np)/2); d2=ones(nd2,1);
%     r=[d1; p; d2];
%     rr=zeros(nr,1); rr(floor(nr/2)+1)=1; % for comparison
    
    % reshape data symbol stream
    num_tx_sym_data1 = [floor((num_subc-num_pilot_subc)/2)-num_data_null_subc_head num_sym];
    num_tx_sym_data2 = [num_pilot_subc num_sym-num_pilot_sym];
    num_tx_sym_data3 = [ceil((num_subc-num_pilot_subc)/2)-num_data_null_subc_tail num_sym];
    tx_sym_data1 = tx_sym_data(1:prod(num_tx_sym_data1));
    tx_sym_data2 = tx_sym_data(prod(num_tx_sym_data1)+1:prod(num_tx_sym_data1)+prod(num_tx_sym_data2));
    tx_sym_data3 = tx_sym_data(prod(num_tx_sym_data1)+prod(num_tx_sym_data2)+1:end);
    tx_sym_data1_dd = reshape(tx_sym_data1, num_tx_sym_data1(1), []);
    tx_sym_data2_dd = reshape(tx_sym_data2, num_tx_sym_data2(1), []);
    tx_sym_data3_dd = reshape(tx_sym_data3, num_tx_sym_data3(1), []);
    
    % generate pilot symbols with guard
    tx_sym_pilot_dd = zeros(num_pilot_subc, num_pilot_sym);
    tx_sym_pilot_dd(floor(num_pilot_subc/2)+1, floor(num_pilot_sym/2)+1) = sqrt(num_pilot_subc*num_pilot_sym);
    
    % generate null symbols
    tx_sym_null_head_dd = zeros(num_data_null_subc_head, num_sym);
    tx_sym_null_tail_dd = zeros(num_data_null_subc_tail, num_sym);
    
    % map data and pilot
    tx_sym_dd_circshift = [...
        tx_sym_null_head_dd;
        tx_sym_data1_dd;
        tx_sym_pilot_dd, tx_sym_data2_dd;
        tx_sym_data3_dd;
        tx_sym_null_tail_dd];
    tx_sym_dd = circshift(tx_sym_dd_circshift, floor((num_sym-num_pilot_sym)/2), 2);
else
    % calculate num. qam symbols per packet
    num_qam_per_pkt = num_subc*num_sym;
    
    % generate bit stream
    tx_bit = randi([0 qam_size-1], num_qam_per_pkt, 1);
    
    % modulate bit stream
    tx_sym_data = qammod(tx_bit, qam_size, 'UnitAveragePower', true);
    
    % reshape data symbol stream
    tx_sym_data_dd = reshape(tx_sym_data, num_subc, []);
    
    % map data and pilot
    tx_sym_dd = tx_sym_data_dd;
end

% 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
tx_sym_tf = sqrt(num_sym/num_subc) * fft(ifft(tx_sym_dd, [], 2), [], 1);

% map to tf domain
idx_ofdm_map = [floor((num_ofdm_subc-num_subc)/2), floor((num_ofdm_sym-num_sym)/2)];
tx_sym_map_shift_tf = zeros(num_ofdm_subc, num_ofdm_sym);
tx_sym_map_shift_tf(idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(2)+1:idx_ofdm_map(2)+num_sym) = tx_sym_tf;
tx_sym_map_tf = fftshift(tx_sym_map_shift_tf, 1);

% ofdm modulate
tx_ofdm_sym = sqrt(num_ofdm_subc) * ifft(tx_sym_map_tf, [], 1);

% add cp (cyclic prefix)
tx_ofdm_sym_cp = [tx_ofdm_sym(end-len_cp+1:end, :); tx_ofdm_sym];       % cyclic-prefix

% serialize
tx_sig = tx_ofdm_sym_cp(:);

% calculate papr
papr_db = 10*log10(max(abs(tx_sig).^2)/mean(abs(tx_sig).^2));

% pass signal through channel
[tx_sig_faded, ch_path_gain] = fading_ch(tx_sig);       % ch_path_gain: normally constant per path, vary when doppler exists
ch_info = info(fading_ch);
ch_filter_coeff = ch_info.ChannelFilterCoefficients;    % ch_filter_coeff: channel filter coefficient per path
num_path = length(fading_ch.PathDelays);                % num_path: number of path

% reproduce real channel (time saving, option 1)
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

% compensate cfo (to observe impact of frequency shift)
cfo_sample = 0:length(rx_sig)-1;
cfo_vec = cfo_sample.'*cfo_norm/(num_ofdm_subc+len_cp);
rx_sig_cfo = rx_sig(:) .* exp(-1i*2*pi*cfo_vec);

% synchronization
rx_sig_synch = circshift(rx_sig_cfo, -len_cp-test_synch);

% reshape
rx_ofdm_sym_cp = reshape(rx_sig_synch, num_ofdm_subc+len_cp, []);

% remove cp (with synch: to observe impact of time shift)
rx_ofdm_sym = rx_ofdm_sym_cp(1:num_ofdm_subc, :);

% ofdm demodulate the symbol
rx_sym_map_tf = (1/sqrt(num_ofdm_subc)) * fft(rx_ofdm_sym, [], 1);

% demap
rx_sym_map_shift_tf = fftshift(rx_sym_map_tf, 1);
rx_sym_tf = rx_sym_map_shift_tf(idx_ofdm_map(1)+1:idx_ofdm_map(1)+num_subc, idx_ofdm_map(2)+1:idx_ofdm_map(2)+num_sym);

% 2d inverse sfft
rx_sym_dd = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf, [], 1), [], 2);

% extract diagonal elements of real channel
ch_est_tf_real = zeros(num_subc, num_sym);
for idx_sym = 1:num_sym
    ch_est_tf_real(:, idx_sym) = diag(ch_real_mat_tf(:, :, idx_sym));  % diagonal term only
end

% estimate channel
if strcmp(test_chest, 'dd')
    % demap pilot
    rx_sym_dd_circshift = circshift(rx_sym_dd, -floor((num_sym-num_pilot_sym)/2), 2);
    rx_sym_pilot_dd = rx_sym_dd_circshift(num_data_null_subc_head+num_tx_sym_data1(1)+1:num_data_null_subc_head+num_tx_sym_data1(1)+num_pilot_subc, 1:num_pilot_sym);
    
    % estimate dd channel
    scale_pilot = sqrt((num_subc*num_sym)/(num_pilot_subc*num_pilot_sym));
    ch_est_dd_circshift = zeros(num_subc, num_sym);
    ch_est_dd_circshift(floor(num_pilot_guard_subc/2)+1:num_pilot_subc-floor(num_pilot_guard_subc/2), floor(num_pilot_guard_sym/2)+1:num_pilot_sym-floor(num_pilot_guard_sym/2)) = ...
        scale_pilot*rx_sym_pilot_dd(floor(num_pilot_guard_subc/2)+1:num_pilot_subc-floor(num_pilot_guard_subc/2), floor(num_pilot_guard_sym/2)+1:num_pilot_sym-floor(num_pilot_guard_sym/2));
    ch_est_dd = circshift(ch_est_dd_circshift, [-floor(num_pilot_subc/2), -floor(num_pilot_sym/2)]);
    
    % estimate tf channel
    ch_est_tf = sqrt(num_sym/num_subc) * fft(ifft(ch_est_dd, [], 2), [], 1);
elseif strcmp(test_chest, 'perfect')
    % estimate tf channel
    ch_est_tf_perfect = rx_sym_tf ./ tx_sym_tf;     % tf domain
    ch_est_tf = ch_est_tf_perfect;
    
    % estimate dd channel
    ch_est_dd = sqrt(num_subc/num_sym) * fft(ifft(ch_est_tf, [], 1), [], 2);
elseif strcmp(test_chest, 'real')
    % get real channel
    ch_est_tf = ch_est_tf_real;
    
    % estimate dd channel
    ch_est_dd = sqrt(num_subc/num_sym) * fft(ifft(ch_est_tf, [], 1), [], 2);
else
    error('test_chest value must be one of these: {''dd'', ''perfect'', ''real''}')
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
    
    % observe symbols in dd domain
    rx_sym_dd_eq = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf_eq, [], 1), [], 2);
elseif strcmp(test_cheq, 'tf_mmse')
    % equalize channel in tf domain
    ch_est_tf_mmse = conj(ch_est_tf) ./ (noise_var+abs(ch_est_tf).^2);
    rx_sym_tf_eq = rx_sym_tf .* ch_est_tf_mmse;
    
    % observe symbols in dd domain
    rx_sym_dd_eq = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf_eq, [], 1), [], 2);
elseif strcmp(test_cheq, 'dd')
    % generate block circular channel matrix
    len_ch = [num_subc num_sym];
    len_ch_pad_head = len_ch/2;
    len_ch_pad_tail = len_ch/2;
    len_ch_pad_zero = len_ch;
    ch_est_dd_shift = fftshift(fftshift(ch_est_dd, 1), 2);
    base_ch = repmat(ch_est_dd_shift(end:-1:1, end:-1:1), 2, 2);
    sub_ch = zeros(len_ch(1)+len_ch_pad_zero(1), len_ch(1)+len_ch_pad_zero(1), len_ch(2)+len_ch_pad_zero(2));
    for j = 1 : len_ch(2)+len_ch_pad_zero(2)
        for k = 1 : len_ch(1)+len_ch_pad_zero(1)
            sub_ch(k, :, j) = circshift(base_ch(:, j), k-len_ch(1));
        end
    end
    sub_ch_crop = sub_ch(len_ch_pad_head(1)+1 : end-len_ch_pad_tail(1), 1:len_ch(1), :);
    sub_ch_reshape = reshape(sub_ch_crop, len_ch(1), len_ch(1)*(len_ch(2)+len_ch_pad_zero(2)));
    new_ch_full = zeros(len_ch(1)*(len_ch(2)+len_ch_pad_zero(2)), len_ch(1)*(len_ch(2)+len_ch_pad_zero(2)));
    for k = 1 : len_ch(2)+len_ch_pad_zero(2)
        new_ch_full(len_ch(1)*(k-1)+1 : len_ch(1)*k, :) = circshift(sub_ch_reshape, (k-len_ch(2))*len_ch(1), 2);
    end
    new_ch = new_ch_full(len_ch_pad_head(2)*len_ch(1)+1 : end-len_ch_pad_tail(2)*len_ch(1), 1:len_ch(2)*len_ch(1)) / sqrt(num_subc*num_sym);

    % equalize channel in dd domain
    rx_sym_dd_eq_vec = new_ch \ rx_sym_dd(:);
    rx_sym_dd_eq = reshape(rx_sym_dd_eq_vec, num_subc, num_sym);
else
    error('test_chest value must be one of these: {tf, tf_mmse, dd}')
end

% demap data symbols
if strcmp(test_chest, 'dd')
    rx_sym_dd_circshift = circshift(rx_sym_dd_eq, -floor((num_sym-num_pilot_sym)/2), 2);
    rx_sym_data1_dd = rx_sym_dd_circshift(num_data_null_subc_head+1:num_data_null_subc_head+num_tx_sym_data1(1), :);
    rx_sym_data2_dd = rx_sym_dd_circshift(num_data_null_subc_head+num_tx_sym_data1(1)+1:num_data_null_subc_head+num_tx_sym_data1(1)+num_tx_sym_data2(1), num_pilot_sym+1:end);
    rx_sym_data3_dd = rx_sym_dd_circshift(num_data_null_subc_head+num_tx_sym_data1(1)+num_tx_sym_data2(1)+1:num_data_null_subc_head+num_tx_sym_data1(1)+num_tx_sym_data2(1)+num_tx_sym_data3(1), :);
    rx_sym_data = [rx_sym_data1_dd(:); rx_sym_data2_dd(:); rx_sym_data3_dd(:)];
else
    rx_sym_data = rx_sym_dd_eq(:);
end

% demodulate qam symbols
rx_bit = qamdemod(rx_sym_data, qam_size, 'UnitAveragePower', true);

% (1) calculate bit errors
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
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_dd)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd_eq))
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_dd)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd))
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_tf)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_tf))
    figure, subplot(1, 2, 1), mesh(1:num_ofdm_sym, 1:num_ofdm_subc, abs(tx_sym_map_tf)), subplot(1, 2, 2), mesh(1:num_ofdm_sym, 1:num_ofdm_subc, abs(rx_sym_map_tf))
    figure, subplot(1, 2, 1), mesh(1:num_ofdm_sym, 1:num_ofdm_subc, abs(tx_ofdm_sym)), subplot(1, 2, 2), mesh(1:num_ofdm_sym, 1:num_ofdm_subc, abs(rx_ofdm_sym))
    figure, subplot(1, 2, 1), mesh(1:num_ofdm_sym, 1:num_ofdm_subc+len_cp, abs(tx_ofdm_sym_cp)), subplot(1, 2, 2), mesh(1:num_ofdm_sym, 1:num_ofdm_subc+len_cp, abs(rx_ofdm_sym_cp))
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(ch_est_tf)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(ch_est_dd))
    
    % channel test
    test_rx_sym_dd_est_tail = (1/sqrt(num_subc*num_sym))*conv2(repmat(ch_est_dd,2,2), tx_sym_dd, 'full');
    test_rx_sym_dd_est = test_rx_sym_dd_est_tail(num_subc+1:2*num_subc, num_sym+1:2*num_sym);
    figure
    subplot(2, 1, 1), plot(real(rx_sym_dd(:)), '-b.'), hold on, plot(real(test_rx_sym_dd_est(:)), ':r.'), hold off, title('real(rx sig) and real(rx sig est) in dd-domain'), grid minor
    subplot(2, 1, 2), plot(imag(rx_sym_dd(:)), '-b.'), hold on, plot(imag(test_rx_sym_dd_est(:)), ':r.'), hold off, title('imag(rx sig) and imag(rx sig est) in dd-domain'), grid minor
    
%     ch_est_tf_perfect = rx_sym_tf ./ tx_sym_tf;
%     ch_est_dd_perfect = sqrt(num_subc/num_sym) * fft(ifft(ch_est_tf_perfect, [], 1), [], 2);
%     figure
%     subplot(2, 1, 1), plot(real(tx_sym_dd(:)), '-b.'), hold on, plot(real(rx_sym_dd_eq(:)), ':r.'), hold off, grid minor
%     subplot(2, 1, 2), plot(imag(tx_sym_dd(:)), '-b.'), hold on, plot(imag(rx_sym_dd_eq(:)), ':r.'), hold off, grid minor
%     figure
%     subplot(2, 1, 1), plot(real(tx_sym_dd(:)), '-b.'), hold on, plot(real(rx_sym_dd(:)), ':r.'), hold off, grid minor
%     subplot(2, 1, 2), plot(imag(tx_sym_dd(:)), '-b.'), hold on, plot(imag(rx_sym_dd(:)), ':r.'), hold off, grid minor
%     figure
%     subplot(2, 1, 1), plot(real(ch_est_dd_perfect(:)), '-b.'), hold on, plot(real(ch_est_dd(:)), ':r.'), hold off, grid minor
%     subplot(2, 1, 2), plot(imag(ch_est_dd_perfect(:)), '-b.'), hold on, plot(imag(ch_est_dd(:)), ':r.'), hold off, grid minor
end

% dump
if test_scope
%     assignin('base', 'tx_sym_data_dd', tx_sym_data_dd);
    assignin('base', 'tx_sym_dd', tx_sym_dd);
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
    assignin('base', 'rx_sym_dd', rx_sym_dd);
    assignin('base', 'rx_sym_dd_eq', rx_sym_dd_eq);
    assignin('base', 'rx_sym_data_dd', rx_sym_data);
    
    assignin('base', 'ch_est_tf', ch_est_tf);
    assignin('base', 'ch_est_dd', ch_est_dd);
    assignin('base', 'test_ch', test_ch);
    
    if strcmp(test_cheq, 'dd')
        assignin('base', 'new_ch', new_ch);
    end
end

end
