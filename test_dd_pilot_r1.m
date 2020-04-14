% test dd pilot
%   - snr_db: snr in db scale
%   - test_pilot: transmit pilot first to get real channel value
%   - test_cp_position = 'prefix' or 'postfix'
%   - test_synch: synch position
%   - test_scope: print and plot results
% created: 2020.03.06
% modified:
%   - 2020.03.02: prefix/postfix option
%   - 2020.03.03: rician channel added
%   - 2020.03.06: pilot estimation added

function [qam_error, num_qam_per_pkt] = test_dd_pilot_r1(snr_db, test_pilot, test_cp_position, test_synch, test_scope, test_seed, test_chest, test_cheq)

% set parameter
num_fft = 64; % 64;
num_subc = 32; % 32;
num_sym = 32; % 32;
num_pilot_subc = 8; % delay spread axis
num_guard_subc = num_pilot_subc; % 8;  % delay spread axis
num_data_subc = num_subc-num_pilot_subc-num_guard_subc;
f_s = 15e3*num_fft;
qam_size = 16;
len_cp = 8; % num_guard_subc/2; % num_fft/4; % 0;
% test_synch = len_cp; % 7; % len_cp in case of prefix
cfo = 0; % 0.00045; % 0.0002;
carrier_freq_mhz = 4000;
velocity_kmh = 120; % 120;
idx_fading = 9; % 9: TDL-A, 12: TDL-D
delay_spread_rms_us = 0.1; %0.1e-6;
test_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, snr_db, delay_spread_rms_us);

% check error
if test_synch > len_cp
    error('test_synch should not be greater than %d.\n', len_cp)
end

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
            'NormalizePathGains', true);
    otherwise
        fading_ch = comm.RayleighChannel(...
            'SampleRate', f_s, ...
            'PathDelays', test_ch.path_delays, ...
            'AveragePathGains', test_ch.average_path_gains, ...
            'NormalizePathGains', true, ...
            'MaximumDopplerShift', test_ch.maximum_doppler_shift, ...
            'DopplerSpectrum', test_ch.doppler_spectrum);
end

% calculate noise variance
noise_var = 10 ^ ((-0.1)*snr_db);

if test_pilot
    % calculate num. qam symbols per packet
    num_qam_per_pkt = num_data_subc*num_sym;
    
    % generate bit stream
    tx_bit = randi([0 qam_size-1], num_qam_per_pkt, 1);
    
    % modulate bit stream
    tx_sym_data = qammod(tx_bit, qam_size, 'UnitAveragePower', true);
    
    % reshape data symbol stream
    tx_sym_data_dd = reshape(tx_sym_data, num_data_subc, []);
    
    % generate pilot symbols
    tx_sym_pilot_dd = zeros(num_pilot_subc, num_sym);
    tx_sym_pilot_dd(1, 1) = sqrt(num_pilot_subc*num_sym);
    
    % generate guard symbols
    tx_sym_guard_dd = zeros(num_guard_subc/2, num_sym);
else
    % calculate num. qam symbols per packet
    num_qam_per_pkt = num_subc*num_sym;
    
    % generate bit stream
    tx_bit = randi([0 qam_size-1], num_qam_per_pkt, 1);
    
    % modulate bit stream
    tx_data_sym = qammod(tx_bit, qam_size, 'UnitAveragePower', true);
    
    % reshape data symbol stream
    tx_sym_data_dd = reshape(tx_data_sym, num_subc, []);
    
    % generate pilot symbols
    tx_sym_pilot_dd = [];
    
    % generate guard symbols
    tx_sym_guard_dd = [];
end
    
% map data and pilot
tx_sym_dd = [...
    tx_sym_guard_dd;
    tx_sym_pilot_dd;
    tx_sym_data_dd;
    tx_sym_guard_dd];

% 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
tx_sym_tf = sqrt(num_sym/num_subc) * fft(ifft(tx_sym_dd, [], 2), [], 1);

% map to tf domain
tx_sym_map_shift_tf = zeros(num_fft, num_sym);
tx_sym_map_shift_tf(num_fft/2-num_subc/2+1:num_fft/2+num_subc/2, :) = tx_sym_tf;
tx_sym_map_tf = fftshift(tx_sym_map_shift_tf, 1);

% ofdm modulate
tx_ofdm_sym = sqrt(num_fft) * ifft(tx_sym_map_tf, [], 1);

% add cp (cyclic postfix)
if strcmp(test_cp_position, 'prefix')
    tx_ofdm_sym_cp = [tx_ofdm_sym(end-len_cp+1:end, :); tx_ofdm_sym];       % cyclic-prefix
else
    tx_ofdm_sym_cp = [tx_ofdm_sym; tx_ofdm_sym(1:len_cp, :)];       % cyclic-postfix
end

% serialize
tx_sig = tx_ofdm_sym_cp(:);

% pass signal through channel
tx_sig_faded = fading_ch(tx_sig);

% add gaussian noise
if test_scope
    rx_sig = tx_sig_faded;
else
    rx_sig = awgn(tx_sig_faded, snr_db, 'measured');
end

% compensate cfo (to observe impact of frequency shift)
% cfo_time = 0:length(rx_sig)-1;
% rx_sig_cfo = rx_sig(:) .* exp(-1i*2*pi*cfo*cfo_time(:));
cfo_time = 0:length(rx_sig)-1;
cfo_delta = 0*cfo;
cfo = cfo+cfo_delta*cfo_time;
rx_sig_cfo = rx_sig(:) .* exp(-1i*2*pi*cfo.'.*cfo_time(:));

% reshape
rx_ofdm_sym_cp = reshape(rx_sig_cfo, num_fft+len_cp, num_sym);

% remove cp (with synch: to observe impact of time shift)
rx_ofdm_sym = rx_ofdm_sym_cp(test_synch+1:test_synch+num_fft, :);

% ofdm demodulate the symbol
rx_sym_map_tf = (1/sqrt(num_fft)) * fft(rx_ofdm_sym, [], 1);

% demap
rx_sym_map_shift_tf = fftshift(rx_sym_map_tf, 1);
rx_sym_tf = rx_sym_map_shift_tf(num_fft/2-num_subc/2+1:num_fft/2+num_subc/2, :);

% 2d inverse sfft
rx_sym_dd = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf, [], 1), [], 2);

% estimate channel
if strcmp(test_chest, 't')
    % crop t signals
    if strcmp(test_cp_position, 'prefix')
        idx_t_pilot_begin = 1;
        idx_t_pilot_end = len_cp+((num_guard_subc/2)+num_pilot_subc)*(num_fft/num_subc);
    else
        idx_t_pilot_begin = 1;
        idx_t_pilot_end = ((num_guard_subc/2)+num_pilot_subc)*(num_fft/num_subc);
    end
    tx_pilot_nfft_t = zeros(num_fft, num_sym);
    tx_pilot_nfft_t(idx_t_pilot_begin:idx_t_pilot_end, :) = tx_ofdm_sym_cp(idx_t_pilot_begin:idx_t_pilot_end,:);
    rx_pilot_nfft_t = zeros(num_fft, num_sym);
    rx_pilot_nfft_t(idx_t_pilot_begin:idx_t_pilot_end, :) = rx_ofdm_sym_cp(idx_t_pilot_begin:idx_t_pilot_end,:);
    
    % transform to tf domain
    tx_pilot_nfft_tf = (1/sqrt(num_fft)) * fft(tx_pilot_nfft_t, [], 1) * sqrt(num_subc/num_pilot_subc);
    rx_pilot_nfft_tf = (1/sqrt(num_fft)) * fft(rx_pilot_nfft_t, [], 1) * sqrt(num_subc/num_pilot_subc);
    tx_pilot_nfft_shift_tf = fftshift(tx_pilot_nfft_tf, 1);
    tx_pilot_tf = tx_pilot_nfft_shift_tf(num_fft/2-num_subc/2+1:num_fft/2+num_subc/2, :);
    rx_pilot_nfft_shift_tf = fftshift(rx_pilot_nfft_tf, 1);
    rx_pilot_tf = rx_pilot_nfft_shift_tf(num_fft/2-num_subc/2+1:num_fft/2+num_subc/2, :);
    
    % estimate tf channel
    ch_est_tf = rx_pilot_tf ./ tx_pilot_tf;
    
    % estimate dd channel
    ch_est_dd = sqrt(num_subc/num_sym) * fft(ifft(ch_est_tf, [], 1), [], 2);
elseif strcmp(test_chest, 'tf')
    % estimate tf channel
    ch_est_tf_perfect = rx_sym_tf ./ tx_sym_tf;     % tf domain
    ch_est_tf = ch_est_tf_perfect;
    
    % estimate dd channel
    ch_est_dd = sqrt(num_subc/num_sym) * fft(ifft(ch_est_tf, [], 1), [], 2);
elseif strcmp(test_chest, 'dd')
    % estimate dd channel
    idx_dd_pilot_begin = (num_guard_subc/2)+1;
    idx_dd_pilot_end = (num_guard_subc/2)+num_pilot_subc;
    ch_est_dd = zeros(num_subc, num_sym);
    ch_est_dd(1:num_pilot_subc, :) = sqrt(num_subc/num_pilot_subc)*rx_sym_dd(idx_dd_pilot_begin:idx_dd_pilot_end, :);
    
    % estimate tf channel
    ch_est_tf = sqrt(num_sym/num_subc) * fft(ifft(ch_est_dd, [], 2), [], 1);
else
    error('test_chest value must be one of these: {t, tf, dd}')
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
if test_pilot
    rx_sym_data_dd = rx_sym_dd_eq((num_guard_subc/2)+num_pilot_subc+1:(num_guard_subc/2)+num_pilot_subc+num_data_subc, :);
else
    rx_sym_data_dd = rx_sym_dd_eq;
end

% demodulate qam symbols
rx_bit = qamdemod(rx_sym_data_dd(:), qam_size, 'UnitAveragePower', true);

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
    subplot(2, 1, 1), plot(real(tx_sym_data_dd(:)), '-b.'), hold on, plot(real(rx_sym_data_dd(:)), ':r.'), hold off, grid minor
    subplot(2, 1, 2), plot(imag(tx_sym_data_dd(:)), '-b.'), hold on, plot(imag(rx_sym_data_dd(:)), ':r.'), hold off, grid minor
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_dd)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd_eq)), grid minor
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_dd)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd))
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_tf)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_tf))
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_fft, abs(tx_sym_map_tf)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_fft, abs(rx_sym_map_tf))
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_fft, abs(tx_ofdm_sym)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_fft, abs(rx_ofdm_sym))
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_fft+len_cp, abs(tx_ofdm_sym_cp)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_fft+len_cp, abs(rx_ofdm_sym_cp))
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(ch_est_tf)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(ch_est_dd))
    
    if strcmp(test_chest, 't')
        figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_fft, abs(tx_pilot_nfft_t)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_fft, abs(rx_pilot_nfft_t))
        figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_fft, abs(tx_pilot_nfft_tf)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_fft, abs(rx_pilot_nfft_tf))
        figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_pilot_tf)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_pilot_tf))
    elseif strcmp(test_chest, 'dd')
        test_rx_sym_dd_est_tail = (1/sqrt(num_subc*num_sym))*conv2(repmat(ch_est_dd,2,2), tx_sym_dd, 'full');
        test_rx_sym_dd_est = test_rx_sym_dd_est_tail(num_subc+1:2*num_subc, num_sym+1:2*num_sym);
        figure
        subplot(2, 1, 1), plot(real(rx_sym_dd(:)), '-b.'), hold on, plot(real(test_rx_sym_dd_est(:)), ':r.'), hold off, title('real(rx sig) and real(rx sig est) in dd-domain')
        subplot(2, 1, 2), plot(imag(rx_sym_dd(:)), '-b.'), hold on, plot(imag(test_rx_sym_dd_est(:)), ':r.'), hold off, title('imag(rx sig) and imag(rx sig est) in dd-domain')
    end
end

% dump
if test_scope
    assignin('base', 'tx_sym_data_dd', tx_sym_data_dd);
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
    assignin('base', 'rx_sym_data_dd', rx_sym_data_dd);
    
    assignin('base', 'ch_est_tf', ch_est_tf);
    assignin('base', 'ch_est_dd', ch_est_dd);
    
    assignin('base', 'test_ch', test_ch);
    
    if strcmp(test_cheq, 'dd')
        assignin('base', 'new_ch', new_ch);
    end
end

end
