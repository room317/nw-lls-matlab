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

function [qam_error_tfeq, qam_error_ddeq, num_qam_per_pkt] = test_dd_pilot(snr_db, test_pilot, test_cp_position, test_synch, test_scope, test_seed)

% set parameter
num_fft = 64; % 64;
num_subc = 32; % 32;
num_sym = 32; % 32;
num_pilot_subc = 16; % delay spread axis
num_guard_subc = 8;  % delay spread axis
f_s = 15e3*num_subc;
qam_size = 16;
len_cp = num_fft/4; % 0;
% test_synch = len_cp; % 7; % len_cp in case of prefix
cfo = 0; % 0.00045; % 0.0002;
carrier_freq_mhz = 4000;
velocity_kmh = 120; % 120;
idx_fading = 9; % 9: TDL-A, 12: TDL-D
delay_spread_rms_us = 0.01; %0.1e-6;
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

% calculate qam symbols per packet
if test_pilot
    num_qam_per_pkt = (num_subc - num_pilot_subc) * num_sym;
else
    num_qam_per_pkt = num_subc * num_sym;
end

% calculate noise variance
noise_var = 10 ^ ((-0.1)*snr_db);

if test_pilot
    % generate bit stream
    tx_bit = randi([0 qam_size-1], (num_subc-num_pilot_subc)*num_sym, 1);
    
    % modulate bit stream
    tx_sym_data = qammod(tx_bit, qam_size, 'UnitAveragePower', true);
    
    % reshape
    tx_sym_data_dd = reshape(tx_sym_data, num_subc-num_pilot_subc, []);
    
    % pilot generation
    tx_sym_pilot_dd = zeros(num_pilot_subc, num_sym);
    tx_sym_pilot_dd((num_pilot_subc/2)+1, (num_sym/2)+1) = sqrt(num_pilot_subc*num_sym);
    
    % map data and pilot
    tx_sym_dd = [...
        tx_sym_data_dd(1:(num_subc-num_pilot_subc)/2, :);
        tx_sym_pilot_dd;
        tx_sym_data_dd((num_subc-num_pilot_subc)/2+1:end, :)];
else
    % generate bit stream
    tx_bit = randi([0 qam_size-1], num_subc*num_sym, 1);
    
    % modulate bit stream
    tx_sym = qammod(tx_bit, qam_size, 'UnitAveragePower', true);
    
    % reshape
    tx_sym_dd = reshape(tx_sym, num_subc, []);
end

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
rx_sig = awgn(tx_sig_faded, snr_db, 'measured');

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

% (1) estimate tf channel
ch_est_tf_perfect = rx_sym_tf ./ tx_sym_tf;     % tf domain

% (1) equalize channel in tf domain
ch_est_tf_mmse_perfect = conj(ch_est_tf_perfect) ./ (noise_var+abs(ch_est_tf_perfect).^2);
rx_sym_tf_tfeq = rx_sym_tf .* ch_est_tf_mmse_perfect;
% rx_sym_tf_tfeq = rx_sym_tf ./ ch_est_tf_perfect;

% (1) observe symbols in dd domain
rx_sym_dd_tfeq = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf_tfeq, [], 1), [], 2);

% (1) demap data symbols
if test_pilot
    rx_sym_data_dd_tfeq = [...
        rx_sym_dd_tfeq(1:(num_subc-num_pilot_subc)/2, :);
        rx_sym_dd_tfeq((num_subc+num_pilot_subc)/2+1:end, :)];
else
    rx_sym_data_dd_tfeq = rx_sym_dd_tfeq;
end

% (1) demodulate qam symbols
rx_bit_tfeq = qamdemod(rx_sym_data_dd_tfeq(:), qam_size, 'UnitAveragePower', true);

% (1) calculate bit errors
if isempty(tx_bit)
    qam_error_tfeq = 0;
else
    qam_error_tfeq = symerr(tx_bit, rx_bit_tfeq);
end

% plot
if test_scope
    figure, scatter(real(tx_sym_dd(:)), imag(tx_sym_dd(:)), 'bx'), hold on, scatter(real(rx_sym_dd_tfeq(:)), imag(rx_sym_dd_tfeq(:)), 'r.'), grid, axis([-2 2 -2 2]), hold off, title('constellation tx symbol and rx symbol after tf-eq(mmse)')
    figure, plot(real(tx_sym_dd(:)), '-b.'), hold on, plot(real(rx_sym_dd_tfeq(:)), ':r.'), hold off, title('real parts of tx and rx symbols after tf-eq(mmse)')
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_dd)), title('abs of tx symbols in dd-domain'), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd_tfeq)), title('abs of rx symbols in dd-domain after tf-eq(mmse)')
    pause
end

% (2) estimate dd channel
if test_pilot
    ch_est_dd_perfect = [...
        zeros((num_subc-(num_pilot_subc-num_guard_subc))/2, num_sym);
        rx_sym_dd((num_subc-(num_pilot_subc-num_guard_subc))/2+1:(num_subc+(num_pilot_subc-num_guard_subc))/2, :);
        zeros((num_subc-(num_pilot_subc-num_guard_subc))/2, num_sym)];
else
    ch_est_dd_perfect_shift = sqrt(num_subc/num_sym) * fft(ifft(ch_est_tf_perfect, [], 1), [], 2);    % dd domain
    ch_est_dd_perfect = fftshift(fftshift(ch_est_dd_perfect_shift, 1), 2);
end

% (2) evaluate dd channel (compare real rx signal and estimated rx signal)
if test_scope
    rx_sym_dd_est_shift = filter2(repmat(ch_est_dd_perfect(end:-1:1, end:-1:1), 2, 2), tx_sym_dd)/sqrt(num_subc*num_sym);
    rx_sym_dd_est = fftshift(fftshift(rx_sym_dd_est_shift, 1), 2);
    
    % plot
    figure, plot(abs(rx_sym_dd(:)), '-b.'), hold on, plot(abs(rx_sym_dd_est(:)), ':r.'), hold off, title('abs of real rx signal and estimated rx signal')
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd)), title('abs of real rx signal in dd-domain'), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd_est)), title('abs of estimated rx signal in dd-domain')
    pause
end

% (3) generate block toeplitz(circular) dd channel matrix
new_tx = tx_sym_dd(:);
new_rx = rx_sym_dd(:);

len_ch = [num_subc num_sym];
len_ch_pad_head = len_ch/2;
len_ch_pad_tail = len_ch/2;
len_ch_pad_zero = len_ch;

base_ch = repmat(ch_est_dd_perfect(end:-1:1, end:-1:1), 2, 2);

% (3) generate block toeplitz(circular) dd channel matrix
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

% (3) evaluate block toeplitz(circular) dd channel matrix
if test_scope
    new_rx_est = new_ch * new_tx;

    % plot
    fprintf('Block Circular Matrix    Rank: %3d    Cond. Num.: %7.2f\n', rank(new_ch), cond(new_ch))
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(ch_est_dd_perfect)), title('abs of channel in dd-domain'), subplot(1, 2, 2), mesh(1:num_sym*num_subc, 1:num_subc*num_sym, abs(new_ch)), title('abs of block circular matrix')
    figure, plot(abs(new_rx(:)), '-b.'), hold on, plot(abs(new_rx_est(:)), ':r.'), hold off, title('abs of real rx signal and rx signal estimated with block circular matrix')
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(reshape(new_rx, num_subc, num_sym))), title('abs of real rx signal in dd-domain'), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(reshape(new_rx_est, num_subc, num_sym))), title('abs of rx signal estimated with block circular matrix')
    figure, plot(abs(new_rx), '-b.'), hold on, plot(abs(new_rx-new_rx_est), ':r.'), hold off, title('abs of real rx signal and estimation error')
    pause
end

% (3) equalize channel in dd domain
rx_sym_dd_ddeq_vec = new_ch \ rx_sym_dd(:);
% rx_sym_dd_ddeq_vec = pinv(new_ch) * rx_sym_dd(:);
rx_sym_dd_ddeq = reshape(rx_sym_dd_ddeq_vec, num_subc, num_sym);

% (3) demap data symbols
if test_pilot
    rx_sym_data_dd_ddeq = [...
        rx_sym_dd_ddeq(1:(num_subc-num_pilot_subc)/2, :);
        rx_sym_dd_ddeq((num_subc+num_pilot_subc)/2+1:end, :)];
else
    rx_sym_data_dd_ddeq = rx_sym_dd_ddeq;
end

% (3) observe symbols in dd domain
rx_bit_ddeq = qamdemod(rx_sym_data_dd_ddeq(:), qam_size, 'UnitAveragePower', true);
if isempty(tx_bit)
    qam_error_ddeq = 0;
else
    qam_error_ddeq = symerr(tx_bit, rx_bit_ddeq);
end

% plot
if test_scope
    figure, scatter(real(tx_sym_dd(:)), imag(tx_sym_dd(:)), 'bx'), hold on, scatter(real(rx_sym_dd_ddeq(:)), imag(rx_sym_dd_ddeq(:)), 'r.'), grid, axis([-2 2 -2 2]), hold off, title('constellation tx symbol and rx symbol after dd-eq')
    figure, plot(real(tx_sym_dd(:)), '-b.'), hold on, plot(real(rx_sym_dd_ddeq(:)), ':r.'), hold off, title('real parts of tx and rx symbols after dd-eq')
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_dd)), title('abs of tx symbols in dd-domain'), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd_ddeq)), title('abs of rx symbols in dd-domain after dd-eq')
end
% figure, plot(real(rx_sym_dd_tfeq(:)), '-b.'), hold on, plot(real(rx_sym_dd_ddeq(:)), ':r.'), hold off, title('real parts of rx symbols after eq'), legend('after tf-eq', 'after dd-eq')
% pause

% dump
if test_scope
    assignin('base', 'ch_est_tf_perfect', ch_est_tf_perfect);
    assignin('base', 'ch_est_dd_perfect', ch_est_dd_perfect);
    assignin('base', 'new_ch', new_ch);
    
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
end

end
