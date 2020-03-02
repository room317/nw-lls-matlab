% test channel convolution
% created: 2020.02.11
% modified:

function [pkt_error_tfeq, pkt_error_ddeq] = test_ch_conv_r1(snr_db, test_pilot, test_signal, test_gauss_ch)

% set parameter
num_fft = 64; % 64;
num_subc = 32; % 32;
num_sym = 32; % 32;
f_s = 15e3*num_subc;
qam_size = 16;
len_cp = num_fft/4; % 0;
synch = 7; % 7;
cfo = 0; % 0.00045; % 0.0002;
carrier_freq_mhz = 4000;
velocity_kmh = 120; % 120;
idx_fading = 9; % 9: TDL-A, 12: TDL-D
delay_spread_rms = 1; %0.1e-6;
test_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, snr_db, delay_spread_rms);

if test_pilot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% estimate real channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % create a rayleigh fading channel object
%     random_seed = rng;
    random_seed = 74;
    fading_ch = comm.RayleighChannel(...
        'SampleRate', f_s, ...
        'PathDelays', test_ch.path_delays, ...
        'AveragePathGains', test_ch.average_path_gains, ...
        'NormalizePathGains', true, ...
        'MaximumDopplerShift', test_ch.maximum_doppler_shift, ...
        'DopplerSpectrum', test_ch.doppler_spectrum);
    %     'RandomStream', 'mt19937ar with seed', ...
    %     'Seed', ch_param.seed);
    %     'PathGainsOutputPort', true);
    %     'Visualization', 'Impulse response');   % enables the impulse response channel visualization

    % generate pilot tone (single tone)
    tx_sym_dd = zeros(num_subc, num_sym);
    tx_sym_dd((num_subc/2)+1, (num_sym/2)+1) = sqrt(num_subc*num_sym);

    % 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
    tx_sym_tf = sqrt(num_sym/num_subc) * fft(ifft(tx_sym_dd, [], 2), [], 1);

    % map to tf domain
    tx_sym_map_shift_tf = zeros(num_fft, num_sym);
    tx_sym_map_shift_tf(num_fft/2-num_subc/2+1:num_fft/2+num_subc/2, :) = tx_sym_tf;
    tx_sym_map_tf = fftshift(tx_sym_map_shift_tf, 1);

    % ofdm modulate
    tx_ofdm_sym = sqrt(num_subc) * ifft(tx_sym_map_tf, [], 1);

    % add cp (cyclic postfix)
    % tx_ofdm_sym_cp = [tx_ofdm_sym(end-len_cp+1:end, :); tx_ofdm_sym];       % cyclic-prefix
    tx_ofdm_sym_cp = [tx_ofdm_sym; tx_ofdm_sym(1:len_cp, :)];       % cyclic-postfix

    % serialize
    tx_sig = tx_ofdm_sym_cp(:);

    % pass signal through channel
    rng(random_seed)
    tx_sig_faded = fading_ch(tx_sig);

    % add gaussian noise
    rx_sig = tx_sig_faded;

    % compensate cfo
    cfo_time = 0:length(rx_sig)-1;
    rx_sig_cfo = rx_sig(:) .* exp(-1i*2*pi*cfo*cfo_time(:));

    % reshape
    rx_ofdm_sym_cp = reshape(rx_sig_cfo, num_fft+len_cp, num_sym);

    % remove cp (with synch)
    rx_ofdm_sym = rx_ofdm_sym_cp(synch+1:synch+num_fft, :);

    % ofdm demodulate the symbol
    rx_sym_map_tf = (1/sqrt(num_subc)) * fft(rx_ofdm_sym, [], 1);

    % demap
    rx_sym_map_shift_tf = fftshift(rx_sym_map_tf, 1);
    rx_sym_tf = rx_sym_map_shift_tf(num_fft/2-num_subc/2+1:num_fft/2+num_subc/2, :);

    % estimate tf channel
    ch_est_tf_perfect = rx_sym_tf ./ tx_sym_tf;     % tf domain
    ch_tf_real = ch_est_tf_perfect;

    % estimate dd channel
    ch_est_dd_perfect_shift = sqrt(num_subc/num_sym) * fft(ifft(ch_tf_real, [], 1), [], 2);    % dd domain
    ch_dd_real_shift = ch_est_dd_perfect_shift;
    ch_dd_real = fftshift(fftshift(ch_dd_real_shift, 1), 2);

    % % test (observe trx signals in various domains)
    % tx_dd = tx_sym_dd;
    % tx_ft = tx_sym_tf;
    % tx_fd = fft(tx_ft, [], 2);
    % tx_dt = tx_ofdm_sym;
    % rx_dd = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf, [], 1), [], 2);
    % rx_ft = rx_sym_tf;
    % rx_fd = fft(rx_ft, [], 2);
    % rx_dt = rx_ofdm_sym;
    % figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_dd)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_dd))
    % figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_ft)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_ft))
    % figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_fd)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_fd))
    % figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_dt)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_dt))
    % pause
else
    ch_tf_real = [];
    ch_dd_real = [];
    ch_dd_real_shift = [];
end

%%%%%%%%%%%%%
%%%%%%%%%% go
%%%%%%%%%%%%%

% create a rayleigh fading channel object
if test_pilot
    rng(random_seed)
end
fading_ch = comm.RayleighChannel(...
    'SampleRate', f_s, ...
    'PathDelays', test_ch.path_delays, ...
    'AveragePathGains', test_ch.average_path_gains, ...
    'NormalizePathGains', true, ...
    'MaximumDopplerShift', test_ch.maximum_doppler_shift, ...
    'DopplerSpectrum', test_ch.doppler_spectrum);
%     'RandomStream', 'mt19937ar with seed', ...
%     'Seed', ch_param.seed);
%     'PathGainsOutputPort', true);
%     'Visualization', 'Impulse response');   % enables the impulse response channel visualization

% calculate noise variance
noise_var = 10 ^ ((-0.1)*snr_db);

% generate signal
switch test_signal
    case 1      % (1) single-tone
        tx_sym_dd = zeros(num_subc, num_sym);
        tx_sym_dd((num_subc/2)+1, (num_sym/2)+1) = sqrt(num_subc*num_sym);
    case 2      % (4) two-tones
        tx_sym_dd = zeros(num_subc, num_sym);
        tx_sym_dd((num_subc/2)+1, (num_sym/2)+1) = sqrt(num_subc*num_sym/2);
        tx_sym_dd((num_subc*3/8)+1, (num_sym/2)+10) = 30;
    case 3      % (5) four-tones
        tx_sym_dd = zeros(num_subc, num_sym);
        tx_sym_dd((num_subc/4)+1, (num_sym/4)+1) = sqrt(num_subc*num_sym/4);
        tx_sym_dd((num_subc/4)+1, (num_sym*3/4)+1) = 20;
        tx_sym_dd((num_subc*3/4)+1, (num_sym/4)+1) = 10;
        tx_sym_dd((num_subc*3/4)+1, (num_sym*3/4)+1) = 5;
    case 4      % (3) guard-banded symbols
        tx_sym_dd = zeros(num_subc, num_sym);
        sym_gap = [10 10];
        tx_sym_dd(sym_gap(1)/2+1:end-sym_gap(1)/2, sym_gap(2)/2+1:end-sym_gap(2)/2) = complex(randn(num_subc-sym_gap(1), num_sym-sym_gap(2)));
    otherwise
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
tx_ofdm_sym = sqrt(num_subc) * ifft(tx_sym_map_tf, [], 1);

% add cp (cyclic postfix)
% tx_ofdm_sym_cp = [tx_ofdm_sym(end-len_cp+1:end, :); tx_ofdm_sym];       % cyclic-prefix
tx_ofdm_sym_cp = [tx_ofdm_sym; tx_ofdm_sym(1:len_cp, :)];       % cyclic-postfix

% serialize
tx_sig = tx_ofdm_sym_cp(:);

% pass signal through channel
if test_pilot
    rng(random_seed)
end
tx_sig_faded = fading_ch(tx_sig);

% add gaussian noise
if ismember(test_signal, [1 2 3 4])
    rx_sig = tx_sig_faded;
else
    rx_sig = awgn(tx_sig_faded, snr_db, 'measured');
end

% compensate cfo
cfo_time = 0:length(rx_sig)-1;
rx_sig_cfo = rx_sig(:) .* exp(-1i*2*pi*cfo*cfo_time(:));

% reshape
rx_ofdm_sym_cp = reshape(rx_sig_cfo, num_fft+len_cp, num_sym);

% remove cp (with synch)
rx_ofdm_sym = rx_ofdm_sym_cp(synch+1:synch+num_fft, :);

% ofdm demodulate the symbol
rx_sym_map_tf = (1/sqrt(num_subc)) * fft(rx_ofdm_sym, [], 1);

% demap
rx_sym_map_shift_tf = fftshift(rx_sym_map_tf, 1);
rx_sym_tf = rx_sym_map_shift_tf(num_fft/2-num_subc/2+1:num_fft/2+num_subc/2, :);

% 2d inverse sfft
rx_sym_dd = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf, [], 1), [], 2);

% (1) estimate tf channel
if test_pilot
    ch_est_tf_perfect = ch_tf_real;     % tf domain
else
    ch_est_tf_perfect = rx_sym_tf ./ tx_sym_tf;     % tf domain
end

if ismember(test_signal, [1 2 3 4])
    pkt_error_tfeq = 0;
else
    % (1) equalize channel in tf domain
    ch_est_tf_mmse_perfect = conj(ch_est_tf_perfect) ./ (noise_var+abs(ch_est_tf_perfect).^2);
    rx_sym_tf_tfeq = rx_sym_tf .* ch_est_tf_mmse_perfect;
    
    % (1) observe symbols in dd domain
    rx_sym_dd_tfeq = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf_tfeq, [], 1), [], 2);
    rx_bit_tfeq = qamdemod(rx_sym_dd_tfeq(:), qam_size, 'UnitAveragePower', true);
    pkt_error_tfeq = symerr(tx_bit, rx_bit_tfeq);
    
    % plot
    figure, scatter(real(tx_sym_dd(:)), imag(tx_sym_dd(:)), 'bo'), hold on, scatter(real(rx_sym_dd_tfeq(:)), imag(rx_sym_dd_tfeq(:)), 'r.'), grid, axis([-2 2 -2 2]), hold off
    figure, plot(real(tx_sym_dd(:)), '-b.'), hold on, plot(real(rx_sym_dd_tfeq(:)), ':r.'), hold off
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_dd)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd_tfeq))
    pause
end

% (2) estimate dd channel
if test_pilot
    ch_est_dd_perfect_shift = ch_dd_real_shift;
    ch_est_dd_perfect = ch_dd_real;
else
    ch_est_dd_perfect_shift = sqrt(num_subc/num_sym) * fft(ifft(ch_est_tf_perfect, [], 1), [], 2);    % dd domain
    ch_est_dd_perfect = fftshift(fftshift(ch_est_dd_perfect_shift, 1), 2);
end

% (2) evaluate dd channel (compare real channel and estimated channel)
if test_pilot
    ch_est_tf_test = rx_sym_tf ./ tx_sym_tf;     % tf domain
    ch_est_dd_test_shift = sqrt(num_subc/num_sym) * fft(ifft(ch_est_tf_test, [], 1), [], 2);    % dd domain
    ch_est_dd_test = fftshift(fftshift(ch_est_dd_test_shift, 1), 2);
    
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(ch_est_dd_perfect)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(ch_est_dd_test))
%     fprintf('Synch.: %3d    Avg. Pwr.: %7.2f\n', synch, sqrt(mean(abs(ch_est_dd_perfect(num_subc/2+1,:)).^2)))
%     fprintf('CFO.  : %3d    Avg. Pwr.: %7.2f\n', cfo, sqrt(mean(abs(ch_est_dd_perfect(:, num_sym/2+1)).^2)))
    pause
end

% (2) evaluate dd channel (compare real rx signal and estimated rx signal)
rx_sym_dd_est_shift = filter2(repmat(ch_est_dd_perfect(end:-1:1, end:-1:1), 2, 2), tx_sym_dd)/sqrt(num_subc*num_sym);
rx_sym_dd_est = fftshift(fftshift(rx_sym_dd_est_shift, 1), 2);
figure, plot(abs(rx_sym_dd(:)), '-b.'), hold on, plot(abs(rx_sym_dd_est(:)), ':r.'), hold off
figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd_est))
pause

% (3) generate block toeplitz dd channel matrix
new_tx = tx_sym_dd(:);
new_rx = rx_sym_dd(:);

len_ch = [num_subc num_sym];
len_ch_pad_head = len_ch/2;
len_ch_pad_tail = len_ch/2;
len_ch_pad_zero = len_ch;

base_ch = repmat(ch_est_dd_perfect(end:-1:1, end:-1:1), 2, 2);

% (3) generate 2d normal distribution channel for test
if test_gauss_ch
    sigma_ch = [1 8];
    mean_ch = len_ch .* [0.75 0.75];
    ch_norm_2d = zeros(size(ch_est_dd_perfect_shift));
    for idx = 1 : len_ch(1)
        ch_norm_2d_scale = 1/(sigma_ch(1)*sqrt(2*pi))*exp(1/(-2)*((idx-mean_ch(1))/sigma_ch(1)).^2);
        ch_norm_2d(idx, :) = ch_norm_2d_scale/(sigma_ch(2)*sqrt(2*pi))*exp(1/(-2)*(((1:len_ch(2))-mean_ch(2))/sigma_ch(2)).^2);
    end
    base_ch(1:len_ch(1), 1:len_ch(2)) = ch_norm_2d(end:-1:1, end:-1:1);
%     figure, mesh(1:num_subc, 1:num_sym, abs(ch_norm_2d))
%     pause
end

% (3) generate block toeplitz dd channel matrix
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

% (3) evaluate block toeplitz dd channel matrix
new_rx_est = new_ch * new_tx;
fprintf('Rank: %3d    Cond. Num.: %7.2f\n', rank(new_ch), cond(new_ch))
figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(ch_est_dd_perfect)), subplot(1, 2, 2), mesh(1:num_sym*num_subc, 1:num_subc*num_sym, abs(new_ch))
figure, plot(abs(new_rx(:)), '-b.'), hold on, plot(abs(new_rx_est(:)), ':r.'), hold off
figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(reshape(new_rx, num_subc, num_sym))), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(reshape(new_rx_est, num_subc, num_sym)))
figure, plot(abs(new_rx), '-b.'), hold on, plot(abs(new_rx-new_rx_est), ':r.'), hold off
pause

if test_signal
    pkt_error_ddeq = 0;
else
    % (3) equalize channel in dd domain
    rx_sym_dd_ddeq_vec = new_ch \ rx_sym_dd(:);
    % rx_sym_dd_ddeq_vec = pinv(new_ch) * rx_sym_dd(:);
    rx_sym_dd_ddeq = reshape(rx_sym_dd_ddeq_vec, num_subc, num_sym);
    
    % (3) observe symbols in dd domain
    rx_bit_ddeq = qamdemod(rx_sym_dd_ddeq(:), qam_size, 'UnitAveragePower', true);
    pkt_error_ddeq = symerr(tx_bit, rx_bit_ddeq);
    
    % plot
    figure, scatter(real(tx_sym_dd(:)), imag(tx_sym_dd(:)), 'bo'), hold on, scatter(real(rx_sym_dd_ddeq(:)), imag(rx_sym_dd_ddeq(:)), 'r.'), grid, axis([-2 2 -2 2]), hold off
    figure, plot(real(tx_sym_dd(:)), '-b.'), hold on, plot(real(rx_sym_dd_ddeq(:)), ':r.'), hold off
    figure, subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, abs(tx_sym_dd)), subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, abs(rx_sym_dd_ddeq))
end

assignin('base', 'ch_est_dd_perfect', ch_est_dd_perfect);
assignin('base', 'ch_est_dd_perfect', ch_est_dd_perfect);
assignin('base', 'new_ch', new_ch);

end
