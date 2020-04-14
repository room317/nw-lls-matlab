% test dd basic
% to observe the relation between rx signal and tx signal in various domain
%   - test_synch: synch position for test
%   - test_scope: print and plot results
%   - test_seed: synch position
% created: 2020.03.23
% modified:
%   - 2020.03.28:

function test_dd_basic(test_synch, test_scope, test_seed)

% set parameter
num_subc = 600; % 32;
num_sym = 14; % 32;
f_s = 15e3*num_subc;
qam_size = 16;
len_cp = 0; % num_fft/4; % 0;
carrier_freq_mhz = 4000;
velocity_kmh = 500; % 120;
idx_fading = 9; % 9: TDL-A, 12: TDL-D
delay_spread_rms_us = 0.1; %0.1e-6;
snr_db = 25;
test_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, snr_db, delay_spread_rms_us);

% check error
if test_synch > len_cp
    error('test_synch should not be greater than %d.\n', len_cp)
end

if test_seed >= 0
    rng(test_seed)
    
    % create a rayleigh fading channel object
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
    
    % dd tone generation
%     tx_pilot_position = [num_subc/2, num_sym/2];
%     tx_pilot_position = [300, 7];
    tx_pilot_position = [0, 0];
    tx_sym_dd = zeros(num_subc, num_sym);
    tx_sym_dd(tx_pilot_position(1)+1, tx_pilot_position(2)+1) = sqrt(num_subc*num_sym);
    
    % 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
    tx_sym_tf = sqrt(num_sym/num_subc) * fft(ifft(tx_sym_dd, [], 2), [], 1);
    
    % map to tf domain
    tx_sym_map_tf = fftshift(tx_sym_tf, 1);
    
    % ofdm modulate
    tx_ofdm_sym = sqrt(num_subc) * ifft(tx_sym_map_tf, [], 1);
    
    % add cp (cyclic postfix)
    tx_ofdm_sym_cp = [tx_ofdm_sym(end-len_cp+1:end, :); tx_ofdm_sym];
    
    % serialize
    tx_sig = tx_ofdm_sym_cp(:);
    
    % pass signal through channel
    rng(test_seed)
    tx_sig_faded = fading_ch(tx_sig);
    
    % add gaussian noise
    rx_sig = tx_sig_faded;
    
    % reshape
    rx_ofdm_sym_cp = reshape(rx_sig, num_subc+len_cp, num_sym);
    
    % remove cp (with synch: to observe impact of time shift)
    rx_ofdm_sym = rx_ofdm_sym_cp(test_synch+1:test_synch+num_subc, :);
    
    % ofdm demodulate the symbol
    rx_sym_map_tf = (1/sqrt(num_subc)) * fft(rx_ofdm_sym, [], 1);
    
    % demap
    rx_sym_tf = fftshift(rx_sym_map_tf, 1);
    
    % 2d inverse sfft
    rx_sym_dd = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf, [], 1), [], 2);
    
    % (1) estimate t channel
    tx_sig_norm = tx_sig ./ sqrt(abs(tx_sig).^2);
    tx_sig_comp = repmat(conj(tx_sig_norm(tx_pilot_position(1)+1:num_subc+len_cp:end, 1).'), num_subc+len_cp, 1);
    ch_est_t_imp_resp_1d = (1/sqrt(num_subc))*tx_sig_comp(:).*circshift(rx_sig, -tx_pilot_position(1), 1);
    ch_est_t_imp_resp = reshape(ch_est_t_imp_resp_1d, num_subc+len_cp, num_sym);
    
    % (2) estimate tf channel
    ch_est_tf_perfect = rx_sym_tf ./ tx_sym_tf;     % tf domain
    ch_est_tf_imp_resp = rx_sym_tf .* conj(tx_sym_tf);     % tf domain
    
    % (3) estimate dd channel
    ch_est_dd_perfect = sqrt(num_subc/num_sym) * fft(ifft(ch_est_tf_perfect, [], 1), [], 2);
    ch_est_dd_imp_resp_nocirc = rx_sym_dd;    % dd domain
    ch_est_dd_imp_resp_circ1d = circshift(ch_est_dd_imp_resp_nocirc, -tx_pilot_position(1), 1);
    ch_est_dd_imp_resp_circ2d = circshift(ch_est_dd_imp_resp_circ1d, -tx_pilot_position(2), 2);
    ch_est_dd_imp_resp = ch_est_dd_imp_resp_circ2d;
end

if test_scope
    % (1) test t channel
    test_rx_sig_t_est = zeros(num_subc+len_cp,num_sym);
    test_rx_sig_t_est_prev = zeros(num_subc+len_cp,1);
    for idx = 1 : num_sym
        test_rx_sig_t_est_tail = conv(ch_est_t_imp_resp(:, idx), tx_ofdm_sym_cp(:, idx));
        test_rx_sig_t_est(:, idx) = test_rx_sig_t_est_prev+test_rx_sig_t_est_tail(1:num_subc+len_cp, :);
        test_rx_sig_t_est_prev = [test_rx_sig_t_est_tail(num_subc+len_cp+1:end, :); 0];
    end
    figure
    subplot(2, 1, 1), plot(real(tx_sig_faded(:)), '-b.'), hold on, plot(real(test_rx_sig_t_est(:)), ':r.'), hold off, title('real(rx sig) and real(rx sig est) in t-domain')
    subplot(2, 1, 2), plot(imag(tx_sig_faded(:)), '-b.'), hold on, plot(imag(test_rx_sig_t_est(:)), ':r.'), hold off, title('imag(rx sig) and real(rx sig est) in t-domain')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc+len_cp, real(ch_est_t_imp_resp)), title('real(t-channel impulse response)')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc+len_cp, imag(ch_est_t_imp_resp)), title('imag(t-channel impulse response)')
    pause
    
    % (2) test tf channel
    figure
    subplot(2, 1, 1), plot(real(ch_est_tf_perfect(:)), '-b.'), hold on, plot(real(ch_est_tf_imp_resp(:)), ':r.'), hold off, title('real(tf-channel) and real(tf-channel impulse response)')
    subplot(2, 1, 2), plot(imag(ch_est_tf_perfect(:)), '-b.'), hold on, plot(imag(ch_est_tf_imp_resp(:)), ':r.'), hold off, title('imag(rx channel) and imag(tf-channel impulse response)')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(ch_est_tf_imp_resp)), title('real(tf-channel)')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(ch_est_tf_imp_resp)), title('imag(tf-channel)')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(fftshift(ch_est_tf_imp_resp, 1))), title('real(fftshift(tf-channel))')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(fftshift(ch_est_tf_imp_resp, 1))), title('imag(fftshift(tf-channel))')
    pause
    
    % (3) test dd channel
    test_rx_sig_dd_est_tail = (1/sqrt(num_subc*num_sym))*conv2(repmat(ch_est_dd_imp_resp,2,2), tx_sym_dd, 'full');
    test_rx_sig_dd_est = test_rx_sig_dd_est_tail(num_subc+1:2*num_subc, num_sym+1:2*num_sym);
    figure
    subplot(2, 1, 1), plot(real(ch_est_dd_perfect(:)), '-b.'), hold on, plot(real(ch_est_dd_imp_resp(:)), ':r.'), hold off, title('real(dd-channel) and real(dd-channel impulse response)')
    subplot(2, 1, 2), plot(imag(ch_est_dd_perfect(:)), '-b.'), hold on, plot(imag(ch_est_dd_imp_resp(:)), ':r.'), hold off, title('imag(dd-channel) and imag(dd-channel impulse response)')
    figure
    subplot(2, 1, 1), plot(real(rx_sym_dd(:)), '-b.'), hold on, plot(real(test_rx_sig_dd_est(:)), ':r.'), hold off, title('real(rx sig) and real(rx sig est) in dd-domain')
    subplot(2, 1, 2), plot(imag(rx_sym_dd(:)), '-b.'), hold on, plot(imag(test_rx_sig_dd_est(:)), ':r.'), hold off, title('imag(rx sig) and imag(rx sig est) in dd-domain')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(ch_est_dd_imp_resp)), title('real(dd-channel impulse response)')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(ch_est_dd_imp_resp)), title('imag(dd-channel impulse response)')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(fftshift(ch_est_dd_perfect, 2))), title('real(fftshift(dd-channel transform))')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(fftshift(ch_est_dd_perfect, 2))), title('imag(fftshift(dd-channel transform))')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(rx_sym_dd)), title('real(fftshift(dd-channel transform))')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(rx_sym_dd)), title('imag(fftshift(dd-channel transform))')
    pause
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

% generate bit stream
tx_bit = randi([0 qam_size-1], num_subc*num_sym, 1);

% modulate bit stream
tx_sym = qammod(tx_bit, qam_size, 'UnitAveragePower', true);

% reshape
tx_sym_dd = reshape(tx_sym, num_subc, []);

% tx_sym_dd = zeros(num_subc, num_sym);
% tx_sym_dd(tx_pilot_position(1)+1, tx_pilot_position(2)+1) = sqrt(num_subc*num_sym);

% 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
tx_sym_tf = sqrt(num_sym/num_subc) * fft(ifft(tx_sym_dd, [], 2), [], 1);

% map to tf domain
tx_sym_map_tf = fftshift(tx_sym_tf, 1);

% ofdm modulate
tx_ofdm_sym = sqrt(num_subc) * ifft(tx_sym_map_tf, [], 1);

% add cp (cyclic postfix)
tx_ofdm_sym_cp = [tx_ofdm_sym(end-len_cp+1:end, :); tx_ofdm_sym];

% serialize
tx_sig = tx_ofdm_sym_cp(:);

% pass signal through channel
if test_seed >= 0
    rng(test_seed)
end
tx_sig_faded = fading_ch(tx_sig);

% add gaussian noise
rx_sig = tx_sig_faded;

% reshape
rx_ofdm_sym_cp = reshape(rx_sig, num_subc+len_cp, num_sym);

% remove cp (with synch: to observe impact of time shift)
rx_ofdm_sym = rx_ofdm_sym_cp(test_synch+1:test_synch+num_subc, :);

% ofdm demodulate the symbol
rx_sym_map_tf = (1/sqrt(num_subc)) * fft(rx_ofdm_sym, [], 1);

% demap
rx_sym_tf = fftshift(rx_sym_map_tf, 1);

% 2d inverse sfft
rx_sym_dd = sqrt(num_subc/num_sym) * fft(ifft(rx_sym_tf, [], 1), [], 2);

if test_scope
    % (1) test t channel
    test_rx_sig_t_est = zeros(num_subc+len_cp,num_sym);
    test_rx_sig_t_est_prev = zeros(num_subc+len_cp,1);
    for idx = 1 : num_sym
        test_rx_sig_t_est_tail = conv(ch_est_t_imp_resp(:, idx), tx_ofdm_sym_cp(:, idx));
        test_rx_sig_t_est(:, idx) = test_rx_sig_t_est_prev+test_rx_sig_t_est_tail(1:num_subc+len_cp, :);
        test_rx_sig_t_est_prev = [test_rx_sig_t_est_tail(num_subc+len_cp+1:end, :); 0];
    end
    figure
    subplot(2, 1, 1), plot(real(rx_ofdm_sym_cp(:)), '-b.'), hold on, plot(real(test_rx_sig_t_est(:)), ':r.'), hold off, title('real(rx sig) and real(rx sig est) in t-domain')
    subplot(2, 1, 2), plot(imag(rx_ofdm_sym_cp(:)), '-b.'), hold on, plot(imag(test_rx_sig_t_est(:)), ':r.'), hold off, title('imag(rx sig) and imag(rx sig est) in t-domain')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc+len_cp, real(rx_ofdm_sym_cp)), title('real(rx sig) in t-domain')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc+len_cp, imag(rx_ofdm_sym_cp)), title('imag(rx sig) in t-domain')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc+len_cp, real(test_rx_sig_t_est)), title('real(rx sig est) in t-domain')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc+len_cp, imag(test_rx_sig_t_est)), title('imag(rx sig est) in t-domain')
    fprintf('MSE of received signal in t-domain: %6.4f\n', sqrt(mean(abs(rx_ofdm_sym_cp(:)-test_rx_sig_t_est(:)).^2)))
    pause
    
    % (2) test tf channel
    test_rx_sig_tf_est = ch_est_tf_perfect .* tx_sym_tf;
    figure
    subplot(2, 1, 1), plot(real(rx_sym_tf(:)), '-b.'), hold on, plot(real(test_rx_sig_tf_est(:)), ':r.'), hold off, title('real(rx sig) and real(rx sig est) in tf-domain')
    subplot(2, 1, 2), plot(imag(rx_sym_tf(:)), '-b.'), hold on, plot(imag(test_rx_sig_tf_est(:)), ':r.'), hold off, title('imag(rx sig) and imag(rx sig est) in tf-domain)')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(rx_sym_tf)), title('real(rx sig) in tf-domain')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(rx_sym_tf)), title('imag(rx sig) in tf-domain')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(test_rx_sig_tf_est)), title('real(rx sig est) in tf-domain')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(test_rx_sig_tf_est)), title('imag(rx sig est) in tf-domain')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(fftshift(rx_sym_tf, 1))), title('real(fftshift(rx sig)) in tf-domain')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(fftshift(rx_sym_tf, 1))), title('imag(fftshift(rx sig)) in tf-domain')
    fprintf('MSE of received signal in tf-domain: %6.4f\n', sqrt(mean(abs(rx_sym_tf(:)-test_rx_sig_tf_est(:)).^2)))
    pause
    
    % (3) test dd channel
    test_rx_sig_dd_est_tail = (1/sqrt(num_subc*num_sym))*conv2(repmat(ch_est_dd_imp_resp,2,2), tx_sym_dd, 'full');
    test_rx_sig_dd_est = test_rx_sig_dd_est_tail(num_subc+1:2*num_subc, num_sym+1:2*num_sym);
    figure
    subplot(2, 1, 1), plot(real(rx_sym_dd(:)), '-b.'), hold on, plot(real(test_rx_sig_dd_est(:)), ':r.'), hold off, title('real(rx sig) and real(rx sig est) in dd-domain')
    subplot(2, 1, 2), plot(imag(rx_sym_dd(:)), '-b.'), hold on, plot(imag(test_rx_sig_dd_est(:)), ':r.'), hold off, title('imag(rx sig) and imag(rx sig est) in dd-domain')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(rx_sym_dd)), title('real(rx sig) in dd-domain')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(rx_sym_dd)), title('imag(rx sig) in dd-domain')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(test_rx_sig_dd_est)), title('real(rx sig est) in dd-domain')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(test_rx_sig_dd_est)), title('imag(rx sig est) in dd-domain')
    figure
    subplot(1, 2, 1), mesh(1:num_sym, 1:num_subc, real(fftshift(rx_sym_dd, 2))), title('real(fftshift(rx sig)) in dd-domain')
    subplot(1, 2, 2), mesh(1:num_sym, 1:num_subc, imag(fftshift(rx_sym_dd, 2))), title('imag(fftshift(rx sig)) in dd-domain')
    fprintf('MSE of received signal in dd-domain: %6.4f\n', sqrt(mean(abs(rx_sym_dd(:)-test_rx_sig_dd_est(:)).^2)))
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
    assignin('base', 'test_rx_sig_t_est', test_rx_sig_t_est);
    assignin('base', 'test_rx_sig_tf_est', test_rx_sig_tf_est);
    assignin('base', 'test_rx_sig_dd_est', test_rx_sig_dd_est);
end

end
