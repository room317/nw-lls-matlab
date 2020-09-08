% test dd pilot_r4 (real channel testing)
%   - snr_db: snr in db scale
%   - test_synch: synch position
%   - test_scope: print and plot results
%   - test_seed: random seed for channel
%   - test_chest: channel estimation option {'dd_tone', 'perfect', 'real'}
%   - test_cheq: channel eq. option {'tf_zf', 'tf_mmse', 'dd'}
% created: 2020.03.06
% modified:
%   - 2020.03.02: prefix/postfix option
%   - 2020.03.03: rician channel added
%   - 2020.03.06: pilot estimation added
%   - 2020.04.29: removed unused options
%   - 2020.09.07: regenerate full-tap real channel (dd- and tf-domain)

function [qam_error, num_qam_usr, ch_est_rmse, papr_db] = test_dd_pilot_r4(test_metric, test_scs_khz, test_bw_mhz, test_num_slot, test_synch, test_scope, test_seed, test_ch_rmse, test_papr, test_chest, test_cheq)

%% parameters

% set subcarrier spacing parameters
scs_khz_list = [15 30 60];
idx_scs = find(scs_khz_list == test_scs_khz, 1);
if isempty(idx_scs)
    error('Subcarrier spacing must be one of these: {15, 30, 60}')
end

% set bandwidth parameters
bw_mhz_list = [5 10 15 20 25 30 40 50 60 80 90 100];
idx_bw = find(bw_mhz_list == test_bw_mhz, 1);
if isempty(idx_bw)
    error('Bandwidth must be one of these: {5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 90, 100}')
end

% set resource block table
% ref.: http://howltestuffworks.blogspot.com/2019/11/5g-nr-resource-blocks.html
% ref.: https://www.sharetechnote.com/html/5G/5G_FR_Bandwidth.html
%     bw (mhz)     |   5 |  10 |  15 |  20 |  25 |  30 |  40 |  50 |  60 |  80 |  90 | 100
%     -------------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+----
%     scs 15 (khz) |  25 |  52 |  79 | 106 | 133 | 160 | 216 | 270 |   0 |   0 |   0 |   0
%     scs 30 (khz) |  11 |  24 |  38 |  51 |  65 |  78 | 106 | 133 | 162 | 217 | 245 | 273
%     scs 60 (khz) |   0 |  11 |  18 |  24 |  31 |  38 |  51 |  65 |  79 | 107 | 121 | 135
rb_tbl = [ ...
    25   52   79  106  133  160  216  270    0    0    0    0
    11   24   38   51   65   78  106  133  162  217  245  273
     0   11   18   24   31   38   51   65   79  107  121  135];

% set numerical parameters
num_rb = rb_tbl(idx_scs, idx_bw);                   % number of rbs
num_subc_rb = 12;                                   % number of subcarriers per rb
num_subc_bw = num_rb*num_subc_rb;                   % number of subcarriers in bandwidth
num_fft = 2^ceil(log2(num_subc_bw));                % fft size
sample_rate = test_scs_khz*1e3*num_fft;             % sampling rate (hz)
num_cp = round(num_fft*144/2048);                   % number of samples in cyclic prefix
num_ofdmsym_slot = 14;                              % number of ofdm symbols per slot
num_ofdmsym = num_ofdmsym_slot*test_num_slot;       % number of ofdm symbols
t_ofdmsym = (num_fft+num_cp)/sample_rate;           % ofdm symbol time

% set user parameters (full span, prb_size: 12*14)
num_rb_usr = num_rb; % 3                            % number of resource blocks per user
num_slot_usr = test_num_slot;                       % number of slots per user
num_subc_usr = num_rb_usr*num_subc_rb;              % number of subcarriers per user
num_ofdmsym_usr = num_ofdmsym_slot*num_slot_usr;    % number of ofdm symbols per user (slot-based)

% set dd parameters
num_delay_usr = num_subc_usr;                       % fixed
num_doppler_usr = num_ofdmsym_usr;                  % fixed
num_delay_pilot = double(strcmp(test_chest, 'dd_tone'))*round(num_delay_usr*0.05)*2;      % fixed (10% of available delay grids)
num_doppler_pilot = double(strcmp(test_chest, 'dd_tone'))*num_doppler_usr;                % fixed (100% of available doppler grids)
num_delay_data = num_delay_usr-num_delay_pilot;
num_doppler_data = num_doppler_usr-num_doppler_pilot;
num_delay_guard = double(num_delay_usr~=num_delay_pilot)*round(num_delay_pilot*0.2);            % 20% of pilot delay grids
num_doppler_guard = double(num_doppler_usr~=num_doppler_pilot)*round(num_doppler_pilot*0.2);    % 20% of pilot doppler grids

% set test parameter
snr_db = test_metric;
qam_size = 16;
cfo_norm = 0;               % normalized cfo. cfo = cfo_norm/t_sym;
carrier_freq_mhz = 4000;
velocity_kmh = 500;
idx_fading = 13;            % 9: TDL-A, 10: TDL-B, 11: TDL-C, 12: TDL-D, 13: TDL-E
delay_spread_rms_us = 0.01;  % 0.1e-6;
test_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, snr_db, delay_spread_rms_us);

% create a rayleigh fading channel object
if test_seed >= 0
    rng(test_seed)
end
switch idx_fading
    case {12, 13}
        fading_ch = comm.RicianChannel(...
            'SampleRate', sample_rate,...
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
            'SampleRate', sample_rate, ...
            'PathDelays', test_ch.path_delays, ...
            'AveragePathGains', test_ch.average_path_gains, ...
            'NormalizePathGains', true, ...
            'MaximumDopplerShift', test_ch.maximum_doppler_shift, ...
            'PathGainsOutputPort', true, ...
            'DopplerSpectrum', test_ch.doppler_spectrum);
end

% calculate noise variance
noise_var = (num_subc_usr/num_fft)*10^((-0.1)*snr_db);

%% transmitter

% calculate num. qam symbols per packet
if strcmp(test_chest, 'dd_tone')
    num_qam_usr = (num_delay_usr*num_doppler_usr)-(num_delay_pilot*num_doppler_pilot);
else
    num_qam_usr = num_delay_usr*num_doppler_usr;
end

% generate bit stream
tx_bit_usr = randi([0 qam_size-1], num_qam_usr, 1);

% modulate bit stream
tx_sym_data_usr = qammod(tx_bit_usr, qam_size, 'UnitAveragePower', true);

% map symbols
if strcmp(test_chest, 'dd_tone')
    % reshape data symbols
    tx_sym_data1 = reshape(tx_sym_data_usr(1:num_delay_pilot*num_doppler_data), num_delay_pilot, []);
    tx_sym_data2 = reshape(tx_sym_data_usr(num_delay_pilot*num_doppler_data+1:end), num_delay_data, []);
    
    % set index
    idx_delay_ctr = floor(num_delay_usr/2)+1;
    idx_doppler_ctr = floor(num_doppler_usr/2)+1;
    idx_delay_pilot_ctr = floor(num_delay_pilot/2)+1;
    idx_doppler_pilot_ctr = floor(num_doppler_pilot/2)+1;
    
    % generate pilot symbols with guard
    tx_sym_pilot = zeros(num_delay_pilot, num_doppler_pilot);
    tx_sym_pilot(idx_delay_pilot_ctr, idx_doppler_pilot_ctr) = sqrt(num_delay_pilot*num_doppler_pilot);
    
    % map data and pilot
    tx_sym_dd_base = [...
        tx_sym_pilot, tx_sym_data1;
        tx_sym_data2];
    tx_sym_dd_usr = circshift(tx_sym_dd_base, [idx_delay_ctr-idx_delay_pilot_ctr, idx_doppler_ctr-idx_doppler_pilot_ctr]);
else
    % reshape data symbols
    tx_sym_data_usr = reshape(tx_sym_data_usr, num_delay_usr, []);
    
    % map data and pilot
    tx_sym_dd_usr = tx_sym_data_usr;
end

% 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
tx_sym_tf_usr = sqrt(num_ofdmsym_usr/num_subc_bw)*fft(ifft(tx_sym_dd_usr, [], 2), [], 1);

% map user block to whole resource blocks (temp)
tx_sym_tf = zeros(num_subc_bw, num_ofdmsym);
tx_sym_tf(1:size(tx_sym_tf_usr, 1), 1:size(tx_sym_tf_usr, 2)) = tx_sym_tf_usr;

% map to fft range
tx_sym_nfft_shift = zeros(num_fft, num_ofdmsym);
tx_sym_nfft_shift((num_fft/2)-(num_subc_bw/2)+1:(num_fft/2)+(num_subc_bw/2), :) = tx_sym_tf;
tx_sym_nfft = fftshift(tx_sym_nfft_shift, 1);

% ofdm modulate
tx_ofdmsym = sqrt(num_fft)*ifft(tx_sym_nfft, [], 1);

% add cp (cyclic prefix)
tx_ofdmsym_cp = tx_ofdmsym([num_fft-num_cp+1:num_fft, 1:num_fft], :);

% serialize
tx_sig = tx_ofdmsym_cp(:);

%% channel

% pass signal through channel
[tx_sig_faded, ch_path_gain] = fading_ch(tx_sig);       % ch_path_gain: normally constant per path, vary when doppler exists

% add gaussian noise
% rx_sig = awgn(tx_sig_faded, snr_db, 'measured');
rx_sig = tx_sig_faded;

% regenerate real channel
if test_scope || test_ch_rmse || strcmp(test_chest, 'real')
    [ch_real_onetap_tf, ch_real_mat_tf, ch_real_eff_tf, ch_real_eff_dd] = gen_real_ch(fading_ch, ch_path_gain, num_fft, num_cp, num_subc_bw, num_ofdmsym);
else
    ch_real_onetap_tf = [];
    ch_real_mat_tf = [];
    ch_real_eff_tf = [];
    ch_real_eff_dd = [];
end

%% receiver

% compensate cfo (to observe impact of frequency shift)
cfo_sample = 0:length(rx_sig)-1;
cfo_vec = cfo_sample.'*cfo_norm/(num_fft+num_cp);
rx_sig_cfo = rx_sig(:).*exp(-1i*2*pi*cfo_vec);

% synchronization
rx_sig_synch = circshift(rx_sig_cfo, -num_cp-test_synch);

% reshape
rx_ofdmsym_cp = reshape(rx_sig_synch, num_fft+num_cp, []);

% remove cp (with synch: to observe impact of time shift)
rx_ofdmsym = rx_ofdmsym_cp(1:num_fft, :);

% ofdm demodulate the symbol
rx_sym_nfft = (1/sqrt(num_fft))*fft(rx_ofdmsym, [], 1);

% demap resource plane
rx_sym_nfft_shift = fftshift(rx_sym_nfft, 1);
rx_sym_tf = rx_sym_nfft_shift((num_fft/2)-(num_subc_bw/2)+1:(num_fft/2)+(num_subc_bw/2), :);

% demap user block
rx_sym_tf_usr = rx_sym_tf(1:num_subc_usr, 1:num_ofdmsym_usr);

% 2d inverse sfft
rx_sym_dd_usr = sqrt(num_subc_bw/num_ofdmsym_usr)*fft(ifft(rx_sym_tf_usr, [], 1), [], 2);

% estimate channel
if strcmp(test_chest, 'dd_tone')
    % demap pilot and remove guard
    rx_sym_dd_base = circshift(rx_sym_dd_usr, [-(idx_delay_ctr-idx_delay_pilot_ctr), -(idx_doppler_ctr-idx_doppler_pilot_ctr)]);
    rx_sym_pilot = rx_sym_dd_base(num_delay_guard+1:num_delay_pilot-num_delay_guard, num_doppler_guard+1:num_doppler_pilot-num_doppler_guard);
    
    % estimate channel
    ch_est_dd_base = zeros(num_delay_usr, num_doppler_usr);
    ch_est_dd_base(num_delay_guard+1:num_delay_pilot-num_delay_guard, num_doppler_guard+1:num_doppler_pilot-num_doppler_guard) = ...
        rx_sym_pilot*sqrt((num_delay_usr*num_doppler_usr)/(num_delay_pilot*num_doppler_pilot));
    ch_est_dd = circshift(ch_est_dd_base, [-idx_delay_pilot_ctr+1, -idx_doppler_pilot_ctr+1]);
    
    % transform channel
    ch_est_tf = sqrt(num_ofdmsym_usr/num_subc_bw)*fft(ifft(ch_est_dd, [], 2), [], 1);
elseif strcmp(test_chest, 'perfect')
    % estimate tf channel
    ch_est_tf = rx_sym_tf_usr./tx_sym_tf_usr;     % tf domain
    
    % estimate dd channel
    ch_est_dd = sqrt(num_subc_bw/num_ofdmsym_usr)*fft(ifft(ch_est_tf, [], 1), [], 2);
elseif strcmp(test_chest, 'real')
    % get real tf channel
    ch_est_tf = ch_real_onetap_tf;
    
    % estimate dd channel
    ch_est_dd = sqrt(num_subc_bw/num_ofdmsym_usr)*fft(ifft(ch_est_tf, [], 1), [], 2);
else
    error('test_chest value must be one of these: {''dd_tone'', ''perfect'', ''real''}')
end
    
% generate block circular channel matrix (to check channel error)
if test_scope || test_ch_rmse || strcmp(test_cheq, 'dd')
    ch_eff_dd = gen_eff_ch(ch_est_dd);
else
    ch_eff_dd = [];
end

% equalize channel
if strcmp(test_cheq, 'tf_zf')
    % equalize channel in tf domain
    if strcmp(test_chest, 'real')
        % full-tap equalization (matrix inversion)
        rx_sym_tf_eq_vec = ch_real_eff_tf\rx_sym_tf_usr(:);
        rx_sym_tf_eq = reshape(rx_sym_tf_eq_vec, num_subc_usr, num_ofdmsym_usr);
    else
        % one-tap equalization (element-wise)
        rx_sym_tf_eq = rx_sym_tf_usr./ch_est_tf;
    end
    
    % observe symbols in dd domain
    rx_sym_dd_eq = sqrt(num_subc_bw/num_ofdmsym_usr)*fft(ifft(rx_sym_tf_eq, [], 1), [], 2);
elseif strcmp(test_cheq, 'tf_mmse')
    % equalize channel in tf domain
    if strcmp(test_chest, 'real')
        % full-tap equalization (matrix inversion)
        ch_real_eff_tf_mmse = ch_real_eff_tf'/(ch_real_eff_tf*ch_real_eff_tf'+noise_var*eye(num_subc_usr*num_ofdmsym_usr));
        rx_sym_tf_eq_vec = ch_real_eff_tf_mmse*rx_sym_tf_usr(:);
        rx_sym_tf_eq = reshape(rx_sym_tf_eq_vec, num_subc_usr, num_ofdmsym_usr);
    else
        % one-tap equalization (element-wise)
        ch_est_tf_mmse = conj(ch_est_tf)./(noise_var+abs(ch_est_tf).^2);
        rx_sym_tf_eq = rx_sym_tf_usr.*ch_est_tf_mmse;
    end
    
    % observe symbols in dd domain
    rx_sym_dd_eq = sqrt(num_subc_bw/num_ofdmsym_usr)*fft(ifft(rx_sym_tf_eq, [], 1), [], 2);
elseif strcmp(test_cheq, 'dd_zf')
    % equalize channel in dd domain
    if strcmp(test_chest, 'real')
        % full-tap equalization (matrix inversion)
        rx_sym_dd_eq_vec = ch_real_eff_dd\rx_sym_dd_usr(:);
        rx_sym_dd_eq = reshape(rx_sym_dd_eq_vec, num_delay_usr, num_doppler_usr);
    else
        % one-tap equalization (matrix inversion)
        rx_sym_dd_eq_vec = ch_eff_dd\rx_sym_dd_usr(:);
        rx_sym_dd_eq = reshape(rx_sym_dd_eq_vec, num_delay_usr, num_doppler_usr);
    end
elseif strcmp(test_cheq, 'dd_mmse')
    % equalize channel in dd domain
    if strcmp(test_chest, 'real')
        % full-tap equalization (matrix inversion)
        ch_real_eff_dd_mmse = ch_real_eff_dd'/(ch_real_eff_dd*ch_real_eff_dd'+noise_var*eye(num_delay_usr*num_doppler_usr));
        rx_sym_dd_eq_vec = ch_real_eff_dd_mmse*rx_sym_dd_usr(:);
        rx_sym_dd_eq = reshape(rx_sym_dd_eq_vec, num_delay_usr, num_doppler_usr);
    else
        % one-tap equalization (matrix inversion)
        ch_eff_dd_mmse = ch_eff_dd'/(ch_eff_dd*ch_eff_dd'+noise_var*eye(num_delay_usr*num_doppler_usr));
        rx_sym_dd_eq_vec = ch_eff_dd_mmse*rx_sym_dd_usr(:);
        rx_sym_dd_eq = reshape(rx_sym_dd_eq_vec, num_delay_usr, num_doppler_usr);
    end
else
    error('test_chest value must be one of these: {tf_zf, tf_mmse, dd_zf, dd_mmse}')
end

% demap data symbols
if strcmp(test_chest, 'dd_tone')
    rx_sym_dd_base = circshift(rx_sym_dd_eq, [-(idx_delay_ctr-idx_delay_pilot_ctr), -(idx_doppler_ctr-idx_doppler_pilot_ctr)]);
    rx_sym_data1 = rx_sym_dd_base(1:num_delay_pilot, num_doppler_pilot+1:end);
    rx_sym_data2 = rx_sym_dd_base(num_delay_pilot+1:end, :);
    rx_sym_data_usr = [rx_sym_data1(:); rx_sym_data2(:)];
else
    rx_sym_data_usr = rx_sym_dd_eq(:);
end

% demodulate qam symbols
rx_bit_usr = qamdemod(rx_sym_data_usr, qam_size, 'UnitAveragePower', true);

% calculate bit errors
if isempty(tx_bit_usr)
    qam_error = 0;
else
    qam_error = symerr(tx_bit_usr, rx_bit_usr);
end

% calculate channel rmse
if test_ch_rmse
    ch_est_rmse = sqrt(mean(abs(ch_real_eff_dd-ch_eff_dd).^2, 'all'));
else
    ch_est_rmse = [];
end

% calculate papr
if test_papr
    papr_db = 10*log10(max(abs(tx_sig).^2)/mean(abs(tx_sig).^2));
else
    papr_db = [];
end

% test scope
if test_scope
    % check bits
    figure
    plot(tx_bit_usr, '-b.'), hold on, plot(rx_bit_usr, ':r.'), hold off, grid minor
    
    % check qam
    figure 
    subplot(2, 1, 1), plot(real(tx_sym_data_usr(:)), '-b.'), hold on, plot(real(rx_sym_data_usr(:)), ':r.'), hold off, grid minor
    subplot(2, 1, 2), plot(imag(tx_sym_data_usr(:)), '-b.'), hold on, plot(imag(rx_sym_data_usr(:)), ':r.'), hold off, grid minor
    
    % check dd sym
    figure
    subplot(1, 3, 1), mesh(1:num_doppler_usr, 1:num_delay_usr, abs(tx_sym_dd_usr))
    subplot(1, 3, 2), mesh(1:num_doppler_usr, 1:num_delay_usr, abs(rx_sym_dd_usr))
    subplot(1, 3, 3), mesh(1:num_doppler_usr, 1:num_delay_usr, abs(rx_sym_dd_eq))
    
    % check tf sym
    figure
    subplot(1, 2, 1), mesh(1:num_ofdmsym_usr, 1:num_subc_usr, abs(tx_sym_tf_usr))
    subplot(1, 2, 2), mesh(1:num_ofdmsym_usr, 1:num_subc_usr, abs(rx_sym_tf_usr))
    
    % check fft sym
    figure
    subplot(1, 2, 1), mesh(1:num_ofdmsym, 1:num_fft, abs(tx_sym_nfft))
    subplot(1, 2, 2), mesh(1:num_ofdmsym, 1:num_fft, abs(rx_sym_nfft))
    
    % check time sig
    figure
    subplot(1, 2, 1), mesh(1:num_ofdmsym, 1:num_fft, abs(tx_ofdmsym))
    subplot(1, 2, 2), mesh(1:num_ofdmsym, 1:num_fft, abs(rx_ofdmsym))
    
    % check time sig with cp
    figure
    subplot(1, 2, 1), mesh(1:num_ofdmsym, 1:num_fft+num_cp, abs(tx_ofdmsym_cp))
    subplot(1, 2, 2), mesh(1:num_ofdmsym, 1:num_fft+num_cp, abs(rx_ofdmsym_cp))
    
    % check one-tap tf channel
    figure
    subplot(1, 2, 1), mesh((0:num_ofdmsym_usr-1)*t_ofdmsym*1e3, (0:num_subc_usr-1)*test_scs_khz*1e3*1e-6, abs(ch_real_onetap_tf)), xlabel('ofdm symbols (ms)'), ylabel('subcarriers (khz)'), title('real channel (one-tap)')
    subplot(1, 2, 2), mesh((0:num_ofdmsym_usr-1)*t_ofdmsym*1e3, (0:num_subc_usr-1)*test_scs_khz*1e3*1e-6, abs(ch_est_tf)), xlabel('ofdm symbols (ms)'), ylabel('subcarriers (khz)'), title('estimated channel (one-tap)')
    
    % check one-tap dd channel
    ch_real_onetap_dd = sqrt(num_subc_bw/num_ofdmsym_usr)*fft(ifft(ch_real_onetap_tf, [], 1), [], 2);
    figure
    subplot(1, 2, 1), mesh(((0:num_doppler_usr-1)-floor(num_doppler_usr/2))/(t_ofdmsym*num_doppler_usr*1e3), ((0:num_delay_usr-1)-floor(num_delay_usr/2))/(num_delay_usr*test_scs_khz*1e3*1e-6), abs(fftshift(fftshift(ch_real_onetap_dd, 1), 2))), xlabel('doppler (khz)'), ylabel('delay (us)'), title('real channel (one-tap)')
    subplot(1, 2, 2), mesh(((0:num_doppler_usr-1)-floor(num_doppler_usr/2))/(t_ofdmsym*num_doppler_usr*1e3), ((0:num_delay_usr-1)-floor(num_delay_usr/2))/(num_delay_usr*test_scs_khz*1e3*1e-6), abs(fftshift(fftshift(ch_est_dd, 1), 2))), xlabel('doppler (khz)'), ylabel('delay (us)'), title('estimated channel (one-tap)')
    
    % check effective tf channel
    rx_sym_tf_usr_regen = ch_real_eff_tf*tx_sym_tf_usr(:);
    fprintf('rx symbol rmse in tf: %10.4f\n', sqrt(mean(abs(rx_sym_tf_usr(:)-rx_sym_tf_usr_regen).^2)))
    figure
    subplot(2, 1, 1), plot(real(rx_sym_tf_usr(:)), '-b.'), hold on, plot(real(rx_sym_tf_usr_regen), ':r.'), hold off
    subplot(2, 1, 2), plot(imag(rx_sym_tf_usr(:)), '-b.'), hold on, plot(imag(rx_sym_tf_usr_regen), ':r.'), hold off
    
    % check effective dd channel
    rx_sym_dd_usr_regen = ch_real_eff_dd*tx_sym_dd_usr(:);
    fprintf('rx symbol rmse in dd: %10.4f\n', sqrt(mean(abs(rx_sym_dd_usr(:)-rx_sym_dd_usr_regen).^2)))
    figure
    subplot(2, 1, 1), plot(real(rx_sym_dd_usr(:)), '-b.'), hold on, plot(real(rx_sym_dd_usr_regen), ':r.'), hold off
    subplot(2, 1, 2), plot(imag(rx_sym_dd_usr(:)), '-b.'), hold on, plot(imag(rx_sym_dd_usr_regen), ':r.'), hold off
end

% dump
if test_scope
    assignin('base', 'num_fft', num_fft);
    assignin('base', 'sample_rate', sample_rate);
    assignin('base', 'num_subc_bw', num_subc_bw);
    assignin('base', 'num_cp', num_cp);
    assignin('base', 'num_ofdmsym', num_ofdmsym);
    assignin('base', 'num_subc_usr', num_subc_usr);
    assignin('base', 'num_ofdmsym_usr', num_ofdmsym_usr);
    assignin('base', 'num_delay_usr', num_delay_usr);
    assignin('base', 'num_doppler_usr', num_doppler_usr);
    assignin('base', 'num_delay_pilot', num_delay_pilot);
    assignin('base', 'num_doppler_pilot', num_doppler_pilot);
    assignin('base', 'num_delay_data', num_delay_data);
    assignin('base', 'num_doppler_data', num_doppler_data);
    assignin('base', 'num_delay_guard', num_delay_guard);
    assignin('base', 'num_doppler_guard', num_doppler_guard);
    
    assignin('base', 'tx_bit_usr', tx_bit_usr);
    assignin('base', 'rx_bit_usr', rx_bit_usr);
    assignin('base', 'tx_sym_data_usr', tx_sym_data_usr);
    assignin('base', 'rx_sym_data_usr', rx_sym_data_usr);
    assignin('base', 'tx_sym_dd_usr', tx_sym_dd_usr);
    assignin('base', 'rx_sym_dd_usr', rx_sym_dd_usr);
    assignin('base', 'rx_sym_dd_eq', rx_sym_dd_eq);
    assignin('base', 'tx_sym_tf_usr', tx_sym_tf_usr);
    assignin('base', 'rx_sym_tf_usr', rx_sym_tf_usr);
    assignin('base', 'tx_sym_tf', tx_sym_tf);
    assignin('base', 'rx_sym_tf', rx_sym_tf);
    assignin('base', 'tx_sym_nfft', tx_sym_nfft);
    assignin('base', 'rx_sym_nfft', rx_sym_nfft);
    assignin('base', 'tx_ofdmsym', tx_ofdmsym);
    assignin('base', 'rx_ofdmsym', rx_ofdmsym);
    assignin('base', 'tx_ofdmsym_cp', tx_ofdmsym_cp);
    assignin('base', 'rx_ofdmsym_cp', rx_ofdmsym_cp);
    assignin('base', 'tx_sig', tx_sig);
    assignin('base', 'tx_sig_faded', tx_sig_faded);
    assignin('base', 'rx_sig', rx_sig);
    assignin('base', 'rx_sig_cfo', rx_sig_cfo);
    assignin('base', 'rx_sig_synch', rx_sig_synch);
    assignin('base', 'ch_real_onetap_tf', ch_real_onetap_tf);
    assignin('base', 'ch_real_onetap_dd', ch_real_onetap_dd);
    assignin('base', 'ch_real_mat_tf', ch_real_mat_tf);
    assignin('base', 'ch_real_eff_tf', ch_real_eff_tf);
    assignin('base', 'ch_real_eff_dd', ch_real_eff_dd);
    assignin('base', 'ch_est_tf', ch_est_tf);
    assignin('base', 'ch_est_dd', ch_est_dd);
    if strcmp(test_cheq, 'dd'), assignin('base', 'ch_eff_dd', ch_eff_dd); end
end

end
