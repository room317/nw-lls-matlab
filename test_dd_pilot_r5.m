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
%   - 2020.09.15: multi-user simulation

function [qam_error, num_qam_usr, ch_est_rmse, papr_db] = test_dd_pilot_r5(snr_db, scs_khz, bw_mhz, num_slot, test_synch, test_scope, test_seed, test_ch_rmse, test_papr, test_chest, test_cheq, test_usr)

%% parameters

% set subcarrier spacing parameters
scs_khz_list = [15 30 60];
idx_scs = find(scs_khz_list == scs_khz, 1);
if isempty(idx_scs)
    error('Subcarrier spacing must be one of these: {15, 30, 60}')
end

% set bandwidth parameters
bw_mhz_list = [5 10 15 20 25 30 40 50 60 80 90 100];
idx_bw = find(bw_mhz_list == bw_mhz, 1);
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
sample_rate = scs_khz*1e3*num_fft;             % sampling rate (hz)
num_cp = round(num_fft*144/2048);                   % number of samples in cyclic prefix
num_ofdmsym_slot = 14;                              % number of ofdm symbols per slot
num_ofdmsym = num_ofdmsym_slot*num_slot;       % number of ofdm symbols
t_ofdmsym = (num_fft+num_cp)/sample_rate;           % ofdm symbol time

% set user parameters (full span, prb_size: 12*14)
if strcmp(test_usr, 'su')                   % full spreading
    num_rb_usr = num_rb;                    % number of resource blocks per user
    num_slot_usr = num_slot;                % number of slots per user
    list_usr = 1;                           % user index list
elseif strcmp(test_usr, 'mu')
    num_rb_usr = 3;                         % number of resource blocks per user
    num_slot_usr = 1;                       % number of slots per user
    list_usr = [1 8];                  % user index list
else
    error('''test_usr'' shall be either ''su'' or ''mu''.\n')
end
num_subc_usr = num_rb_usr*num_subc_rb;              % number of subcarriers per user
num_ofdmsym_usr = num_ofdmsym_slot*num_slot_usr;    % number of ofdm symbols per user (slot-based)
max_usr_rb = floor(num_rb/num_rb_usr);              % max. user rb index
max_usr_slot = floor(num_slot/num_slot_usr);        % max. user slot index
max_usr = max_usr_rb*max_usr_slot;                  % max. number of users
num_usr = length(list_usr);                         % number of users
if sum(double(list_usr > max_usr)) > 0
    error('max. user index shall be %d.\n', max_usr)
end

% set dd parameters
delay_pilot_ratio = 0.4;                % num_delay_pilot/num_delay_data (40%)
doppler_pilot_ratio = 1.0;              % num_delay_pilot/num_delay_data (100%)
delay_guard_ratio = 0.2;                % num_delay_guard/num_delay_pilot (20%)
doppler_guard_ratio = 0.2;              % num_doppler_guard/num_doppler_pilot (20%)
num_delay_usr = num_subc_usr;                       % fixed
num_doppler_usr = num_ofdmsym_usr;                  % fixed
num_delay_pilot = double(strcmp(test_chest, 'dd_tone'))*round(num_delay_usr*(delay_pilot_ratio/2))*2;       % fixed (40% of available delay grids, even number)
num_doppler_pilot = double(strcmp(test_chest, 'dd_tone'))*round(num_doppler_usr*(doppler_pilot_ratio/2))*2; % fixed (100% of available doppler grids)
num_delay_data = num_delay_usr-num_delay_pilot;
num_doppler_data = num_doppler_usr-num_doppler_pilot;
num_delay_guard = double(num_delay_usr~=num_delay_pilot)*round(num_delay_pilot*delay_guard_ratio);            % 20% of pilot delay grids
num_doppler_guard = double(num_doppler_usr~=num_doppler_pilot)*round(num_doppler_pilot*doppler_guard_ratio);    % 20% of pilot doppler grids

% set test parameter
qam_size = 16;
cfo_norm = 0;               % normalized cfo. cfo = cfo_norm/t_sym;
carrier_freq_mhz = 4000;
velocity_kmh = 500;
idx_fading = 11;            % 9: TDL-A, 10: TDL-B, 11: TDL-C, 12: TDL-D, 13: TDL-E
delay_spread_rms_us = 0.1;  % 0.1e-6;
snr_db_adj = snr_db-(10*log10((num_rb*num_slot)/(num_usr*num_rb_usr*num_slot_usr)));    % adjust time-domain snr wrt resource size
noise_var = 10^((-0.1)*snr_db)*(num_subc_bw/num_fft);                                   % calculate noise variance at dd resource block
test_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, delay_spread_rms_us);

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

%% transmitter

% calculate num. qam symbols per packet
if strcmp(test_chest, 'dd_tone')
    num_qam_usr = (num_delay_usr*num_doppler_usr)-(num_delay_pilot*num_doppler_pilot);
else
    num_qam_usr = num_delay_usr*num_doppler_usr;
end

% simulate per user
tx_bit_usr = zeros(num_qam_usr, num_usr);
tx_sym_data_usr = zeros(num_qam_usr, num_usr);
tx_sym_dd_usr = zeros(num_delay_usr, num_doppler_usr, num_usr);
tx_sym_tf_usr = zeros(num_subc_usr, num_ofdmsym_usr, num_usr);
tx_sym_tf = zeros(num_subc_bw, num_ofdmsym);
for idx_usr = 1:num_usr

    % generate bit stream
    tx_bit_usr(:, idx_usr) = randi([0 qam_size-1], num_qam_usr, 1);
    
    % modulate bit stream
    tx_sym_data_usr(:, idx_usr) = qammod(tx_bit_usr(:, idx_usr), qam_size, 'UnitAveragePower', true);
    
    % map symbols
    if strcmp(test_chest, 'dd_tone')
        % reshape data symbols
        tx_sym_data1 = reshape(tx_sym_data_usr(1:num_delay_pilot*num_doppler_data, idx_usr), num_delay_pilot, []);
        tx_sym_data2 = reshape(tx_sym_data_usr(num_delay_pilot*num_doppler_data+1:end, idx_usr), num_delay_data, []);
        
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
        tx_sym_dd_usr(:, :, idx_usr) = circshift(tx_sym_dd_base, [idx_delay_ctr-idx_delay_pilot_ctr, idx_doppler_ctr-idx_doppler_pilot_ctr]);
    else
        % reshape and map data symbols
        tx_sym_dd_usr(:, :, idx_usr) = reshape(tx_sym_data_usr(:, idx_usr), num_delay_usr, []);
    end
    
    % 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
    tx_sym_tf_usr(:, :, idx_usr) = sqrt(num_ofdmsym_usr/num_subc_usr)*fft(ifft(tx_sym_dd_usr(:, :, idx_usr), [], 2), [], 1);
    
    % map user block to whole resource blocks
    usr_id = list_usr(idx_usr);
    idx_rb_usr = (ceil(usr_id/max_usr_slot)-1)*num_rb_usr+1;
    idx_slot_usr = mod(usr_id-1, max_usr_slot)*num_slot_usr+1;
    list_subc_usr = (idx_rb_usr-1)*num_subc_rb+1:(idx_rb_usr-1)*num_subc_rb+num_subc_usr;
    list_ofdmsym_usr = (idx_slot_usr-1)*num_ofdmsym_slot+1:(idx_slot_usr-1)*num_ofdmsym_slot+num_ofdmsym_usr;
    tx_sym_tf(list_subc_usr, list_ofdmsym_usr) = tx_sym_tf_usr(:, :, idx_usr);
end

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
if test_seed >= 0
    rng(test_seed)
end
[tx_sig_faded, ch_path_gain] = fading_ch(tx_sig);       % ch_path_gain: normally constant per path, vary when doppler exists

% add gaussian noise
rx_sig = awgn(tx_sig_faded, snr_db_adj, 'measured');

% regenerate real channel
if test_scope || test_ch_rmse || strcmp(test_chest, 'real')
    [ch_real_mat_tf, ~, ~, ~] = gen_real_ch_r1(fading_ch, ch_path_gain, num_fft, num_cp, num_subc_bw, num_ofdmsym, [], [], false);
else
    ch_real_mat_tf = [];
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

% simulate per user
ch_real_onetap_usr_tf = zeros(num_subc_usr, num_ofdmsym_usr, num_usr);
ch_real_onetap_usr_dd = zeros(num_delay_usr, num_doppler_usr, num_usr);
ch_real_eff_usr_tf = zeros(num_subc_usr*num_ofdmsym_usr, num_subc_usr*num_ofdmsym_usr, num_usr);
ch_real_eff_usr_dd = zeros(num_subc_usr*num_ofdmsym_usr, num_subc_usr*num_ofdmsym_usr, num_usr);
ch_est_tf = zeros(num_subc_usr, num_ofdmsym_usr, num_usr);
ch_est_dd = zeros(num_delay_usr, num_doppler_usr, num_usr);
rx_sym_tf_usr = zeros(num_subc_usr, num_ofdmsym_usr, num_usr);
rx_sym_dd_usr = zeros(num_delay_usr, num_doppler_usr, num_usr);
rx_sym_dd_eq = zeros(num_delay_usr, num_doppler_usr, num_usr);
rx_sym_data_usr = zeros(num_qam_usr, num_usr);
rx_bit_usr = zeros(num_qam_usr, num_usr);
for idx_usr = 1:num_usr
    
    % demap user block
    usr_id = list_usr(idx_usr);
    idx_rb_usr = (ceil(usr_id/max_usr_slot)-1)*num_rb_usr+1;
    idx_slot_usr = mod(usr_id-1, max_usr_slot)*num_slot_usr+1;
    list_subc_usr = (idx_rb_usr-1)*num_subc_rb+1:(idx_rb_usr-1)*num_subc_rb+num_subc_usr;
    list_ofdmsym_usr = (idx_slot_usr-1)*num_ofdmsym_slot+1:(idx_slot_usr-1)*num_ofdmsym_slot+num_ofdmsym_usr;
    rx_sym_tf_usr(:, :, idx_usr) = rx_sym_tf(list_subc_usr, list_ofdmsym_usr);
    
    % 2d inverse sfft
    rx_sym_dd_usr(:, :, idx_usr) = sqrt(num_subc_usr/num_ofdmsym_usr)*fft(ifft(rx_sym_tf_usr(:, :, idx_usr), [], 1), [], 2);
    
    % generate real one-tap time-frequency channel matrix
    if test_scope
        % extract diagonal elements of real channel
        for idx_sym = 1:num_ofdmsym_usr
            ch_real_onetap_usr_tf(:, idx_sym, idx_usr) = diag(ch_real_mat_tf(list_subc_usr, list_subc_usr, list_ofdmsym_usr(idx_sym)));  % diagonal term only
        end
        ch_real_onetap_usr_dd(:, :, idx_usr) = sqrt(num_subc_usr/num_ofdmsym_usr)*fft(ifft(ch_real_onetap_usr_tf(:, :, idx_usr), [], 1), [], 2);
    else
        ch_real_onetap_usr_tf = [];
        ch_real_onetap_usr_dd = [];
    end
    
    % generate real effective time-frequency channel matrix
    if test_scope || test_ch_rmse || strcmp(test_chest, 'real')
        for idx_sym = 1:num_ofdmsym_usr
            ch_real_eff_usr_tf((idx_sym-1)*num_subc_usr+1:idx_sym*num_subc_usr, (idx_sym-1)*num_subc_usr+1:idx_sym*num_subc_usr, idx_usr) = ...
                ch_real_mat_tf(list_subc_usr, list_subc_usr, list_ofdmsym_usr(idx_sym));      % demap user block
        end
    else
        ch_real_eff_usr_tf = [];
    end
    
    % generate real effective delay-doppler channel matrix
    if test_scope || ((test_ch_rmse || strcmp(test_chest, 'real')) && (strcmp(test_cheq, 'dd_zf') || strcmp(test_cheq, 'dd_mmse')))
        % generate sfft matrix
        % idft_column = kron(eye(nsym), conj(dftmtx(nbw))/sqrt(nbw));
        % dft_row = kron(dftmtx(nsym)/sqrt(nsym), eye(nbw));
        % sfft_mtx = dft_row*idft_column;
        sfft_mtx = kron(dftmtx(num_ofdmsym_usr), conj(dftmtx(num_subc_usr))/num_subc_usr);      % kron(A, B)*kron(C, D) = kron(AC, BD)
        isfft_mtx = kron(conj(dftmtx(num_ofdmsym_usr))/num_ofdmsym_usr, dftmtx(num_subc_usr));  % inv(kron(A, B)) = kron(inv(A), inv(B))
        
        % generate effective delay-doppler channel matrix
        % ch_eff_dd = sfft_mtx*ch_eff_tf/sfft_mtx;
        ch_real_eff_usr_dd(:, :, idx_usr) = sfft_mtx*ch_real_eff_usr_tf(:, :, idx_usr)*isfft_mtx;
    else
        ch_real_eff_usr_dd = [];
    end
    
    % estimate channel
    if strcmp(test_chest, 'dd_tone')
        % demap pilot and remove guard
        rx_sym_dd_base = circshift(rx_sym_dd_usr(:, :, idx_usr), [-(idx_delay_ctr-idx_delay_pilot_ctr), -(idx_doppler_ctr-idx_doppler_pilot_ctr)]);
        rx_sym_pilot = rx_sym_dd_base(num_delay_guard+1:num_delay_pilot-num_delay_guard, num_doppler_guard+1:num_doppler_pilot-num_doppler_guard);
        
        % estimate channel
        ch_est_dd_base = zeros(num_delay_usr, num_doppler_usr);
        ch_est_dd_base(num_delay_guard+1:num_delay_pilot-num_delay_guard, num_doppler_guard+1:num_doppler_pilot-num_doppler_guard) = ...
            rx_sym_pilot*sqrt((num_delay_usr*num_doppler_usr)/(num_delay_pilot*num_doppler_pilot));
        ch_est_dd(:, :, idx_usr) = circshift(ch_est_dd_base, [-idx_delay_pilot_ctr+1, -idx_doppler_pilot_ctr+1]);
        
        % transform channel
        ch_est_tf(:, :, idx_usr) = sqrt(num_ofdmsym_usr/num_subc_usr)*fft(ifft(ch_est_dd(:, :, idx_usr), [], 2), [], 1);
        
    elseif strcmp(test_chest, 'perfect')
        % estimate tf channel
        ch_est_tf(:, :, idx_usr) = rx_sym_tf_usr(:, :, idx_usr)./tx_sym_tf_usr(:, :, idx_usr);     % tf domain
        
        % estimate dd channel
        ch_est_dd(:, :, idx_usr) = sqrt(num_subc_usr/num_ofdmsym_usr)*fft(ifft(ch_est_tf(:, :, idx_usr), [], 1), [], 2);
    elseif strcmp(test_chest, 'real')
        ch_est_tf = ch_real_onetap_usr_tf;
        ch_est_dd = ch_real_onetap_usr_dd;
    else
        error('test_chest value must be one of these: {''dd_tone'', ''perfect'', ''real''}')
    end
    
    % generate effective dd channel matrix (to check channel error)
    if test_scope || strcmp(test_cheq, 'dd_zf') || strcmp(test_cheq, 'dd_mmse')
        if strcmp(test_chest, 'real')
            ch_eff_dd = ch_real_eff_usr_dd(:, :, idx_usr);
        else
            ch_eff_dd = gen_eff_ch(ch_est_dd(:, :, idx_usr));
        end
    else
        ch_eff_dd = [];
    end
    
    % generate effective tf channel matrix (to check channel error)
    if test_scope || (strcmp(test_cheq, 'tf_zf') || strcmp(test_cheq, 'tf_mmse'))
        if strcmp(test_chest, 'real')
            ch_eff_tf = ch_real_eff_usr_tf(:, :, idx_usr);
        else
            ch_eff_tf = diag(reshape(ch_est_tf(:, :, idx_usr), [], 1));
        end
    else
        ch_eff_tf = [];
    end
    
    % calculate channel rmse
    if test_ch_rmse
        if strcmp(test_cheq, 'dd_zf') || strcmp(test_cheq, 'dd_mmse')
            ch_est_rmse = sqrt(mean(abs(ch_real_eff_usr_dd(:, :, idx_usr)-ch_eff_dd).^2, 'all'));
        else
            ch_est_rmse = sqrt(mean(abs(ch_real_eff_usr_tf(:, :, idx_usr)-ch_eff_tf).^2, 'all'));
        end
    else
        ch_est_rmse = [];
    end
    
    % equalize channel
    if strcmp(test_cheq, 'tf_zf')
        % equalize channel in tf domain
        if strcmp(test_chest, 'real')
            % full-tap equalization (matrix inversion)
            rx_sym_tf_eq_vec = ch_eff_tf\reshape(rx_sym_tf_usr(:, :, idx_usr), [], 1);
            rx_sym_tf_eq = reshape(rx_sym_tf_eq_vec, num_subc_usr, num_ofdmsym_usr);
        else
            % one-tap equalization (element-wise)
            rx_sym_tf_eq = rx_sym_tf_usr(:, :, idx_usr)./ch_est_tf(:, :, idx_usr);
        end
        
        % observe symbols in dd domain
        rx_sym_dd_eq(:, :, idx_usr) = sqrt(num_subc_usr/num_ofdmsym_usr)*fft(ifft(rx_sym_tf_eq, [], 1), [], 2);
    elseif strcmp(test_cheq, 'tf_mmse')
        % equalize channel in tf domain
        if strcmp(test_chest, 'real')
            % full-tap equalization (matrix inversion)
            ch_real_eff_usr_tf_mmse = ch_eff_tf'/(ch_eff_tf*ch_eff_tf'+noise_var*eye(num_subc_usr*num_ofdmsym_usr));
            rx_sym_tf_eq_vec = ch_real_eff_usr_tf_mmse*reshape(rx_sym_tf_usr(:, :, idx_usr), [], 1);
            rx_sym_tf_eq = reshape(rx_sym_tf_eq_vec, num_subc_usr, num_ofdmsym_usr);
        else
            % one-tap equalization (element-wise)
            ch_est_tf_mmse = conj(ch_est_tf(:, :, idx_usr))./(noise_var+abs(ch_est_tf(:, :, idx_usr)).^2);
            rx_sym_tf_eq = rx_sym_tf_usr(:, :, idx_usr).*ch_est_tf_mmse;
        end
        
        % observe symbols in dd domain
        rx_sym_dd_eq(:, :, idx_usr) = sqrt(num_subc_usr/num_ofdmsym_usr)*fft(ifft(rx_sym_tf_eq, [], 1), [], 2);
    elseif strcmp(test_cheq, 'dd_zf')
        % equalize channel in dd domain
        if strcmp(test_chest, 'real')
            % full-tap equalization (matrix inversion)
            rx_sym_dd_eq_vec = ch_eff_dd\reshape(rx_sym_dd_usr(:, :, idx_usr), [], 1);
            rx_sym_dd_eq(:, :, idx_usr) = reshape(rx_sym_dd_eq_vec, num_delay_usr, num_doppler_usr);
        else
            % one-tap equalization (matrix inversion)
            rx_sym_dd_eq_vec = ch_eff_dd\reshape(rx_sym_dd_usr(:, :, idx_usr), [], 1);
            rx_sym_dd_eq(:, :, idx_usr) = reshape(rx_sym_dd_eq_vec, num_delay_usr, num_doppler_usr);
        end
    elseif strcmp(test_cheq, 'dd_mmse')
        % equalize channel in dd domain
        if strcmp(test_chest, 'real')
            % full-tap equalization (matrix inversion)
            ch_real_eff_dd_mmse = ch_eff_dd'/(ch_eff_dd*ch_eff_dd'+noise_var*eye(num_delay_usr*num_doppler_usr));
            rx_sym_dd_eq_vec = ch_real_eff_dd_mmse*reshape(rx_sym_dd_usr(:, :, idx_usr), [], 1);
            rx_sym_dd_eq(:, :, idx_usr) = reshape(rx_sym_dd_eq_vec, num_delay_usr, num_doppler_usr);
        else
            % one-tap equalization (matrix inversion)
            ch_eff_dd_mmse = ch_eff_dd'/(ch_eff_dd*ch_eff_dd'+noise_var*eye(num_delay_usr*num_doppler_usr));
            rx_sym_dd_eq_vec = ch_eff_dd_mmse*reshape(rx_sym_dd_usr(:, :, idx_usr), [], 1);
            rx_sym_dd_eq(:, :, idx_usr) = reshape(rx_sym_dd_eq_vec, num_delay_usr, num_doppler_usr);
        end
    else
        error('test_chest value must be one of these: {tf_zf, tf_mmse, dd_zf, dd_mmse}')
    end
    
    % demap data symbols
    if strcmp(test_chest, 'dd_tone')
        rx_sym_dd_base = circshift(rx_sym_dd_eq(:, :, idx_usr), [-(idx_delay_ctr-idx_delay_pilot_ctr), -(idx_doppler_ctr-idx_doppler_pilot_ctr)]);
        rx_sym_data1 = rx_sym_dd_base(1:num_delay_pilot, num_doppler_pilot+1:end);
        rx_sym_data2 = rx_sym_dd_base(num_delay_pilot+1:end, :);
        rx_sym_data_usr(:, idx_usr) = [rx_sym_data1(:); rx_sym_data2(:)];
    else
        rx_sym_data_usr(:, idx_usr) = reshape(rx_sym_dd_eq(:, :, idx_usr), [], 1);
    end
    
    % demodulate qam symbols
    rx_bit_usr(:, idx_usr) = qamdemod(rx_sym_data_usr(:, idx_usr), qam_size, 'UnitAveragePower', true);
end

% calculate bit errors
if isempty(tx_bit_usr)
    qam_error = 0;
else
    qam_error = symerr(tx_bit_usr, rx_bit_usr, 'column-wise');
end

% calculate papr
if test_papr
    papr_db = 10*log10(max(abs(tx_sig).^2)/mean(abs(tx_sig).^2));
else
    papr_db = [];
end

% test scope
if test_scope
    % check fft sym
    figure
    subplot(1, 2, 1), mesh(1:num_ofdmsym, 1:num_fft, abs(tx_sym_nfft)), title('tx tf symbols (nfft)')
    subplot(1, 2, 2), mesh(1:num_ofdmsym, 1:num_fft, abs(rx_sym_nfft)), title('rx tf symbols (nfft)')
    
    % check time sig
    figure
    subplot(1, 2, 1), mesh(1:num_ofdmsym, 1:num_fft, abs(tx_ofdmsym)), title('tx t ofdm symbols')
    subplot(1, 2, 2), mesh(1:num_ofdmsym, 1:num_fft, abs(rx_ofdmsym)), title('rx t ofdm symbols')
    
    % check time sig with cp
    figure
    subplot(1, 2, 1), mesh(1:num_ofdmsym, 1:num_fft+num_cp, abs(tx_ofdmsym_cp)), title('tx t ofdm symbols with cp')
    subplot(1, 2, 2), mesh(1:num_ofdmsym, 1:num_fft+num_cp, abs(rx_ofdmsym_cp)), title('rx t ofdm symbols with cp')
    
    for idx_usr = 1:num_usr
        % check bits
        figure
        plot(tx_bit_usr(:, idx_usr), '-b.'), hold on, plot(rx_bit_usr(:, idx_usr), ':r.'), hold off, grid minor, legend('tx data bits', 'rx data bits')
        
        % check qam
        figure 
        subplot(2, 1, 1), plot(real(tx_sym_data_usr(:, idx_usr)), '-b.'), hold on, plot(real(rx_sym_data_usr(:, idx_usr)), ':r.'), hold off, grid minor, legend('tx qam symbols', 'rx qam symbols'), title('real')
        subplot(2, 1, 2), plot(imag(tx_sym_data_usr(:, idx_usr)), '-b.'), hold on, plot(imag(rx_sym_data_usr(:, idx_usr)), ':r.'), hold off, grid minor, legend('tx qam symbols', 'rx qam symbols'), title('imag')
        
        % check dd sym
        figure
        subplot(1, 3, 1), mesh(1:num_doppler_usr, 1:num_delay_usr, abs(tx_sym_dd_usr(:, :, idx_usr))), title('abs(tx dd qam sym)')
        subplot(1, 3, 2), mesh(1:num_doppler_usr, 1:num_delay_usr, abs(rx_sym_dd_usr(:, :, idx_usr))), title('abs(faded rx dd qam sym)')
        subplot(1, 3, 3), mesh(1:num_doppler_usr, 1:num_delay_usr, abs(rx_sym_dd_eq(:, :, idx_usr))), title('abs(equalized rx dd qam sym)')
        
        % check tf sym
        figure
        subplot(1, 2, 1), mesh(1:num_ofdmsym_usr, 1:num_subc_usr, abs(tx_sym_tf_usr(:, :, idx_usr))), title('abs(tx tf sym)')
        subplot(1, 2, 2), mesh(1:num_ofdmsym_usr, 1:num_subc_usr, abs(rx_sym_tf_usr(:, :, idx_usr))), title('abs(faded rx tf sym)')
        
        % check one-tap tf channel
        figure
        subplot(1, 2, 1), mesh((0:num_ofdmsym_usr-1)*t_ofdmsym*1e3, (0:num_subc_usr-1)*scs_khz*1e3*1e-6, abs(ch_real_onetap_usr_tf(:, :, idx_usr))), xlabel('ofdm symbols (ms)'), ylabel('subcarriers (khz)'), title('abs(real tf channel (one-tap))')
        subplot(1, 2, 2), mesh((0:num_ofdmsym_usr-1)*t_ofdmsym*1e3, (0:num_subc_usr-1)*scs_khz*1e3*1e-6, abs(ch_est_tf(:, :, idx_usr))), xlabel('ofdm symbols (ms)'), ylabel('subcarriers (khz)'), title('abs(estimated tf channel (one-tap))')
        
        % check one-tap dd channel
        figure
        subplot(1, 2, 1), mesh(((0:num_doppler_usr-1)-floor(num_doppler_usr/2))/(t_ofdmsym*num_doppler_usr*1e3), ((0:num_delay_usr-1)-floor(num_delay_usr/2))/(num_delay_usr*scs_khz*1e3*1e-6), abs(fftshift(fftshift(ch_real_onetap_usr_dd(:, :, idx_usr), 1), 2))), xlabel('doppler (khz)'), ylabel('delay (us)'), title('abs(real dd channel (one-tap))')
        subplot(1, 2, 2), mesh(((0:num_doppler_usr-1)-floor(num_doppler_usr/2))/(t_ofdmsym*num_doppler_usr*1e3), ((0:num_delay_usr-1)-floor(num_delay_usr/2))/(num_delay_usr*scs_khz*1e3*1e-6), abs(fftshift(fftshift(ch_est_dd(:, :, idx_usr), 1), 2))), xlabel('doppler (khz)'), ylabel('delay (us)'), title('abs(estimated dd channel (one-tap))')
        
        % check effective tf channel
        rx_sym_tf_usr_regen = ch_real_eff_usr_tf(:, :, idx_usr)*reshape(tx_sym_tf_usr(:, :, idx_usr), [], 1);
        fprintf('regenerated rx symbol rmse in tf: %10.4f\n', sqrt(mean(abs(reshape(rx_sym_tf_usr(:, :, idx_usr), [], 1)-rx_sym_tf_usr_regen).^2)))
        figure
        subplot(2, 1, 1), plot(real(reshape(rx_sym_tf_usr(:, :, idx_usr), [], 1)), '-b.'), hold on, plot(real(rx_sym_tf_usr_regen), ':r.'), hold off, legend('original rx tf signal', 'regenerated rx tf signal'), title('real')
        subplot(2, 1, 2), plot(imag(reshape(rx_sym_tf_usr(:, :, idx_usr), [], 1)), '-b.'), hold on, plot(imag(rx_sym_tf_usr_regen), ':r.'), hold off, legend('original rx tf signal', 'regenerated rx tf signal'), title('imag')
        
        % check effective dd channel
        rx_sym_dd_usr_regen = ch_real_eff_usr_dd(:, :, idx_usr)*reshape(tx_sym_dd_usr(:, :, idx_usr), [], 1);
        fprintf('regenerated rx symbol rmse in dd: %10.4f\n', sqrt(mean(abs(reshape(rx_sym_dd_usr(:, :, idx_usr), [], 1)-rx_sym_dd_usr_regen).^2)))
        figure
        subplot(2, 1, 1), plot(real(reshape(rx_sym_dd_usr(:, :, idx_usr), [], 1)), '-b.'), hold on, plot(real(rx_sym_dd_usr_regen), ':r.'), hold off, legend('original rx dd signal', 'regenerated dd tf signal'), title('real')
        subplot(2, 1, 2), plot(imag(reshape(rx_sym_dd_usr(:, :, idx_usr), [], 1)), '-b.'), hold on, plot(imag(rx_sym_dd_usr_regen), ':r.'), hold off, legend('original rx dd signal', 'regenerated dd tf signal'), title('imag')
        
        fprintf('user %d result plotted.\n', list_usr(idx_usr))
        pause
    end
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
    assignin('base', 'ch_real_mat_tf', ch_real_mat_tf);
    assignin('base', 'ch_real_eff_usr_tf', ch_real_eff_usr_tf);
    assignin('base', 'ch_real_eff_usr_dd', ch_real_eff_usr_dd);
    assignin('base', 'ch_est_tf', ch_est_tf);
    assignin('base', 'ch_est_dd', ch_est_dd);
    assignin('base', 'ch_eff_dd', ch_eff_dd);
end

end
