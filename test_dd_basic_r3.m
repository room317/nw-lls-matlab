% test dd basic_r3 (real channel testing) is to test
%   1) real channel reproduction in all domains
%   2) system models in all domains
%
% options:
%   - scs_khz: subcarrier spacing (khz)
%   - bw_mhz: bandwidth (mhz)
%   - num_slot: number of slots
%   - test_synch: synch position for test (0 for perfect synch)
%   - test_pilot: 1) data, 2) impulse, 3) zadoff-chu sequence, 4) complementary golay sequences ('data', 'impulse', 'zc', 'golay')
%   - test_pilot_pos: vector for pilot position (2 elements for 'impulse' or 'zc', and 4 elements for 'golay')
%   - test_golay_len: complementary golay sequence length (32, 64, 128)
%   - test_scope: print and plot results
%   - test_seed: seed for random channel (negative scalar for random seed)
%   - test_partial: partial reception test (number of symbols received)
% created: 2020.03.23
% modified:
%   - 2020.03.28:z
%   - 2020.04.29:
%   - 2020.06.30: walsh_hadamard spreading (not working as intended, removed)
%   - 2020.06.30: papr calculation (removed)
%   - 2020.06.30: symbol scrambling (not working as intended, removed)
%   - 2020.07.19: latency test with partial reception (test_partial)
%   - 2020.11.24: system model check in all domains
%   - 2020.11.25: golay complementary sequence pilot test
% memo
%   - https://kr.mathworks.com/help/comm/ref/comm.rayleighchannel-system-object.html#d120e191883

function [ch_rmse_dd, ch_rmse_tf] = test_dd_basic_r3(snr_db, scs_khz, bw_mhz, num_slot, test_synch, test_pilot, test_pilot_pos, test_seq_len, test_scope, test_seed, test_partial)

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
     5   11   18   24   31   38   51   65   79  107  121  135];

% set numerical parameters
num_rb = rb_tbl(idx_scs, idx_bw);                   % number of rbs
num_subc_rb = 12;                                   % number of subcarriers per rb
num_subc_bw = num_rb*num_subc_rb;                   % number of subcarriers in bandwidth
num_fft = 2^ceil(log2(num_subc_bw));                % fft size
sample_rate = scs_khz*1e3*num_fft;             % sampling rate (hz)
num_cp = round(num_fft*144/2048);                   % number of samples in cyclic prefix
num_ofdmsym_slot = 14;                              % number of ofdm symbols per slot
num_ofdmsym = num_ofdmsym_slot*num_slot;       % number of ofdm symbols
% t_ofdmsym = (num_fft+num_cp)/sample_rate;           % ofdm symbol time

% set user parameters (full span, prb_size: 12*14)
num_rb_usr = 4; % num_rb; % num_rb; % 3;                   % number of resource blocks per user
num_slot_usr = 1; % num_slot; % num_slot; % 1;               % number of slots per user
usr_id = 1;                                 % user id
num_subc_usr = num_rb_usr*num_subc_rb;              % number of subcarriers per user
num_ofdmsym_usr = num_ofdmsym_slot*num_slot_usr;    % number of ofdm symbols per user (slot-based)
max_usr_slot = floor(num_slot/num_slot_usr);        % max. user slot index

% set dd parameters
num_delay_usr = num_subc_usr;                       % fixed
num_doppler_usr = num_ofdmsym_usr;                  % fixed

% set test parameter
qam_size = 16;
cfo_norm = 0;               % normalized cfo. cfo = cfo_norm/t_sym;
carrier_freq_mhz = 4000;
velocity_kmh = 500;
idx_fading = 11;            % 9: TDL-A, 10: TDL-B, 11: TDL-C, 12: TDL-D, 13: TDL-E
delay_spread_rms_us = 0.2;  % 0.1e-6;
test_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, delay_spread_rms_us);
test_option.gpu_flag = false;

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

% % check error
% if test_synch > len_cp
%     error('test_synch should not be greater than %d.\n', len_cp)
% end

%% transmitter

if strcmp(test_pilot, 'data')
    % calculate num. qam symbols per packet
    num_qam_usr = num_delay_usr*num_doppler_usr;
    
    % generate bit stream
    tx_bit_usr = randi([0 qam_size-1], num_qam_usr, 1);
    
    % modulate bit stream
    tx_sym_data_usr = qammod(tx_bit_usr, qam_size, 'UnitAveragePower', true);
    
    % reshape and map data symbols
    tx_sym_usr_dd = reshape(tx_sym_data_usr, num_delay_usr, []);
elseif strcmp(test_pilot, 'impulse')
    % check position
    if ~(iscell(test_pilot_pos) && length(test_pilot_pos) == 1 && length(test_pilot_pos{1}) == 2)
        error('''test_pilot_pos'' shall be a cell with a vector which has 2 elements.')
    end
    
    % map delay-doppler impulse
    tx_sym_usr_dd = zeros(num_delay_usr, num_doppler_usr);
    tx_sym_usr_dd(test_pilot_pos{1}(1)+1, test_pilot_pos{1}(2)+1) = sqrt(num_delay_usr*num_doppler_usr);
elseif strcmp(test_pilot, 'zc')
    % check position
    if ~(iscell(test_pilot_pos) && length(test_pilot_pos) == 1 && length(test_pilot_pos{1}) == 2)
        error('''test_pilot_pos'' shall be a cell with a vector which has 2 elements.')
    end
    
    % generate complementary golay sequences
    if test_seq_len > num_delay_usr
        error('''test_seq_len'' shall be equal or smaller than %d.', num_delay_usr)
    else
        zc_seq = zadoffChuSeq(1, test_seq_len);
    end
    
    % map sequences
    tx_sym_usr_dd = zeros(num_delay_usr, num_doppler_usr);
    tx_sym_usr_dd(test_pilot_pos{1}(1)+1:test_pilot_pos{1}(1)+test_seq_len, test_pilot_pos{1}(2)+1) = zc_seq*sqrt(num_delay_usr*num_doppler_usr/test_seq_len);
elseif strcmp(test_pilot, 'golay')
    % check position
    if ~(iscell(test_pilot_pos) && length(test_pilot_pos) == 2 && length(test_pilot_pos{1}) == 2 && length(test_pilot_pos{2}) == 2)
        error('''test_pilot_pos'' shall be a cell with 2 vectors, each of which has 2 elements.')
    elseif (test_pilot_pos{1}(1) < 0) || (test_pilot_pos{2}(1) < 0)
        error('''test_pilot_pos'' delay index shall be equal or greater than 0.')
    elseif (test_pilot_pos{1}(1)+test_seq_len > num_delay_usr) || (test_pilot_pos{2}(1) > num_delay_usr)
        error('''test_pilot_pos'' delay index is too high for the sequences to fit in the resource.')
    elseif (test_pilot_pos{1}(2) < 0) || (test_pilot_pos{2}(2) < 0)
        error('''test_pilot_pos'' doppler index shall be equal or greater than 0.')
    elseif ~(test_pilot_pos{1}(2) < num_doppler_usr) || ~(test_pilot_pos{2}(2) < num_doppler_usr)
        error('''test_pilot_pos'' doppler index is too high for the sequences to fit in the resource.')
    end
    
    % generate complementary golay sequences
    if ~ismember(test_seq_len, [32, 64 128])
        error('''test_seq_len'' shall be one of thes: {32, 64, 128}.')
    elseif test_seq_len > num_delay_usr
        error('''test_seq_len'' shall be equal or smaller than %d.', num_delay_usr)
    else
        [Ga, Gb] = wlanGolaySequence(test_seq_len);
    end
    
    % map sequences
    tx_sym_usr_dd = zeros(num_delay_usr, num_doppler_usr);
    tx_sym_usr_dd(test_pilot_pos{1}(1)+1:test_pilot_pos{1}(1)+test_seq_len, test_pilot_pos{1}(2)+1) = Ga*sqrt(num_delay_usr*num_doppler_usr/test_seq_len/2);
    tx_sym_usr_dd(test_pilot_pos{2}(1)+1:test_pilot_pos{2}(1)+test_seq_len, test_pilot_pos{2}(2)+1) = Gb*sqrt(num_delay_usr*num_doppler_usr/test_seq_len/2);
else
    error('test_pilot shall be one of these: {''data'', ''impulse'', ''golay''}.')
end

% 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
tx_sym_usr_tf = sqrt(num_ofdmsym_usr/num_subc_usr)*fft(ifft(tx_sym_usr_dd, [], 2), [], 1);

% map user block to whole resource blocks
idx_rb_usr = (ceil(usr_id/max_usr_slot)-1)*num_rb_usr+1;
idx_slot_usr = mod(usr_id-1, max_usr_slot)*num_slot_usr+1;
list_subc_usr = (idx_rb_usr-1)*num_subc_rb+1:(idx_rb_usr-1)*num_subc_rb+num_subc_usr;
list_ofdmsym_usr = (idx_slot_usr-1)*num_ofdmsym_slot+1:(idx_slot_usr-1)*num_ofdmsym_slot+num_ofdmsym_usr;
tx_sym_tf = zeros(num_subc_bw, num_ofdmsym);
tx_sym_tf(list_subc_usr, list_ofdmsym_usr) = tx_sym_usr_tf;

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
if isempty(snr_db)
    rx_sig = tx_sig_faded;
%     rx_sig = tx_sig;
else
    rx_sig = awgn(tx_sig_faded, snr_db, 'measured');
end

%% receiver

% compensate cfo (to observe impact of frequency shift)
cfo_sample = 0:length(rx_sig)-1;
cfo_vec = cfo_sample.'*cfo_norm/(num_fft+num_cp);
rx_sig_cfo = rx_sig(:) .* exp(-1i*2*pi*cfo_vec);

% synchronization
rx_sig_synch = circshift(rx_sig_cfo, -num_cp-test_synch);

% reshape
rx_ofdmsym_cp = reshape(rx_sig_synch, num_fft+num_cp, []);

% remove cp (with synch: to observe impact of time shift)
rx_ofdmsym = rx_ofdmsym_cp(1:num_fft, :);

% remove parts of received symbols
if isscalar(test_partial) && ismember(test_partial, 1:num_ofdmsym-1)
    rx_ofdmsym_partial = [rx_ofdmsym(:, 1:test_partial) zeros(size(rx_ofdmsym, 1), num_ofdmsym-test_partial)];
else
    rx_ofdmsym_partial = rx_ofdmsym;
end

% ofdm demodulate the symbol
rx_sym_nfft = (1/sqrt(num_fft))*fft(rx_ofdmsym_partial, [], 1);

% demap resource plane
rx_sym_nfft_shift = fftshift(rx_sym_nfft, 1);
rx_sym_tf = rx_sym_nfft_shift((num_fft/2)-(num_subc_bw/2)+1:(num_fft/2)+(num_subc_bw/2), :);

% demap user block
rx_sym_usr_tf = rx_sym_tf(list_subc_usr, list_ofdmsym_usr);

% 2d inverse sfft
rx_sym_usr_dd = sqrt(num_subc_usr/num_ofdmsym_usr)*fft(ifft(rx_sym_usr_tf, [], 1), [], 2);

%% channel regeneration

% regenerate real channel
[ch_mat_t, ~, ch_real_onetap_usr_tf, ch_real_eff_usr_tf, ch_real_eff_usr_dd] = ...
    gen_real_ch_r1(fading_ch, ch_path_gain, num_fft, num_cp, num_subc_bw, num_ofdmsym, list_subc_usr, list_ofdmsym_usr, true, test_option);

% regenerate one-tap delay-doppler channel
if strcmp(test_pilot, 'data')
    ch_est_onetap_usr_tf = ch_real_onetap_usr_tf;
    ch_est_onetap_usr_dd = sqrt(num_subc_usr/num_ofdmsym_usr)*fft(ifft(ch_real_onetap_usr_tf, [], 1), [], 2);
elseif strcmp(test_pilot, 'impulse')
    ch_est_onetap_usr_dd = circshift(rx_sym_usr_dd, (-1)*test_pilot_pos{1});
    ch_est_onetap_usr_tf = sqrt(num_ofdmsym_usr/num_subc_usr)*fft(ifft(ch_est_onetap_usr_dd, [], 2), [], 1);
    
    if test_scope
        % rx signal
        rx_dd = circshift(rx_sym_usr_dd, (-1)*test_pilot_pos{1});
        rx_tf = rx_sym_usr_tf;
        
        % real channel
        ch_dd = sqrt(num_subc_usr/num_ofdmsym_usr)*fft(ifft(ch_real_onetap_usr_tf, [], 1), [], 2);
        ch_tf = ch_real_onetap_usr_tf;
        
        figure
        subplot(1, 2, 1), mesh(1:size(ch_dd, 2), 1:size(ch_dd, 1), abs(ch_dd).^2), title('dd channel')
        subplot(1, 2, 2), mesh(1:size(rx_dd, 2), 1:size(rx_dd, 1), abs(rx_dd).^2), title('dd rx signal')
        
        figure
        subplot(1, 2, 1), mesh(1:size(ch_tf, 2), 1:size(ch_tf, 1), abs(ch_tf).^2), title('tf channel')
        subplot(1, 2, 2), mesh(1:size(rx_tf, 2), 1:size(rx_tf, 1), abs(rx_tf).^2), title('tf rx signal')
        pause
    end
elseif strcmp(test_pilot, 'zc')
    % circular shift resources
    rx_sym_usr_cir_dd = circshift(rx_sym_usr_dd, (-1)*test_pilot_pos{1});
    
    % xcorr sequences
    ch_est_usr_seq_dd = zeros(2*num_delay_usr-1, num_doppler_usr);
    for idx_sym = 1:num_doppler_usr
        ch_est_usr_seq_dd(:, idx_sym) = xcorr(rx_sym_usr_cir_dd(:, idx_sym), zc_seq/sqrt(test_seq_len));
    end
    
    % extract valid resource
    ch_est_onetap_usr_valid_dd = ch_est_usr_seq_dd(ceil(num_delay_usr/2):ceil(num_delay_usr/2)+num_delay_usr-1, :);
    
    % shift pilot to origin
    ch_est_onetap_usr_dd = circshift(ch_est_onetap_usr_valid_dd, (-1)*[ceil(num_delay_usr/2), 0]);
    
    % 2d sfft channel
    ch_est_onetap_usr_tf = sqrt(num_ofdmsym_usr/num_subc_usr)*fft(ifft(ch_est_onetap_usr_dd, [], 2), [], 1);
    
    if test_scope
        % tx signal in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(tx_sym_usr_dd, 2), 1:size(tx_sym_usr_dd, 1), real(tx_sym_usr_dd)), title('tx real')
        subplot(1, 2, 2), mesh(1:size(tx_sym_usr_dd, 2), 1:size(tx_sym_usr_dd, 1), imag(tx_sym_usr_dd)), title('tx imag')
        
        % rx signal in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(rx_sym_usr_dd, 2), 1:size(rx_sym_usr_dd, 1), real(rx_sym_usr_dd)), title('rx real')
        subplot(1, 2, 2), mesh(1:size(rx_sym_usr_dd, 2), 1:size(rx_sym_usr_dd, 1), imag(rx_sym_usr_dd)), title('rx imag')
        
        % shifted rx sequence in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(rx_sym_usr_cir_dd, 2), 1:size(rx_sym_usr_cir_dd, 1), real(rx_sym_usr_cir_dd)), title('xcorr real')
        subplot(1, 2, 2), mesh(1:size(rx_sym_usr_cir_dd, 2), 1:size(rx_sym_usr_cir_dd, 1), imag(rx_sym_usr_cir_dd)), title('xcorr imag')
        
        % rx sequence after xcorr in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(ch_est_usr_seq_dd, 2), 1:size(ch_est_usr_seq_dd, 1), real(ch_est_usr_seq_dd)), title('xcorr real')
        subplot(1, 2, 2), mesh(1:size(ch_est_usr_seq_dd, 2), 1:size(ch_est_usr_seq_dd, 1), imag(ch_est_usr_seq_dd)), title('xcorr imag')
        
        % channel estimation after circular shift in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(ch_est_onetap_usr_dd, 2), 1:size(ch_est_onetap_usr_dd, 1), real(ch_est_onetap_usr_dd)), title('real')
        subplot(1, 2, 2), mesh(1:size(ch_est_onetap_usr_dd, 2), 1:size(ch_est_onetap_usr_dd, 1), imag(ch_est_onetap_usr_dd)), title('imag'), pause
        
        % one-tap channel estimation vs real one-tap channel for comparison in dd domain
        ch_real_onetap_usr_dd = sqrt(num_subc_usr/num_ofdmsym_usr)*fft(ifft(ch_real_onetap_usr_tf, [], 1), [], 2);
        figure
        subplot(1, 2, 1), mesh(1:size(ch_real_onetap_usr_dd, 2), 1:size(ch_real_onetap_usr_dd, 1), abs(fftshift(fftshift(ch_real_onetap_usr_dd, 1), 2)).^2), title('real channel'), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power')
        subplot(1, 2, 2), mesh(1:size(ch_est_onetap_usr_dd, 2), 1:size(ch_est_onetap_usr_dd, 1), abs(fftshift(fftshift(ch_est_onetap_usr_dd, 1), 2)).^2), title('channel estimation'), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power'), pause
        
%         % test plot for documents
%         figure
%         subplot(1, 3, 1), mesh(1:size(ch_est_usr_seq1_dd, 2), 1:size(ch_est_usr_seq1_dd, 1), abs(ch_est_usr_seq1_dd).^2), title({'Autocorrelation of'; 'Ga Area'}), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power')
%         subplot(1, 3, 2), mesh(1:size(ch_est_usr_seq2_dd, 2), 1:size(ch_est_usr_seq2_dd, 1), abs(ch_est_usr_seq2_dd).^2), title({'Autocorrelation of'; 'Gb Area'}), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power')
%         subplot(1, 3, 3), mesh(1:size(ch_est_usr_seq_dd, 2), 1:size(ch_est_usr_seq_dd, 1), abs(ch_est_usr_seq_dd).^2), title({'Sum of'; 'two autocorrelations'}), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power'), pause
    end
    
elseif strcmp(test_pilot, 'golay')
    % set index (make seq1 to be located at lower left)
    [~, seq1] = min([test_pilot_pos{1}(2), test_pilot_pos{2}(2)]);
    if seq1 ~= 1
        seq2 = 1;
        golay_seq1 = Gb/sqrt(test_seq_len*2);
        golay_seq2 = Ga/sqrt(test_seq_len*2);
    else
        seq2 = 2;
        golay_seq1 = Ga/sqrt(test_seq_len*2);
        golay_seq2 = Gb/sqrt(test_seq_len*2);
    end
    
    % calculate valid sequence resources
    num_rsc_left = floor(mod(test_pilot_pos{seq1}(2)-test_pilot_pos{seq2}(2), num_ofdmsym_usr)/2);
    num_rsc_right = floor(mod(test_pilot_pos{seq2}(2)-test_pilot_pos{seq1}(2), num_ofdmsym_usr)/2);
    num_rsc_down = min(test_pilot_pos{seq1}(1), test_pilot_pos{seq2}(1));
    num_rsc_up = min(num_subc_usr-test_pilot_pos{seq1}(1), num_subc_usr-test_pilot_pos{seq2}(1));
    
    % set pilot size
    pilot_shift = [num_rsc_down, min(num_rsc_left, num_rsc_right)];
    pilot_size = [num_rsc_down+num_rsc_up, 2*pilot_shift(2)+1];
    
    % extract valid sequence resources
    rx_sym_usr_circ1_dd = circshift(rx_sym_usr_dd, (-1)*test_pilot_pos{seq1}+pilot_shift);
    rx_sym_usr_circ2_dd = circshift(rx_sym_usr_dd, (-1)*test_pilot_pos{seq2}+pilot_shift);
    rx_sym_usr_seq1_dd = rx_sym_usr_circ1_dd(1:pilot_size(1), 1:pilot_size(2));
    rx_sym_usr_seq2_dd = rx_sym_usr_circ2_dd(1:pilot_size(1), 1:pilot_size(2));
    
    % xcorr sequences
    pilot_xcorr_size = [2*max(pilot_size(1), test_seq_len)-1, pilot_size(2)];
    ch_est_usr_seq1_dd = zeros(pilot_xcorr_size);
    ch_est_usr_seq2_dd = zeros(pilot_xcorr_size);
    for idx_sym = 1:pilot_xcorr_size(2)
        ch_est_usr_seq1_dd(:, idx_sym) = xcorr(rx_sym_usr_seq1_dd(:, idx_sym), golay_seq1);
        ch_est_usr_seq2_dd(:, idx_sym) = xcorr(rx_sym_usr_seq2_dd(:, idx_sym), golay_seq2);
    end
    ch_est_usr_seq_dd = ch_est_usr_seq1_dd+ch_est_usr_seq2_dd;
    
    % demap valid resource area
    if size(ch_est_usr_seq_dd, 1) > num_subc_usr
        impulse_subc_pos = abs(pilot_size(1)-test_seq_len)+num_rsc_down+test_seq_len-1;
        impulse_subc_pos_new = min(floor(num_subc_usr/2), impulse_subc_pos);
        impulse_subc_shift = impulse_subc_pos-impulse_subc_pos_new;
        ch_est_usr_rcs_dd = ch_est_usr_seq_dd(impulse_subc_shift+1:impulse_subc_shift+num_subc_usr, :);
        impulse_shift = [impulse_subc_pos_new, pilot_shift(2)];
    else
        ch_est_usr_rcs_dd = ch_est_usr_seq_dd;
        impulse_shift = [abs(pilot_size(1)-test_seq_len)+num_rsc_down+test_seq_len-1, pilot_shift(2)];
    end
    
    % reallocate pilots and shift pilot to origin
    ch_est_onetap_usr_circ_dd = zeros(num_subc_usr, num_ofdmsym_usr);
    ch_est_onetap_usr_circ_dd(1:size(ch_est_usr_rcs_dd, 1), 1:size(ch_est_usr_rcs_dd, 2)) = ch_est_usr_rcs_dd;
    ch_est_onetap_usr_dd = circshift(ch_est_onetap_usr_circ_dd, (-1)*impulse_shift);
    
    % 2d sfft channel
    ch_est_onetap_usr_tf = sqrt(num_ofdmsym_usr/num_subc_usr)*fft(ifft(ch_est_onetap_usr_dd, [], 2), [], 1);
    
    if test_scope
        % tx signal in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(tx_sym_usr_dd, 2), 1:size(tx_sym_usr_dd, 1), real(tx_sym_usr_dd)), title('real')
        subplot(1, 2, 2), mesh(1:size(tx_sym_usr_dd, 2), 1:size(tx_sym_usr_dd, 1), imag(tx_sym_usr_dd)), title('imag')
        
        % rx signal in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(rx_sym_usr_dd, 2), 1:size(rx_sym_usr_dd, 1), real(rx_sym_usr_dd)), title('real')
        subplot(1, 2, 2), mesh(1:size(rx_sym_usr_dd, 2), 1:size(rx_sym_usr_dd, 1), imag(rx_sym_usr_dd)), title('imag')
        
        % two rx sequences in dd domain
        figure
        subplot(2, 2, 1), mesh(1:size(rx_sym_usr_seq1_dd, 2), 1:size(rx_sym_usr_seq1_dd, 1), real(rx_sym_usr_seq1_dd)), title('seq1 real')
        subplot(2, 2, 2), mesh(1:size(rx_sym_usr_seq2_dd, 2), 1:size(rx_sym_usr_seq2_dd, 1), real(rx_sym_usr_seq2_dd)), title('seq2 real')
        subplot(2, 2, 3), mesh(1:size(rx_sym_usr_seq1_dd, 2), 1:size(rx_sym_usr_seq1_dd, 1), imag(rx_sym_usr_seq1_dd)), title('seq1 imag')
        subplot(2, 2, 4), mesh(1:size(rx_sym_usr_seq2_dd, 2), 1:size(rx_sym_usr_seq2_dd, 1), imag(rx_sym_usr_seq2_dd)), title('seq2 imag')
        
        % two rx sequences after xcorr in dd domain
        figure
        subplot(2, 2, 1), mesh(1:size(ch_est_usr_seq1_dd, 2), 1:size(ch_est_usr_seq1_dd, 1), real(ch_est_usr_seq1_dd)), title('seq1 real')
        subplot(2, 2, 2), mesh(1:size(ch_est_usr_seq2_dd, 2), 1:size(ch_est_usr_seq2_dd, 1), real(ch_est_usr_seq2_dd)), title('seq2 real')
        subplot(2, 2, 3), mesh(1:size(ch_est_usr_seq1_dd, 2), 1:size(ch_est_usr_seq1_dd, 1), imag(ch_est_usr_seq1_dd)), title('seq1 imag')
        subplot(2, 2, 4), mesh(1:size(ch_est_usr_seq2_dd, 2), 1:size(ch_est_usr_seq2_dd, 1), imag(ch_est_usr_seq2_dd)), title('seq2 imag')
        
        % sum of two rx sequences after xcorr in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(ch_est_usr_seq_dd, 2), 1:size(ch_est_usr_seq_dd, 1), real(ch_est_usr_seq_dd)), title('real')
        subplot(1, 2, 2), mesh(1:size(ch_est_usr_seq_dd, 2), 1:size(ch_est_usr_seq_dd, 1), imag(ch_est_usr_seq_dd)), title('imag')
        
        % demap valid channel impulse response area in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(ch_est_usr_rcs_dd, 2), 1:size(ch_est_usr_rcs_dd, 1), real(ch_est_usr_rcs_dd)), title('real')
        subplot(1, 2, 2), mesh(1:size(ch_est_usr_rcs_dd, 2), 1:size(ch_est_usr_rcs_dd, 1), imag(ch_est_usr_rcs_dd)), title('imag')
        
        % channel estimation before circular shift in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(ch_est_onetap_usr_circ_dd, 2), 1:size(ch_est_onetap_usr_circ_dd, 1), real(ch_est_onetap_usr_circ_dd)), title('real')
        subplot(1, 2, 2), mesh(1:size(ch_est_onetap_usr_circ_dd, 2), 1:size(ch_est_onetap_usr_circ_dd, 1), imag(ch_est_onetap_usr_circ_dd)), title('imag')
        
        % channel estimation after circular shift in dd domain
        figure
        subplot(1, 2, 1), mesh(1:size(ch_est_onetap_usr_dd, 2), 1:size(ch_est_onetap_usr_dd, 1), real(ch_est_onetap_usr_dd)), title('real')
        subplot(1, 2, 2), mesh(1:size(ch_est_onetap_usr_dd, 2), 1:size(ch_est_onetap_usr_dd, 1), imag(ch_est_onetap_usr_dd)), title('imag'), pause
        
        % one-tap channel estimation vs real one-tap channel for comparison in dd domain
        ch_real_onetap_usr_dd = sqrt(num_subc_usr/num_ofdmsym_usr)*fft(ifft(ch_real_onetap_usr_tf, [], 1), [], 2);
        figure
        subplot(1, 2, 1), mesh(1:size(ch_real_onetap_usr_dd, 2), 1:size(ch_real_onetap_usr_dd, 1), abs(fftshift(fftshift(ch_real_onetap_usr_dd, 1), 2)).^2, 'EdgeColor', 'black'), title('Real channel'), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power')
        subplot(1, 2, 2), mesh(1:size(ch_est_onetap_usr_dd, 2), 1:size(ch_est_onetap_usr_dd, 1), abs(fftshift(fftshift(ch_est_onetap_usr_dd, 1), 2)).^2, 'EdgeColor', 'black'), title('Channel estimation'), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power'), pause
        
        assignin('base', 'ch_real_onetap_usr_dd', ch_real_onetap_usr_dd)
        assignin('base', 'ch_est_onetap_usr_dd', ch_est_onetap_usr_dd)
        assignin('base', 'ch_est_usr_seq1_dd', ch_est_usr_seq1_dd)
        assignin('base', 'ch_est_usr_seq2_dd', ch_est_usr_seq2_dd)
        assignin('base', 'ch_est_usr_seq_dd', ch_est_usr_seq_dd)
        
        
        % test plot for documents
        figure
%         subplot(1, 3, 1), mesh(1:size(ch_est_usr_seq1_dd, 2), 1:size(ch_est_usr_seq1_dd, 1), abs(ch_est_usr_seq1_dd).^2), title({'Autocorrelation of'; 'Ga Area'}), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power')
%         subplot(1, 3, 2), mesh(1:size(ch_est_usr_seq2_dd, 2), 1:size(ch_est_usr_seq2_dd, 1), abs(ch_est_usr_seq2_dd).^2), title({'Autocorrelation of'; 'Gb Area'}), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power')
%         subplot(1, 3, 3), mesh(1:size(ch_est_usr_seq_dd, 2), 1:size(ch_est_usr_seq_dd, 1), abs(ch_est_usr_seq_dd).^2), title({'Sum of'; 'two autocorrelations'}), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power'), pause
        subplot(2, 2, 1), mesh(1:size(ch_est_usr_seq1_dd, 2), 1:size(ch_est_usr_seq1_dd, 1), abs(ch_est_usr_seq1_dd).^2, 'EdgeColor', 'black'), title({'Autocorrelation of'; '\itG_a Area'}), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power')
        subplot(2, 2, 2), mesh(1:size(ch_est_usr_seq2_dd, 2), 1:size(ch_est_usr_seq2_dd, 1), abs(ch_est_usr_seq2_dd).^2, 'EdgeColor', 'black'), title({'Autocorrelation of'; '\itG_b Area'}), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power')
        subplot('Position', [0.35 0.1 0.35 0.35]), mesh(1:size(ch_est_usr_seq_dd, 2), 1:size(ch_est_usr_seq_dd, 1), abs(ch_est_usr_seq_dd).^2, 'EdgeColor', 'black'), title('Sum of two autocorrelations'), xlabel('Doppler'), ylabel('Delay'), zlabel('channel power'), pause
    end
end

%% output

% calculate delay-doppler channel error (full vs. one-tap)
ch_rmse_dd = sqrt(mean(abs(ch_real_eff_usr_dd-gen_eff_ch(ch_est_onetap_usr_dd)).^2, 'all'));

% calculate time-freq channel error (full vs. one-tap)
ch_rmse_tf = sqrt(mean(abs(ch_real_eff_usr_tf-diag(ch_est_onetap_usr_tf(:))).^2, 'all'));

%% regenerate rx signals

% regenerate time domain received signal (time variant channel)
if test_scope
    ch_info = info(fading_ch);
    ch_filter_coeff = ch_info.ChannelFilterCoefficients;    % ch_filter_coeff: channel filter coefficient per path
    num_path = length(fading_ch.PathDelays);                % num_path: number of path
    ch_delayed = zeros(size(tx_sig, 1), num_path);
    for idx_path = 1:num_path
        ch_delayed(:, idx_path) = filter(ch_filter_coeff(idx_path, :), 1, tx_sig);
    end
    tx_sig_faded_regen_tv = sum(ch_path_gain.*ch_delayed, 2);
    
    figure
    subplot(2, 1, 1), plot(real(tx_sig_faded), '-b.'), hold on, plot(real(tx_sig_faded_regen_tv), ':r.'), hold off, grid minor, title('time domaim')
    subplot(2, 1, 2), plot(imag(tx_sig_faded), '-b.'), hold on, plot(imag(tx_sig_faded_regen_tv), ':r.'), hold off, grid minor
    pause
end

% regenerate time domain received signal (time invariant channel: filter/conv operation)
if test_scope
    rx_ofdmsym_regen_ti = zeros(size(rx_ofdmsym));
    for idx_sym = 1:size(tx_ofdmsym, 2)
        rx_ofdmsym_regen_ti(:, idx_sym) = filter(ch_mat_t(:, 1, idx_sym), 1, tx_ofdmsym(:, idx_sym));
    end
    
    figure
    subplot(2, 1, 1), plot(real(rx_ofdmsym(:)), '-b.'), hold on, plot(real(rx_ofdmsym_regen_ti(:)), ':r.'), hold off, grid minor, title('time domaim: block fading assumed')
    subplot(2, 1, 2), plot(imag(rx_ofdmsym(:)), '-b.'), hold on, plot(imag(rx_ofdmsym_regen_ti(:)), ':r.'), hold off, grid minor
    fprintf('rx signal rmse in t domain: %6.4f\n', sqrt(mean(abs(rx_ofdmsym(:)-rx_ofdmsym_regen_ti(:)).^2, 'all')))
    pause
end

% regenerate time-frequency domain rx signals (time variant channel)
if test_scope
    rx_sym_usr_regen_tv_tf = ch_real_eff_usr_tf*tx_sym_usr_tf(:);
    
    figure
    subplot(2, 1, 1), plot(real(rx_sym_usr_tf(:)), '-b.'), hold on, plot(real(rx_sym_usr_regen_tv_tf), ':r.'), hold off, grid minor, title('time-freq domaim')
    subplot(2, 1, 2), plot(imag(rx_sym_usr_tf(:)), '-b.'), hold on, plot(imag(rx_sym_usr_regen_tv_tf), ':r.'), hold off, grid minor
    pause
end

% regenerate time-frequency domain rx signals (time invariant channel: hadamard product)
if test_scope
    rx_sym_usr_regen_ti_tf = ch_est_onetap_usr_tf.*tx_sym_usr_tf;
    
    figure
    subplot(2, 1, 1), plot(real(rx_sym_usr_tf(:)), '-b.'), hold on, plot(real(rx_sym_usr_regen_ti_tf(:)), ':r.'), hold off, grid minor, title('time-freq domaim: block fading assumed')
    subplot(2, 1, 2), plot(imag(rx_sym_usr_tf(:)), '-b.'), hold on, plot(imag(rx_sym_usr_regen_ti_tf(:)), ':r.'), hold off, grid minor
    fprintf('rx signal rmse in tf domain: %6.4f\n', sqrt(mean(abs(rx_sym_usr_tf(:)-rx_sym_usr_regen_ti_tf(:)).^2, 'all')))
    pause
end

% regenerate delay-doppler domain rx signals (time variant channel)
if test_scope
    rx_sym_usr_regen_tv_dd = ch_real_eff_usr_dd*tx_sym_usr_dd(:);
    
    figure
    subplot(2, 1, 1), plot(real(rx_sym_usr_dd(:)), '-b.'), hold on, plot(real(rx_sym_usr_regen_tv_dd), ':r.'), hold off, grid minor, title('delay-doppler domaim')
    subplot(2, 1, 2), plot(imag(rx_sym_usr_dd(:)), '-b.'), hold on, plot(imag(rx_sym_usr_regen_tv_dd), ':r.'), hold off, grid minor
    pause
end

% regenerate delay-doppler domain rx signals (time invariant channel: 2d filter/conv operation)
if test_scope
    rx_sym_usr_regen_ti_cir_dd = conv2(repmat(ch_est_onetap_usr_dd, 2, 2), tx_sym_usr_dd);
    rx_sym_usr_regen_ti_dd = rx_sym_usr_regen_ti_cir_dd(num_subc_usr+1:2*num_subc_usr, num_ofdmsym_usr+1:2*num_ofdmsym_usr)/sqrt(num_subc_usr*num_ofdmsym_usr);
    
    figure
    subplot(2, 1, 1), plot(real(rx_sym_usr_dd(:)), '-b.'), hold on, plot(real(rx_sym_usr_regen_ti_dd(:)), ':r.'), hold off, grid minor, title('delay-doppler domaim: block fading assumed')
    subplot(2, 1, 2), plot(imag(rx_sym_usr_dd(:)), '-b.'), hold on, plot(imag(rx_sym_usr_regen_ti_dd(:)), ':r.'), hold off, grid minor
    fprintf('rx signal rmse in dd domain: %6.4f\n', sqrt(mean(abs(rx_sym_usr_dd(:)-rx_sym_usr_regen_ti_dd(:)).^2, 'all')))
    pause
end

%% dump

if test_scope
    assignin('base', 'tx_sym_usr_dd', tx_sym_usr_dd);
    assignin('base', 'tx_sym_usr_tf', tx_sym_usr_tf);
    assignin('base', 'tx_sym_tf', tx_sym_tf);
    assignin('base', 'tx_sym_nfft', tx_sym_nfft);
    assignin('base', 'tx_ofdmsym', tx_ofdmsym);
    assignin('base', 'tx_ofdmsym_cp', tx_ofdmsym_cp);
    assignin('base', 'tx_sig', tx_sig);
    assignin('base', 'rx_sig', rx_sig);
    assignin('base', 'rx_ofdmsym_cp', rx_ofdmsym_cp);
    assignin('base', 'rx_ofdmsym', rx_ofdmsym);
    assignin('base', 'rx_sym_nfft', rx_sym_nfft);
    assignin('base', 'rx_sym_tf', rx_sym_tf);
    assignin('base', 'rx_sym_usr_tf', rx_sym_usr_tf);
    assignin('base', 'rx_sym_usr_dd', rx_sym_usr_dd);
end

end
