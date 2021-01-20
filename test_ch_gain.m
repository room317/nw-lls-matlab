% channel gain test

% test parameter
num_cnt = 1;

% parameters
num_rb = 11;                            % number of rbs
num_subc_rb = 12;                       % number of subcarriers per rb
num_subc_bw = num_rb*num_subc_rb;       % number of subcarriers in bandwidth
num_fft = 2^ceil(log2(num_subc_bw));    % fft size
num_slot = 4;                           % number of slots
num_cp = round(num_fft*144/2048);       % number of samples in cyclic prefix
num_ofdmsym_slot = 14;                  % number of ofdm symbols per slot
num_ofdmsym = num_ofdmsym_slot*num_slot;   % number of ofdm symbols

% user parameter
usr_id = 1;
max_usr_slot = 1;
idx_rb = (ceil(usr_id/max_usr_slot)-1)*num_rb+1;
idx_slot = mod(usr_id-1, max_usr_slot)*num_slot+1;
list_subc = (idx_rb-1)*num_subc_rb+1:(idx_rb-1)*num_subc_rb+num_subc_bw;
list_ofdmsym = (idx_slot-1)*num_ofdmsym_slot+1:(idx_slot-1)*num_ofdmsym_slot+num_ofdmsym;

% channel
scs_khz = 15;
sample_rate = scs_khz*1e3*num_fft;      % sampling rate (hz)
carrier_freq_mhz = 4000;
velocity_kmh = 3;
idx_fading = 11;                        % 9: TDL-A, 10: TDL-B, 11: TDL-C, 12: TDL-D, 13: TDL-E
delay_spread_rms_us = 0.1;              % 0.1e-6;
test_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, delay_spread_rms_us);

% test options
test_option.gpu_flag = false;

% test channel
ch_gain = zeros(num_subc_bw, num_cnt);
ch_gain_diag = zeros(num_subc_bw, num_cnt);
tic
for i = 1:num_cnt
    % generate channel
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
    tx_sig = randn((num_fft+num_cp)*num_ofdmsym, 1);
    [~, ch_pg] = fading_ch(tx_sig);
    [ch_full, ch_onetap_tf, ~, ch_eff_dd] = gen_real_ch_r1(fading_ch, ch_pg, num_fft, num_cp, num_subc_bw, num_ofdmsym, list_subc, list_ofdmsym, true, test_option);
    ch_onetap_dd = sqrt(num_subc_bw/num_ofdmsym)*fft(ifft(ch_onetap_tf, [], 1), [], 2);
    
    % mean channel gain
    ch_gain(:, i) = mean(sum(abs(ch_full).^2, 2), 3);
    ch_gain_diag_sym = zeros(num_subc_bw, num_ofdmsym);
    for j = 1:num_ofdmsym
        ch_gain_diag_sym(:, j) = abs(diag(ch_full(:, :, j))).^2;
    end
    ch_gain_diag(:, i) = mean(ch_gain_diag_sym, 2);
    
    % test
%     figure, plot(10*log10(mean(abs(ch_pg).^2, 1))), grid minor
%     figure(10), stem(mean(abs(ch_pg).^2, 1)), grid minor
%     pause
end
toc

% check result
ch_gain_avg = mean(ch_gain, 2);
ch_gain_diag_avg = mean(ch_gain_diag, 2);
figure, plot(1:num_subc_bw, ch_gain_avg, '-b.'), grid minor, xlabel('subcarriers'), ylabel('average channel gain'), title('Average channel gain per subcarrier')
figure, plot(1:num_subc_bw, ch_gain_diag_avg, '-b.'), grid minor, xlabel('subcarriers'), ylabel('average channel gain'), title('Average channel gain per subcarrier')

% dump
assignin('base', 'fading_ch', fading_ch)
assignin('base', 'ch_pg', ch_pg)
assignin('base', 'ch_full', ch_full)
assignin('base', 'ch_onetap_tf', ch_onetap_tf)
assignin('base', 'ch_onetap_dd', ch_onetap_dd)
assignin('base', 'ch_eff_dd', ch_eff_dd)
