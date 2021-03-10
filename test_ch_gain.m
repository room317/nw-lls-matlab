% channel gain test

% simulation
num_sim = 100;

% channel
carrier_freq_mhz = 4000;
velocity_kmh = 500;
delay_spread_rms_us = 0.01; % us
idx_fading = 13;                        % 9: TDL-A, 10: TDL-B, 11: TDL-C, 12: TDL-D, 13: TDL-E

% numerology
bw_mhz = 10;
scs_khz = 15;
num_slot = 4;
len_tb_bit = 1748;
mcs = 8;
chest_option = 'real';
waveform = 'ofdm';
usr_option = 'mu';

% set test simulation option
sim_option.override = true;     % override simulation options when true
sim_option.num_rb = 11;         % number of total rbs
sim_option.scs_khz = 60;        % subcarrier spacing (khz)
sim_option.num_slot = 4;        % number of total slots
sim_option.Qm = 4;              % number of bits per qam symbol
sim_option.coderate = 666;      % number of information bits per 1024 bits

% generate parameter
num = nw_num_prm_r1(scs_khz, bw_mhz, num_slot, waveform, chest_option, usr_option, sim_option);
ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, delay_spread_rms_us);

% user parameter
idx_usr = 1;
usr_id = num.list_usr(idx_usr);                                         % user id
idx_rb_usr = (ceil(usr_id/num.max_usr_slot)-1)*num.num_rb_usr+1;        % user rb index
idx_slot_usr = mod(usr_id-1, num.max_usr_slot)*num.num_slot_usr+1;      % user slot index
list_subc_usr = (idx_rb_usr-1)*num.num_subc_rb+1:(idx_rb_usr-1)*num.num_subc_rb+num.num_subc_usr;                       % user subcarrier list
list_ofdmsym_usr = (idx_slot_usr-1)*num.num_ofdmsym_slot+1:(idx_slot_usr-1)*num.num_ofdmsym_slot+num.num_ofdmsym_usr;   % user ofdm symbol list

% test options
test_option.gpu_flag = false;

% initialization
ch_gain_fulltap = zeros(num.num_subc_bw, num_sim);
ch_gain_fulltap_usr = zeros(num.num_subc_usr, num_sim);
ch_gain_onetap = zeros(num.num_subc_bw, num_sim);
ch_gain_onetap_usr = zeros(num.num_subc_usr, num_sim);

% generate channel
switch idx_fading
    case {12, 13}
%         fading_ch = comm.RicianChannel(...
%             'SampleRate', num.sample_rate,...
%             'PathDelays', ch.path_delays,...
%             'AveragePathGains', ch.average_path_gains,...
%             'KFactor', ch.k_factor,...
%             'DirectPathDopplerShift', ch.maximum_doppler_shift,...
%             'MaximumDopplerShift', ch.maximum_doppler_shift,...
%             'DopplerSpectrum', ch.doppler_spectrum,...
%             'PathGainsOutputPort', true, ...
%             'NormalizePathGains', true);
        
        fading_ch = nrTDLChannel( ...
            'DelayProfile', 'TDL-E', ...
            'DelaySpread', delay_spread_rms_us*1e-6, ...
            'MaximumDopplerShift', ch.maximum_doppler_shift, ...
            'NumTransmitAntennas', 1, ...
            'NumReceiveAntennas', 1, ...
            'SampleRate', num.sample_rate, ...
            'KFactorScaling', true, ...
            'KFactor', 10*log10(ch.k_factor),...
            'NormalizePathGains', true);
    otherwise
%         fading_ch = comm.RayleighChannel(...
%             'SampleRate', num.sample_rate, ...
%             'PathDelays', ch.path_delays, ...
%             'AveragePathGains', ch.average_path_gains, ...
%             'NormalizePathGains', true, ...
%             'MaximumDopplerShift', ch.maximum_doppler_shift, ...
%             'DopplerSpectrum', ch.doppler_spectrum, ...
%             'PathGainsOutputPort', true);

        fading_ch = nrTDLChannel( ...
            'DelayProfile', 'TDL-C', ...
            'DelaySpread', delay_spread_rms_us*1e-6, ...
            'MaximumDopplerShift', ch.maximum_doppler_shift, ...
            'NumTransmitAntennas', 1, ...
            'NumReceiveAntennas', 1, ...
            'SampleRate', num.sample_rate, ...
            'NormalizePathGains', true);

%         fading_ch = nrTDLChannel( ...
%             'DelayProfile', 'Custom', ...
%             'AveragePathGains', ch.average_path_gains,...
%             'PathDelays', ch.path_delays, ...
%             'MaximumDopplerShift', ch.maximum_doppler_shift, ...
%             'NumTransmitAntennas', 1, ...
%             'NumReceiveAntennas', 1, ...
%             'SampleRate', num.sample_rate, ...
%             'NormalizePathGains', true);
end

tic
% test channel
for i = 1:num_sim
    tx_ofdmsym_serial = randn((num.num_fft+num.num_cp)*num.num_ofdmsym, 1);
    [~, ch_pg] = fading_ch(tx_ofdmsym_serial);
    [~, ch_fulltap_tf, ch_onetap_usr_tf, ch_eff_usr_tf, ch_eff_usr_dd] = gen_real_ch_r1(fading_ch, ch_pg, num.num_fft, num.num_cp, num.num_subc_bw, num.num_ofdmsym, [], [], false, test_option);
    ch_onetap_usr_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(ch_onetap_usr_tf, [], 1), [], 2);
    
    % mean channel gain
    ch_gain_fulltap(:, i) = mean(sum(abs(ch_fulltap_tf).^2, 2), 3);
    ch_gain_fulltap_usr(:, i) = mean(sum(abs(ch_fulltap_tf(list_subc_usr, list_subc_usr, list_ofdmsym_usr)).^2, 2), 3);
    
    % mean one-tap channel gain
    ch_gain_onetap_ofdmsym = zeros(num.num_subc_bw, num.num_ofdmsym);
    for j = 1:num.num_ofdmsym
        ch_gain_onetap_ofdmsym(:, j) = abs(diag(ch_fulltap_tf(:, :, j))).^2;
    end
    ch_gain_onetap(:, i) = mean(ch_gain_onetap_ofdmsym, 2);
    
    ch_gain_onetap_usr_ofdmsym = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
    for j = list_ofdmsym_usr
        ch_gain_onetap_usr_ofdmsym(:, j) = abs(diag(ch_fulltap_tf(list_subc_usr, list_subc_usr, j))).^2;
    end
    ch_gain_onetap_usr(:, i) = mean(ch_gain_onetap_usr_ofdmsym, 2);
end
toc

% check result
ch_gain_fulltap_avg = mean(ch_gain_fulltap, 2);
ch_gain_onetap_avg = mean(ch_gain_onetap, 2);
figure, plot(1:num.num_subc_bw, ch_gain_fulltap_avg, '-b.'), hold on, plot(1:num.num_subc_bw, ch_gain_onetap_avg, '-r.'), hold off, grid minor, xlabel('subcarriers'), ylabel('average channel gain'), title('Average channel gain per subcarrier'), legend('full-tap', 'one-tap')

% check result
ch_gain_fulltap_usr_avg = mean(ch_gain_fulltap_usr, 2);
ch_gain_onetap_usr_avg = mean(ch_gain_onetap_usr, 2);
figure, plot(list_subc_usr, ch_gain_fulltap_usr_avg, '-b.'), hold on, plot(list_subc_usr, ch_gain_onetap_usr_avg, '-r.'), hold off, grid minor, xlabel('subcarriers'), ylabel('average channel gain'), title('Average channel gain per subcarrier'), legend('full-tap', 'one-tap')

% % dump
% assignin('base', 'fading_ch', fading_ch)
% assignin('base', 'ch_pg', ch_pg)
% assignin('base', 'ch_full', ch_fulltap_tf)
% assignin('base', 'ch_onetap_tf', ch_onetap_usr_tf)
% assignin('base', 'ch_onetap_dd', ch_onetap_usr_dd)
% assignin('base', 'ch_eff_dd', ch_eff_usr_dd)
