% plot sum_ch_pwr

% load sum_ch_pwr first!
if ~test_option.ch_mse
    error('Set ''test_option.ch_mse'' to ''true''!')
elseif ~exist('sum_ch_pwr', 'var') || ~exist('nw_num', 'var') || ...
        ~exist('ch_est_rbs_dd', 'var') || ~exist('ch_real_eff_tf', 'var')
    error('Load variables first!')
end

% set options
option_idx_scale = false;
option_plot_thr = false;

% set parameters (manually)
df = scs_khz*1e3;
nfft = nw_num.num_fft;
nsubc = nw_num.num_subc_usr;
nsym = nw_num.num_ofdmsym_usr;
nusr = nw_num.num_usr;                         % number of users
nrb = nw_rm.N_RB;                              % number of resource blocks
p_delay = nw_num.num_delay_pilot_usr;       % number of pilot delay grids
p_doppler = nw_num.num_doppler_pilot_usr;   % number of pilot doppler grids
idx_p_sym = 4:7:nsym;                       % indices for ofdm pilot symbols
thr_value = 0.0025;                         % power threshold value

% set parameters (automatically)
fs = df*nfft;
Ts = 1/fs;
ncp = nfft*72/1024;
Tsym = Ts*(nfft+ncp);
fbw = df*nsubc;
T = Tsym*nsym;
ddelay = 1/fbw;
ddoppler = 1/T;
rdelay = ddelay*nsubc;
rdoppler = ddoppler*nsym;

% set index (scale)
if option_idx_scale
    doppler_axis = ((0:nsym-1)-(nsym/2))*ddoppler;
    delay_axis = ((0:nsubc-1)-(nsubc/2))*ddelay;
    time_axis = ((0:nsym-1)-(nsym/2))*Tsym;
    freq_axis = ((0:nsubc-1)-(nsubc/2))*df;
else
    doppler_axis = (0:nsym-1)-(nsym/2);
    delay_axis = (0:nsubc-1)-(nsubc/2);
    time_axis = (0:nsym-1)-(nsym/2);
    freq_axis = (0:nsubc-1)-(nsubc/2);
end

% plot channel impulse response
if option_idx_scale
    figure
    mesh(doppler_axis*1e-3, delay_axis*1e6, sqrt(fftshift(fftshift(sum_ch_pwr, 1), 2)))
    xlabel('Doppler (kHz)'), ylabel('Delay (us)'), zlabel('Channel Power'), title('Channel Impulse Response'), axis tight, grid minor
else
    figure
    mesh(doppler_axis, delay_axis, sqrt(fftshift(fftshift(sum_ch_pwr, 1), 2)))
    xlabel('Doppler'), ylabel('Delay'), zlabel('Channel Power'), title('Channel Impulse Response'), axis tight, grid minor
end

% plot area above threshold
if option_plot_thr
    thr = thr_value*max(sum_ch_pwr, [], 'all');     % 0.001 of max power 
    thr_ch_pwr = double(sum_ch_pwr > thr);
    
    % plot
    if option_idx_scale
        figure
        mesh(doppler_axis*1e-3, delay_axis*1e6, sqrt(fftshift(fftshift(thr_ch_pwr, 1), 2)))
        xlabel('Doppler (kHz)'), ylabel('Delay (us)'), zlabel('Channel Power > Threshold'), title('Channel Power > Threshold'), axis tight, grid minor
    else
        figure
        mesh(doppler_axis, delay_axis, sqrt(fftshift(fftshift(thr_ch_pwr, 1), 2)))
        xlabel('Doppler'), ylabel('Delay'), zlabel('Channel Power > Threshold'), title('Channel Power > Threshold'), axis tight, grid minor
    end
end

% plot channel estimation
ch_est_rbs_tf = sqrt(nsym/nsubc)*fft(ifft(ch_est_rbs_dd, [], 2), [], 1);
ch_real_rbs_tf = zeros(nsubc, nsym, nrb, nusr);

for idx_usr = 1:nusr
    for idx_rb = 1:nrb
        for idx_sym = 1:nsym
            ch_real_rbs_tf(:, idx_sym, idx_rb, idx_usr) = ...
                diag(ch_real_eff_tf(nsubc*(idx_sym-1)+1:nsubc*idx_sym, nsubc*(idx_sym-1)+1:nsubc*idx_sym, idx_rb, idx_usr));  % diagonal term only
        end
        % ch_real_eff_usr_dd = nw_num.sfft_mtx*ch_real_eff_usr_tf*nw_num.isfft_mtx;
        % ch_real_rbs_dd = gen_conv_ch(ch_real_eff_usr_dd, nsubc);
    end
end
ch_real_rbs_dd = sqrt(nsubc/nsym)*fft(ifft(ch_real_rbs_tf, [], 1), [], 2);

% calculate error
fprintf('channel estimation error: %10.4f\n', sqrt(mean(abs(ch_real_rbs_dd-ch_est_rbs_dd).^2, 'all')))

if option_idx_scale
    for idx_usr = 1:nusr
        for idx_rb = 1:nrb
            figure
            mesh(doppler_axis*1e-3, delay_axis*1e6, abs(fftshift(fftshift(ch_real_rbs_dd(:, :, idx_rb, idx_usr), 1), 2)))
            title({'Channel Response in Delay-Doppler Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
            xlabel('Doppler (kHz)'), ylabel('Delay (us)'), zlabel('Perfect Channel Amplitude'), axis tight, grid minor
            
            figure
            mesh(doppler_axis*1e-3, delay_axis*1e6, abs(fftshift(fftshift(ch_est_rbs_dd(:, :, idx_rb, idx_usr), 1), 2)))
            title({'Channel Estimation in Delay-Doppler Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
            xlabel('Doppler (kHz)'), ylabel('Delay (us)'), zlabel('Estimated Channel Amplitude'), axis tight, grid minor
            
            figure
            mesh(doppler_axis*1e3, delay_axis*1e-6, abs(ch_real_rbs_tf(:, :, idx_rb, idx_usr)))
            title({'Channel Response in Frequency-Time Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
            xlabel('Time (ms)'), ylabel('Frequency (MHz)'), zlabel('Perfect Channel Amplitude'), axis tight, grid minor
            
            figure
            mesh(doppler_axis*1e3, delay_axis*1e-6, abs(ch_est_rbs_tf(:, :, idx_rb, idx_usr)))
            title({'Channel Estimation in Frequency-Time Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
            xlabel('Time (ms)'), ylabel('Frequency (MHz)'), zlabel('Estimated Channel Amplitude'), axis tight, grid minor
        end
    end
else
    for idx_usr = 1:nusr
        for idx_rb = 1:nrb
            figure
            mesh(doppler_axis, delay_axis, abs(fftshift(fftshift(ch_real_rbs_dd(:, :, idx_rb, idx_usr), 1), 2)))
            title({'Channel Response in Delay-Doppler Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
            xlabel('Doppler'), ylabel('Delay'), zlabel('Perfect Channel Amplitude'), axis tight, grid minor
            
            figure
            mesh(doppler_axis, delay_axis, abs(fftshift(fftshift(ch_est_rbs_dd(:, :, idx_rb, idx_usr), 1), 2)))
            title({'Channel Estimation in Delay-Doppler Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
            xlabel('Doppler'), ylabel('Delay'), zlabel('Estimated Channel Amplitude'), axis tight, grid minor
            
            figure
            mesh(doppler_axis, delay_axis, abs(ch_real_rbs_tf(:, :, idx_rb, idx_usr)))
            title({'Channel Response in Frequency-Time Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
            xlabel('Symbols'), ylabel('Subcarriers'), zlabel('Perfect Channel Amplitude'), axis tight, grid minor
            
            figure
            mesh(doppler_axis, delay_axis, abs(ch_est_rbs_tf(:, :, idx_rb, idx_usr)))
            title({'Channel Estimation in Frequency-Time Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
            xlabel('Symbols'), ylabel('Subcarriers'), zlabel('Estimated Channel Amplitude'), axis tight, grid minor
        end
    end
end

% use 'xlim' and 'ylim' to save axes
