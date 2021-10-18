% plot sum_ch_pwr

% load sum_ch_pwr first!
if ~exist('sum_ch_pwr', 'var') && ~exist('ch_est_rbs_dd', 'var')
    error('Load ''sum_ch_pwr'' and ''ch_est_rbs_dd'' first!')
end

% set options
option_idx_scale = false;
option_plot_thr = false;
option_chest = 'dd';

% set parameters (manually)
df = scs_khz*1e3;
nfft = nw_num.num_fft;

% set parameters (automatically)
fs = df*nfft;
Ts = 1/fs;
ncp = nfft*72/1024;
Tsym = Ts*(nfft+ncp);
nsubc = size(sum_ch_pwr, 1);
nsym = size(sum_ch_pwr, 2);
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

% plot
if option_idx_scale
    figure
    mesh(doppler_axis*1e-3, delay_axis*1e6, sqrt(fftshift(fftshift(sum_ch_pwr, 1), 2)))
    xlabel('Doppler (kHz)'), ylabel('Delay (us)'), zlabel('Channel Power'), title('Channel Impulse Response'), axis tight, grid minor
else
    figure
    mesh(doppler_axis, delay_axis, sqrt(fftshift(fftshift(sum_ch_pwr, 1), 2)))
    xlabel('Doppler'), ylabel('Delay'), zlabel('Channel Power'), title('Channel Impulse Response'), axis tight, grid minor
end

% set threshold
if option_plot_thr
    thr = 0.0025*max(sum_ch_pwr, [], 'all');     % 0.001 of max power 
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

% set dd pilot
if strcmp(option_chest, 'dd')
    p_delay = 40;       % number of pilot delay grids
    p_doppler = 14;     % number of pilot doppler grids
    
    chest_dd = zeros(size(ch_est_rbs_dd));
    chest_dd(1:(p_delay/2), 1:(p_doppler/2)) = ch_est_rbs_dd(1:(p_delay/2), 1:(p_doppler/2));
    chest_dd(1:(p_delay/2), nsym-(p_doppler/2)+1:nsym) = ch_est_rbs_dd(1:(p_delay/2), nsym-(p_doppler/2)+1:nsym);
    chest_dd(nsubc-(p_delay/2)+1:nsubc, 1:(p_doppler/2)) = ch_est_rbs_dd(nsubc-(p_delay/2)+1:nsubc, 1:(p_doppler/2));
    chest_dd(nsubc-(p_delay/2)+1:nsubc, nsym-(p_doppler/2)+1:nsym) = ch_est_rbs_dd(nsubc-(p_delay/2)+1:nsubc, nsym-(p_doppler/2)+1:nsym);
    
    ch_est_rbs_tf = sqrt(nsym/nsubc)*fft(ifft(ch_est_rbs_dd, [], 2), [], 1);
    chest_tf = sqrt(nsym/nsubc)*fft(ifft(chest_dd, [], 2), [], 1);
elseif strcmp(option_chest, 'tf')
    idx_p_sym = 4:7:nsym;
    
    ch_est_rbs_tf = sqrt(nsym/nsubc)*fft(ifft(ch_est_rbs_dd, [], 2), [], 1);
    chest_tf = zeros(size(ch_est_rbs_tf));
    for i = 1:nsubc
        chest_tf(i, :) = interp1(idx_p_sym, ch_est_rbs_tf(i, idx_p_sym), 1:nsym, 'linear', 'extrap');
    end
    chest_dd = sqrt(nsubc/nsym)*fft(ifft(chest_tf, [], 1), [], 2);
else
    ch_est_rbs_tf = zeros(size(ch_est_rbs_dd));
    chest_tf = zeros(size(ch_est_rbs_dd));
    chest_dd = zeros(size(ch_est_rbs_dd));
end

% calculate error
fprintf('channel estimation error: %10.4f\n', sqrt(mean(abs(ch_est_rbs_tf-chest_tf).^2, 'all')))

% plot
if option_idx_scale
    figure
    mesh(doppler_axis*1e-3, delay_axis*1e6, abs(fftshift(fftshift(ch_est_rbs_dd, 1), 2)))
    title({'Channel Response in Delay-Doppler Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
    xlabel('Doppler (kHz)'), ylabel('Delay (us)'), zlabel('Channel Amplitude'), axis tight, grid minor
    
    figure
    mesh(doppler_axis*1e-3, delay_axis*1e6, abs(fftshift(fftshift(chest_dd, 1), 2)))
    title({'Channel Estimation in Delay-Doppler Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
    xlabel('Doppler (kHz)'), ylabel('Delay (us)'), zlabel('Channel Amplitude'), axis tight, grid minor
    
    figure
    mesh(doppler_axis*1e3, delay_axis*1e-6, abs(ch_est_rbs_tf))
    title({'Channel Response in Frequency-Time Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
    xlabel('Time (ms)'), ylabel('Frequency (MHz)'), zlabel('Channel Amplitude'), axis tight, grid minor
    
    figure
    mesh(doppler_axis*1e3, delay_axis*1e-6, abs(chest_tf))
    title({'Channel Estimation in Frequency-Time Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
    xlabel('Time (ms)'), ylabel('Frequency (MHz)'), zlabel('Channel Amplitude'), axis tight, grid minor
else
    figure
    mesh(doppler_axis, delay_axis, abs(fftshift(fftshift(ch_est_rbs_dd, 1), 2)))
    title({'Channel Response in Delay-Doppler Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
    xlabel('Doppler'), ylabel('Delay'), zlabel('Channel Amplitude'), axis tight, grid minor
    
    figure
    mesh(doppler_axis, delay_axis, abs(fftshift(fftshift(chest_dd, 1), 2)))
    title({'Channel Estimation in Delay-Doppler Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
    xlabel('Doppler'), ylabel('Delay'), zlabel('Channel Amplitude'), axis tight, grid minor
    
    figure
    mesh(doppler_axis, delay_axis, abs(ch_est_rbs_tf))
    title({'Channel Response in Frequency-Time Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
    xlabel('Symbols'), ylabel('Subcarriers'), zlabel('Channel Amplitude'), axis tight, grid minor
    
    figure
    mesh(doppler_axis, delay_axis, abs(chest_tf))
    title({'Channel Estimation in Frequency-Time Domain'; [num2str(nsubc), ' Subcarriers (', num2str(df*1e-3), ' kHz) and ', num2str(nsym/14), ' Slots']})
    xlabel('Symbols'), ylabel('Subcarriers'), zlabel('Channel Amplitude'), axis tight, grid minor
end

% use 'xlim' and 'ylim' to save axes
