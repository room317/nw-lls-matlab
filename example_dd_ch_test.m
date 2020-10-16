% parameter
ndelay = 64;
ndoppler = 14;
ncp = 16;
fsubc = 15e3;
fc_mhz = 4e3;

% channel
delay_tdla_norm = [ ...
    0.0000 0.3819 0.4025 0.5868 0.4610 0.5375 0.6708 0.5750 ...
    0.7618 1.5375 1.8978 2.2242 2.1718 2.4942 2.5119 3.0582 ...
    4.0810 4.4579 4.5695 4.7966 5.0066 5.3043 9.6586];
pwr_tdla_db = [ ...
    -13.4 0 -2.2 -4 -6 -8.2 -9.9 -10.5 ...
    -7.5 -15.9 -6.6 -16.7 -12.4 -15.2 -10.8 -11.3 ...
    -12.7 -16.2 -18.3 -18.9 -16.6 -19.9 -29.7];
v_kmh = 120;
ds_rms = 1e-7;
pd = delay_tdla_norm*ds_rms;
pg_avg = pwr_tdla_db;
doppler_max = v_kmh*fc_mhz/1080;
fading_ch = comm.RayleighChannel( ...
    'SampleRate', fsubc*ndelay, 'PathDelays', pd, ...
    'AveragePathGains', pg_avg, 'NormalizePathGains', true, ...
    'MaximumDopplerShift', doppler_max, 'PathGainsOutputPort', true, ...
    'DopplerSpectrum', doppler('Jakes'));

% tx
tx_sym_dd = zeros(ndelay, ndoppler);
tx_sym_dd(1, 1) = sqrt(ndelay*ndoppler);
tx_sym_tf = sqrt(ndoppler/ndelay)*fft(ifft(tx_sym_dd, [], 2), [], 1);
tx_ofdmsym = sqrt(ndelay)*ifft(fftshift(tx_sym_tf, 1), [], 1);
tx_ofdmsym = tx_ofdmsym([end-ncp+1:end 1:end], :);

% ch
tx_ofdmsym_faded = fading_ch(tx_ofdmsym(:));

% rx
rx_ofdmsym = reshape(tx_ofdmsym_faded, ndelay+ncp, []);
rx_ofdmsym = rx_ofdmsym(ncp+1:end, :);
rx_sym_tf = (1/sqrt(ndelay))*fft(rx_ofdmsym, [], 1);
rx_sym_tf = fftshift(rx_sym_tf, 1);
rx_sym_dd = sqrt(ndelay/ndoppler)*fft(ifft(rx_sym_tf, [], 1), [], 2);

% regen rx
chest_dd = rx_sym_dd;
chest_dd_eff = gen_eff_ch(chest_dd);
chest_tf = sqrt(ndoppler/ndelay)*fft(ifft(chest_dd, [], 2), [], 1);
rx_sym_tf_regen = chest_tf.*tx_sym_tf;
rx_sym_dd_regen = chest_dd_eff*tx_sym_dd(:);
rx_sym_dd_regen = reshape(rx_sym_dd_regen, ndelay, []);

% check
figure
subplot(2, 1, 1), plot(real(rx_sym_tf(:)), '-b.'), hold on, plot(real(rx_sym_tf_regen(:)), ':r.'), hold off, grid minor, title('real(tf-domain rx symbol)'), legend('original rx', 'regenerated rx')
subplot(2, 1, 2), plot(imag(rx_sym_tf(:)), '-b.'), hold on, plot(imag(rx_sym_tf_regen(:)), ':r.'), hold off, grid minor, title('imag(tf-domain rx symbol)'), legend('original rx', 'regenerated rx')
figure
subplot(2, 1, 1), plot(real(rx_sym_dd(:)), '-b.'), hold on, plot(real(rx_sym_dd_regen(:)), ':r.'), hold off, grid minor, title('real(dd-domain rx symbol)'), legend('original rx', 'regenerated rx')
subplot(2, 1, 2), plot(imag(rx_sym_dd(:)), '-b.'), hold on, plot(imag(rx_sym_dd_regen(:)), ':r.'), hold off, grid minor, title('imag(dd-domain rx symbol)'), legend('original rx', 'regenerated rx')
