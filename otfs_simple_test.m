% otfs_simple test code

% time domain noise check
nt = rx_ofdmsym_serial-tx_ofdmsym_faded;
fprintf('time domain noise variance: %8.4f\n', mean(abs(nt).^2, 'all'))

% tf domain signal check
x = tx_sym_rbs_tf(:);
y = rx_sym_rbs_tf(:);
h = ch_real_eff;
hx = h*x;
n = y-hx;
fprintf('tf domain noise variance: %8.4f\n', mean(abs(n).^2, 'all'))
figure, plot(real(y), '-b.'), hold on, plot(real(h*x), ':r.'), grid, legend('real(y)', 'real(hx)'), hold off
figure, plot(real(y), '-b'), hold on, plot(real(hx+n), ':r'), grid, legend('real(y)', 'real(hx+n)'), hold off

% zf signal check (wy = whx+wn)
w = 1./ch_est_rbs_tf(:);
wy = w.*y;
whx = w.*hx;
wn = w.*n;
xest = h\y;
figure, plot(real(wy), '-b'), hold on, plot(real(whx+wn), ':r'), grid, legend('real(wy)', 'real(whx+wn)'), hold off

% tf domain estimation check (wy = whx+wn)
figure
subplot(1, 3, 1), mesh(1:num_doppler, 1:num_delay, abs(reshape(xest, num_delay, num_doppler))), legend('abs(xest)'), axis_ref = axis;
subplot(1, 3, 2), mesh(1:num_doppler, 1:num_delay, abs(reshape(wy, num_delay, num_doppler))), legend('abs(wy)'), axis(axis_ref)
subplot(1, 3, 3), mesh(1:num_doppler, 1:num_delay, abs(reshape(x, num_delay, num_doppler))), legend('abs(x)'), axis(axis_ref)
figure, plot(real(x), '-k'), hold on, plot(real(xest), '-b'), plot(real(wy), ':r'), grid, legend('real(x)', 'real(xest):perfect', 'real(wy):zf'), hold off
fprintf('tx avg sig pwr: %8.4f\n', mean(abs(x).^2, 'all'))
fprintf('zf rx avg sig pwr: %8.4f\n', mean(abs(wy).^2, 'all'))
fprintf('est tx avg sig pwr: %8.4f\n', mean(abs(xest).^2, 'all'))

% dd domain signal check (kswy = kswhx+kswn)
wy_dd = sqrt(132/56)*sfft_mtx*wy;
whx_dd = sqrt(132/56)*sfft_mtx*diag(w)*h*x;
wn_dd = sqrt(132/56)*sfft_mtx*diag(w)*n;
figure
subplot(1, 2, 1), scatter(real(wy_dd), imag(wy_dd)), grid, legend('zf rx signal with noise')
subplot(1, 2, 2), scatter(real(whx_dd), imag(whx_dd)), grid, legend('zf rx signal without noise')
figure, plot(real(wy_dd), '-b'), hold on, plot(real(whx_dd+wn_dd), ':r'), grid, legend('wy(dd)', 'whx(dd)+wn(dd)'), hold off
fprintf('dd signal variance: %8.4f\n', mean(abs(whx_dd).^2, 'all'))
fprintf('dd noise variance: %8.4f\n', mean(abs(wn_dd).^2, 'all'))
figure, plot(abs(whx_dd).^2, '-b'), hold on, plot(abs(wn_dd).^2, ':r'), grid, legend('whx(dd) domain power', 'wn(dd) domain power'), hold off
figure, plot(real(wy_dd), '-b'), hold on, plot(real(whx_dd+wn_dd), ':r'), grid, legend('real(wy(dd))', 'real(whx(dd)+wn(dd))'), hold off

% kswy = kswhx+kswn
y_dd = sqrt(132/56)*sfft_mtx/h*y;
x_dd = sqrt(132/56)*sfft_mtx*x;
n_dd = sqrt(132/56)*sfft_mtx/h*n;
figure, plot(real(y_dd), '-b'), hold on, plot(real(x_dd+n_dd), ':r'), grid, legend('real(y(dd))', 'real(x(dd)+n(dd))'), hold off
figure, plot(abs(x_dd).^2, '-b'), hold on, plot(abs(n_dd).^2, ':r'), grid, legend('abs(x(dd))', 'abs(n(dd))'), hold off

% noise variance check
n_var_ref = abs(n_dd).^2;
n_var_est1 = noise_var*mean(sum(abs(sqrt(132/56)*sfft_mtx/h).^2, 2), 'all');
n_var_est2 = noise_var*mean(sum(abs(inv(h)).^2, 2), 'all');
n_var_est3 = noise_var*mean(abs(w).^2, 'all');
fprintf('noise var ref: %8.4f\n', mean(n_var_ref, 'all'))
fprintf('noise var est1: %8.4f\n', n_var_est1)
fprintf('noise var est2: %8.4f\n', n_var_est2)
fprintf('noise var est3: %8.4f\n', n_var_est3)
wn_est = noise_var*(abs(w).^2);
figure, plot(abs(wn).^2, '-b'), hold on, plot(wn_est, ':r'), grid, legend('wn var', 'wn var est'), hold off
figure, plot(n_var_ref, '-b'), hold on, plot(n_var_est1*ones(num_delay*num_doppler, 1), ':r'), grid, legend('n var ref', 'n var est1'), hold off
figure, plot(n_var_ref, '-b'), hold on, plot(n_var_est2*ones(num_delay*num_doppler, 1), ':r'), grid, legend('n var ref', 'n var est2'), hold off
figure, plot(n_var_ref, '-b'), hold on, plot(n_var_est3*ones(num_delay*num_doppler, 1), ':r'), grid, legend('n var ref', 'n var est3'), hold off
