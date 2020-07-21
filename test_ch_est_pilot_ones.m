%% set parameters
num_delay = 600;
num_doppler = 14;
idx_doppler_ctr = floor(num_doppler/2)+1;
idx_delay_pilot = 268:333;
len_delay_pilot = length(idx_delay_pilot);

%% generate real channel

% ch_real = complex(randn(num_delay, num_doppler), randn(num_delay, num_doppler));
ch_real = ch_real_rbs_dd/sqrt(num_delay*num_doppler);

sub_ch = zeros(num_delay, num_delay, num_doppler);
for i = 1:num_doppler
    sub_ch(:, :, i) = toeplitz(ch_real(:, i), circshift(ch_real(end:-1:1, i), 1));
end
idx_ch_order = reshape(toeplitz(1:num_doppler, circshift(num_doppler:-1:1, 1)), [], 1);
new_ch_reorder1 = reshape(permute(sub_ch(:, :, idx_ch_order), [1 3 2]), num_delay*num_doppler*num_doppler, num_delay);
new_ch_reorder2 = reshape(permute(reshape(new_ch_reorder1, num_delay*num_doppler, num_doppler, num_delay), [1 3 2]), num_delay*num_doppler, num_delay*num_doppler);
new_ch_real = new_ch_reorder2;

%% get real tx signal and extract tx pilot signal

x_real = tx_sym_rbs_dd;
m_real = toeplitz(x_real(:, idx_doppler_ctr), circshift(flipud(x_real(:, idx_doppler_ctr)), 1));   % known

x_pilot = zeros(num_delay, num_doppler); x_pilot(idx_delay_pilot, idx_doppler_ctr) = x_real(idx_delay_pilot, idx_doppler_ctr);                  % known
m_pilot = toeplitz(x_pilot(:, idx_doppler_ctr), circshift(flipud(x_pilot(:, idx_doppler_ctr)), 1));   % known

figure, subplot(2, 1, 1), plot(real(x_real(:)), '-b.'), hold on, plot(real(x_pilot(:)), ':r.'), hold off, subplot(2, 1, 2), plot(imag(x_real(:)), '-b.'), hold on, plot(imag(x_pilot(:)), ':r.'), hold off

figure
subplot(2, 2, 1), mesh(1:size(m_real, 2), 1:size(m_real, 1), real(m_real)), subplot(2, 2, 2), mesh(1:size(m_real, 2), 1:size(m_real, 1), imag(m_real))
subplot(2, 2, 3), mesh(1:size(m_pilot, 2), 1:size(m_pilot, 1), real(m_pilot)), subplot(2, 2, 4), mesh(1:size(m_pilot, 2), 1:size(m_pilot, 1), imag(m_pilot))

%% get real rx signal and compare with generated rx signal

y_real = new_ch_real*x_real(:);                                        % known
y_rx = rx_sym_rbs_dd;

figure, subplot(2, 1, 1), plot(real(y_rx(:)), '-b.'), hold on, plot(real(y_real(:)), ':r.'), hold off, subplot(2, 1, 2), plot(imag(y_rx(:)), '-b.'), hold on, plot(imag(y_real(:)), ':r.'), hold off

%% generate pilot matrix and generate rx signal and compare with real rx signal

% y1 = m*ch(:, idx_doppler_ctr);
y_pilot_est = m_pilot*circshift(ch_real, [0 -(idx_doppler_ctr-1)]);
y_pilot = zeros(size(y_rx)); y_pilot(idx_delay_pilot, :) = y_rx(idx_delay_pilot, :);                  % known
y_pilot_interp = interp2((1:num_doppler)', [1 idx_delay_pilot num_delay], ...                       % reduce to zeros
    y_pilot([1 idx_delay_pilot num_delay], 1:num_doppler), (1:num_doppler)', 1:num_delay, 'linear');

figure
subplot(4, 2, 1), mesh(1:size(y_rx, 2), 1:size(y_rx, 1), real(y_rx)), subplot(4, 2, 2), mesh(1:size(y_rx, 2), 1:size(y_rx, 1), imag(y_rx))
subplot(4, 2, 3), mesh(1:size(y_pilot, 2), 1:size(y_pilot, 1), real(y_pilot)), subplot(4, 2, 4), mesh(1:size(y_pilot, 2), 1:size(y_pilot, 1), imag(y_pilot))
subplot(4, 2, 5), mesh(1:size(y_pilot_interp, 2), 1:size(y_pilot_interp, 1), real(y_pilot_interp)), subplot(4, 2, 6), mesh(1:size(y_pilot_interp, 2), 1:size(y_pilot_interp, 1), imag(y_pilot_interp))
subplot(4, 2, 7), mesh(1:size(y_pilot_est, 2), 1:size(y_pilot_est, 1), real(y_pilot_est)), subplot(4, 2, 8), mesh(1:size(y_pilot_est, 2), 1:size(y_pilot_est, 1), imag(y_pilot_est))

%% estimate channel

n_pilot = pinv(m_pilot);

ch_est = n_pilot*circshift(y_pilot_interp, [0 (idx_doppler_ctr-1)]);
% ch_est = n_pilot*circshift(y_pilot, [0 (idx_doppler_ctr-1)]);

figure
subplot(2, 2, 1), mesh(1:size(ch_real, 2), 1:size(ch_real, 1), real(ch_real)), subplot(2, 2, 2), mesh(1:size(ch_real, 2), 1:size(ch_real, 1), imag(ch_real))
subplot(2, 2, 3), mesh(1:size(ch_est, 2), 1:size(ch_est, 1), real(ch_est)), subplot(2, 2, 4), mesh(1:size(ch_est, 2), 1:size(ch_est, 1), imag(ch_est))

%% generate rx signal with estimated channel

sub_ch = zeros(num_delay, num_delay, num_doppler);
for i = 1:num_doppler
    sub_ch(:, :, i) = toeplitz(ch_est(:, i), circshift(ch_est(end:-1:1, i), 1));
end
idx_ch_order = reshape(toeplitz(1:num_doppler, circshift(num_doppler:-1:1, 1)), [], 1);
new_ch_reorder1 = reshape(permute(sub_ch(:, :, idx_ch_order), [1 3 2]), num_delay*num_doppler*num_doppler, num_delay);
new_ch_reorder2 = reshape(permute(reshape(new_ch_reorder1, num_delay*num_doppler, num_doppler, num_delay), [1 3 2]), num_delay*num_doppler, num_delay*num_doppler);
new_ch_est = new_ch_reorder2;

y_est = new_ch_est*x_real(:);

figure, subplot(2, 1, 1), plot(real(y_rx(:)), '-b.'), hold on, plot(real(y_est), ':r.'), hold off, subplot(2, 1, 2), plot(imag(y_rx(:)), '-b.'), hold on, plot(imag(y_est), ':r.'), hold off

%% reduce pilot matrix size

m_pilot_reduced = m_pilot(idx_delay_pilot, [1:len_delay_pilot num_delay-(len_delay_pilot-1)+1:num_delay]);
n_pilot_reduced = pinv(m_pilot_reduced);
y_pilot_reduced = y_rx(idx_delay_pilot, :);
ch_est_reduced = n_pilot_reduced*circshift(y_pilot_reduced, [0 (idx_doppler_ctr-1)]);

figure
subplot(2, 2, 1), mesh(1:size(ch_est, 2), 1:size(ch_est, 1), real(ch_est)), subplot(2, 2, 2), mesh(1:size(ch_est, 2), 1:size(ch_est, 1), imag(ch_est))
subplot(2, 2, 3), mesh(1:size(ch_est_reduced, 2), 1:size(ch_est_reduced, 1), real(ch_est_reduced)), subplot(2, 2, 4), mesh(1:size(ch_est_reduced, 2), 1:size(ch_est_reduced, 1), imag(ch_est_reduced))



