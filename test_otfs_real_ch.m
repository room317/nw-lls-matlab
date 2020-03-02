% real channel for comparison with estimated channel
%     - compares real channel with estimated channel in tf and dd domain
%     - real channel is calculated in tf domain
%     - real tf channel is converted to real dd channel

function [ch_est_mse_tf, ch_est_mse_dd] = test_otfs_real_ch(ch_in_time_serial, ch_out_time_serial, ch_est_dd_ndft_shift, ch_est_tf_ndft, num)

% demap real channel input
ch_in_time = reshape(ch_in_time_serial, num.nfft+num.num_cp, num.num_ofdmsym_per_subframe);
ch_in_tf_nfft = (1/sqrt(num.nfft)) * fft(ch_in_time(num.num_cp+1:end,:), [], 1);
ch_in_tf_nfft_shift = fftshift(ch_in_tf_nfft, 1);
ch_in_tf_subc = ch_in_tf_nfft_shift(num.nfft/2-num.num_subcarrier/2+1:num.nfft/2+num.num_subcarrier/2, :);
ch_in_tf_ndft = ch_in_tf_subc;

% demap real channel output
ch_out_time = reshape(ch_out_time_serial, num.nfft+num.num_cp, num.num_ofdmsym_per_subframe);
ch_out_tf_nfft = (1/sqrt(num.nfft)) * fft(ch_out_time(num.num_cp+1:end,:), [], 1);
ch_out_tf_nfft_shift = fftshift(ch_out_tf_nfft, 1);
ch_out_tf_subc = ch_out_tf_nfft_shift(num.nfft/2-num.num_subcarrier/2+1:num.nfft/2+num.num_subcarrier/2, :);
ch_out_tf_ndft = ch_out_tf_subc;

% estimate channel
real_ch_tf_ndft = ch_out_tf_ndft ./ ch_in_tf_ndft;  % estimate channel

% 2d inverse sfft
real_ch_dd_ndft = sqrt(num.ndft/num.num_ofdmsym_per_subframe)*fft(ifft(real_ch_tf_ndft, [], 1), [], 2);
real_ch_dd_ndft_shift = fftshift(fftshift(real_ch_dd_ndft, 1), 2);

% channel mse
ch_est_mse_tf = sqrt(mean(abs(real_ch_tf_ndft(:) - ch_est_tf_ndft(:)).^2));
ch_est_mse_dd = sqrt(mean(abs(real_ch_dd_ndft_shift(:) - ch_est_dd_ndft_shift(:)).^2));

% dump variables
assignin('base', 'ch_in_time', ch_in_time);
assignin('base', 'ch_in_tf_nfft', ch_in_tf_nfft);
assignin('base', 'ch_in_tf_subc', ch_in_tf_subc);
assignin('base', 'ch_in_tf_ndft', ch_in_tf_ndft);
assignin('base', 'ch_out_time', ch_out_time);
assignin('base', 'ch_out_tf_nfft', ch_out_tf_nfft);
assignin('base', 'ch_out_tf_ndft', ch_out_tf_ndft);
assignin('base', 'real_ch_tf_ndft', real_ch_tf_ndft);
% assignin('base', 'real_ch_tf_ndft', real_ch_tf_ndft);
% assignin('base', 'real_ch_dd_ndft_shift', real_ch_dd_ndft_shift);
% assignin('base', 'ch_est_tf_ndft', ch_est_tf_ndft);
% assignin('base', 'ch_est_dd_ndft_shift', ch_est_dd_ndft_shift);

% plot results
% figure, mesh(1:14, 1:600, abs(real_ch_tf_ndft))         % tf real channel
% figure, mesh(1:14, 1:600, abs(ch_est_tf_ndft))          % tf est. channel
% figure, mesh(1:14, 1:600, abs(real_ch_dd_ndft_shift))   % dd real channel
% figure, mesh(1:14, 1:600, abs(ch_est_dd_ndft_shift))    % dd est. channel

% figure(1)
% subplot(1, 2, 1), mesh(1:14, 1:600, abs(real_ch_dd_ndft_shift))
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Channel (500 km/h)')
% subplot(1, 2, 2), contour(abs(real_ch_dd_ndft_shift))
% xlabel('Doppler'), ylabel('Delay'), title('Channel (500 km/h)'), grid

figure(1)
subplot(1, 2, 1), mesh(1:14, 1:600, abs(real_ch_dd_ndft_shift))
xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Perfect Channel Estimation (500 km/h)')
subplot(1, 2, 2), mesh(1:14, 1:600, abs(ch_est_dd_ndft_shift))
xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Practical Channel Estimation (500 km/h)')

figure(2)
subplot(1, 2, 1), mesh(1:14, 1:600, abs(real_ch_tf_ndft))
xlabel('Time (Symbols)'), ylabel('Frequency (Subcarriers)'), zlabel('Amplitude'), title('Perfect Channel Estimation (500 km/h)')
subplot(1, 2, 2), mesh(1:14, 1:600, abs(ch_est_tf_ndft))
xlabel('Time (Symbols)'), ylabel('Frequency (Subcarriers)'), zlabel('Amplitude'), title('Practical Channel Estimation (500 km/h)')

assignin('base', 'ch_est_dd_ndft_shift', ch_est_dd_ndft_shift);
assignin('base', 'ch_est_tf_ndft', ch_est_tf_ndft);
pause

end
