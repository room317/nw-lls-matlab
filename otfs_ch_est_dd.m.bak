function ch_est_dd_ndft = otfs_ch_est_dd(rx_sym_dd_ndft, num)

% demap pilot symbols
[~, pilot_sym] = otfs_sym_demap_r1(rx_sym_dd_ndft, num, 2);
pilot_sym_upper = pilot_sym(:, 1:ceil(num.num_pilot_doppler/2));
pilot_sym_lower = pilot_sym(:, end-floor(num.num_pilot_doppler/2)+1:end);

% estimate channels
ch_est_dd_ndft_raw = zeros(num.ndft, num.num_ofdmsym_per_subframe);
scale_ch_est_dd = sqrt((num.ndft*num.num_ofdmsym_per_subframe)/((num.num_pilot_delay+num.num_guard_delay)*num.num_pilot_doppler));
% ch_est_dd_ndft_raw(1:num.num_pilot_delay, 1:ceil(num.num_pilot_doppler/2)) = pilot_sym_upper*scale_ch_est_dd;
% ch_est_dd_ndft_raw(1:num.num_pilot_delay, end-floor(num.num_pilot_doppler/2)+1:end) = pilot_sym_lower*scale_ch_est_dd;
ch_est_buffer = 0;
ch_est_dd_ndft_raw(1:num.num_pilot_delay-ch_est_buffer, 1:ceil(num.num_pilot_doppler/2)) = pilot_sym_upper(1:end-ch_est_buffer, :)*scale_ch_est_dd;
ch_est_dd_ndft_raw(1:num.num_pilot_delay-ch_est_buffer, end-floor(num.num_pilot_doppler/2)+1:end) = pilot_sym_lower(1:end-ch_est_buffer, :)*scale_ch_est_dd;
% ch_est_dd_ndft_raw(1:num.num_guard_delay+num.num_pilot_delay, 1:ceil(num.num_pilot_doppler/2)) = pilot_sym_upper*scale_ch_est_dd;
% ch_est_dd_ndft_raw(1:num.num_guard_delay+num.num_pilot_delay, end-floor(num.num_pilot_doppler/2)+1:end) = pilot_sym_lower*scale_ch_est_dd;

% output
ch_est_dd_ndft = ch_est_dd_ndft_raw;

% % smoothe channel (moving average)
% window_size = [3, 1];
% poly_b = (1/prod(window_size))*ones(window_size);
% % window_coef = [ ...
% %     1 1 1;
% %     1 10 1;
% %     1 1 1];
% % poly_b = window_coef./sum(window_coef(:));
% ch_est_dd_ndft_raw_shift = fftshift(ch_est_dd_ndft_raw, 2);
% ch_est_dd_ndft_smooth_shift = filter2(poly_b, ch_est_dd_ndft_raw_shift);
% ch_est_dd_ndft_smooth = fftshift(ch_est_dd_ndft_smooth_shift, 2);
% 
% % output
% ch_est_dd_ndft = ch_est_dd_ndft_smooth;

% % dump variables
% assignin('base', 'rx_sym_dd_ndft', rx_sym_dd_ndft);
% assignin('base', 'pilot_sym', pilot_sym);
% assignin('base', 'ch_est_dd_ndft_raw', ch_est_dd_ndft_raw);
% assignin('base', 'ch_est_dd_ndft_smooth', ch_est_dd_ndft_smooth);
% 
% % plot
% figure
% subplot(1, 2, 1), mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft, real(ch_est_dd_ndft_raw)), title('real(raw channel estimation) in dd-domain')
% subplot(1, 2, 2), mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft, imag(ch_est_dd_ndft_raw)), title('imag(raw channel estimation) in dd-domain')
% figure
% subplot(1, 2, 1), mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft, real(ch_est_dd_ndft_smooth)), title('real(smoothed channel estimation) in dd-domain')
% subplot(1, 2, 2), mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft, imag(ch_est_dd_ndft_smooth)), title('imag(smoothed channel estimation) in dd-domain')
% pause

end
