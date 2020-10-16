function [rx_sym_dd_ndft_eq, ch_est_dd_ndft, ch_est_tf_ndft] = otfs_ch_comp_dd_mmse_r1(rx_sym_dd_ndft, rx_sym_tf_ndft, idx_pilot_sym, num, noise_var, chest_option, ch_dd_real)

% demap pilot symbols
[~, pilot_sym] = otfs_sym_demap_r1(rx_sym_dd_ndft, num);
pilot_sym_upper = pilot_sym(:, 1:ceil(num.num_pilot_doppler/2));
pilot_sym_lower = pilot_sym(:, end-floor(num.num_pilot_doppler/2)+1:end);

% estimate channels
ch_est_dd_ndft_raw = zeros(num.ndft, num.num_ofdmsym_per_subframe);
scale_ch_est_dd = sqrt((num.ndft*num.num_ofdmsym_per_subframe)/(num.num_pilot_delay*num.num_pilot_doppler));
ch_est_dd_ndft_raw(1:num.num_pilot_delay, 1:ceil(num.num_pilot_doppler/2)) = pilot_sym_upper*scale_ch_est_dd;
ch_est_dd_ndft_raw(1:num.num_pilot_delay, end-floor(num.num_pilot_doppler/2)+1:end) = pilot_sym_lower*scale_ch_est_dd;

% smoothe channel (moving average)
window_size = [7, 3];
ch_est_dd_ndft_raw_shift = fftshift(ch_est_dd_ndft_raw, 2);
poly_b = (1/prod(window_size))*ones(window_size);
ch_est_dd_ndft_smooth_shift = filter2(poly_b, ch_est_dd_ndft_raw_shift);
ch_est_dd_ndft_smooth = fftshift(ch_est_dd_ndft_smooth_shift, 2);









% offset pilot (to prevent data penetration)
pilot_offset = [20 1]; % [20 1];

% estimate channels
    ch_est_dd_ndft = zeros(num.ndft, num.num_ofdmsym_per_subframe);
    ch_est_scale = sqrt((num.ndft*num.num_ofdmsym_per_subframe)/(num.num_pilot_delay*num.num_pilot_doppler));
    ch_est_dd_ndft(idx_pilot_sym(1)+1+pilot_offset(1):idx_pilot_sym(1)+num.num_pilot_delay-pilot_offset(1), idx_pilot_sym(2)+1+pilot_offset(2):idx_pilot_sym(2)+num.num_pilot_doppler-pilot_offset(2)) = ch_est_scale * pilot_sym;
    
    % linear extrapolation along symbols
    idx_interp_subc = idx_pilot_sym(1)+1+pilot_offset(1):idx_pilot_sym(1)+num.num_pilot_delay-pilot_offset(1);
    idx_interp_lower_sym_a = idx_pilot_sym(2)+1+pilot_offset(2);
    idx_interp_lower_sym_b = idx_pilot_sym(2)+1+pilot_offset(2)+1;
    interp_lower_sym_a = ch_est_dd_ndft(idx_interp_subc, idx_interp_lower_sym_a);
    interp_lower_sym_b = ch_est_dd_ndft(idx_interp_subc, idx_interp_lower_sym_b);
    ch_est_dd_ndft(idx_interp_subc, 1:idx_interp_lower_sym_a-1) = lin_interp2([interp_lower_sym_a interp_lower_sym_b], [idx_interp_lower_sym_a idx_interp_lower_sym_b], 1:idx_interp_lower_sym_a-1, 2);
    
    idx_interp_upper_sym_a = idx_pilot_sym(2)+num.num_pilot_doppler-pilot_offset(2)-1;
    idx_interp_upper_sym_b = idx_pilot_sym(2)+num.num_pilot_doppler-pilot_offset(2);
    interp_lower_sym_a = ch_est_dd_ndft(idx_interp_subc, idx_interp_upper_sym_a);
    interp_lower_sym_b = ch_est_dd_ndft(idx_interp_subc, idx_interp_upper_sym_b);
    ch_est_dd_ndft(idx_interp_subc, idx_interp_upper_sym_b+1:end) = lin_interp2([interp_lower_sym_a interp_lower_sym_b], [idx_interp_upper_sym_a idx_interp_upper_sym_b], idx_interp_upper_sym_b+1:num.num_ofdmsym_per_subframe, 2);
    
    % extrapolation along subcarriers
%     % 1. linear interpolation (extend the last two samples)
%     idx_interp_lower_subc_a = idx_pilot_sym(1)+1+pilot_offset(1);
%     idx_interp_lower_subc_b = idx_pilot_sym(1)+1+pilot_offset(1)+1;
%     interp_lower_subc_a = ch_est_dd_ndft(idx_interp_lower_subc_a, :);
%     interp_lower_subc_b = ch_est_dd_ndft(idx_interp_lower_subc_b, :);
%     ch_est_dd_ndft(1:idx_interp_lower_subc_a-1, :) = lin_interp2([interp_lower_subc_a; interp_lower_subc_b], [idx_interp_lower_subc_a idx_interp_lower_subc_b], 1:idx_interp_lower_subc_a-1, 1);
%     
%     idx_interp_upper_subc_a = idx_pilot_sym(1)+num.num_pilot_delay-pilot_offset(1)-1;
%     idx_interp_upper_subc_b = idx_pilot_sym(1)+num.num_pilot_delay-pilot_offset(1);
%     interp_upper_subc_a = ch_est_dd_ndft(idx_interp_upper_subc_a, :);
%     interp_upper_subc_b = ch_est_dd_ndft(idx_interp_upper_subc_b, :);
%     ch_est_dd_ndft(idx_interp_upper_subc_b+1:end, :) = lin_interp2([interp_upper_subc_a; interp_upper_subc_b], [idx_interp_upper_subc_a idx_interp_upper_subc_b], idx_interp_upper_subc_b+1:num.ndft, 1);
    
    % 2. linear extrapolation (converge to zero)
    idx_interp_lower_subc_a = 1;
    idx_interp_lower_subc_b = idx_pilot_sym(1)+1+pilot_offset(1);
    interp_lower_subc_a = ch_est_dd_ndft(idx_interp_lower_subc_a, :);
    interp_lower_subc_b = ch_est_dd_ndft(idx_interp_lower_subc_b, :);
    ch_est_dd_ndft(1:idx_interp_lower_subc_b-1, :) = lin_interp2([interp_lower_subc_a; interp_lower_subc_b], [idx_interp_lower_subc_a idx_interp_lower_subc_b], 1:idx_interp_lower_subc_b-1, 1);
    
    idx_interp_upper_subc_a = idx_pilot_sym(1)+num.num_pilot_delay-pilot_offset(1);
    idx_interp_upper_subc_b = num.ndft;
    interp_upper_subc_a = ch_est_dd_ndft(idx_interp_upper_subc_a, :);
    interp_upper_subc_b = ch_est_dd_ndft(idx_interp_upper_subc_b, :);
    ch_est_dd_ndft(idx_interp_upper_subc_a+1:end, :) = lin_interp2([interp_upper_subc_a; interp_upper_subc_b], [idx_interp_upper_subc_a idx_interp_upper_subc_b], idx_interp_upper_subc_a+1:num.ndft, 1);
    
%     % 3. hybrid extrapolation (dft after linear interpolation)
%     interp_step = 8;
%     idx_interp_lower_subc_a = idx_pilot_sym(1)+1+pilot_offset(1);
%     idx_interp_lower_subc = 1 : interp_step : num.ndft/2;
%     if mod(length(idx_interp_lower_subc), 2) ~= 0
%         idx_interp_lower_subc(end) = [];
%     end
%     interp_lower_subc = ch_est_dd_ndft(idx_interp_lower_subc, :);
%     lower_subc = dft_interp(interp_lower_subc, interp_step*length(idx_interp_lower_subc), 1, 0);
%     ch_est_dd_ndft(1:idx_interp_lower_subc_a-1, :) = lower_subc(1:idx_interp_lower_subc_a-1, :);
%     
%     idx_interp_upper_subc_a = idx_pilot_sym(1)+num.num_pilot_delay-pilot_offset(1);
%     idx_interp_upper_subc_b = num.ndft-idx_interp_upper_subc_a;
%     idx_interp_upper_subc_rev = num.ndft-interp_step+1 : -interp_step : num.ndft/2+1;
%     idx_interp_upper_subc = idx_interp_upper_subc_rev(end:-1:1);
%     if mod(length(idx_interp_upper_subc), 2) ~= 0
%         idx_interp_upper_subc(1) = [];
%     end
%     interp_upper_subc = ch_est_dd_ndft(idx_interp_upper_subc, :);
%     upper_subc = dft_interp(interp_upper_subc, interp_step*length(idx_interp_upper_subc), 1, 0);
%     ch_est_dd_ndft(idx_interp_upper_subc_a+1:end, :) = upper_subc(end-idx_interp_upper_subc_b+1:end, :);

ch_est_dd_ndft_shift = fftshift(fftshift(ch_est_dd_ndft, 1), 2);

% % 4. exponential extrapolation
% smooth_filter = [1 2 3 2 1]/sqrt(19);  % remove white noise
% idx_subc_lower_inter = 3 : idx_pilot_sym(1)+1+pilot_offset(1);
% idx_subc_lower_extra = 3 : num.ndft/2;
% for idx_sym = 1 : num.num_pilot_doppler
%     conv()
% end
% aaa = ch_est_scale*fftshift(fftshift(rx_sym_dd_ndft, 1), 2);
% assignin('base', 'aaa', aaa);

% 2d sfft (otfs transform, from delay-doppler domain to freq-time domain)
ch_est_tf_ndft = sqrt(num.num_ofdmsym_per_subframe/num.ndft)*fft(ifft(ch_est_dd_ndft_shift, [], 2), [], 1);
ch_est_tf_ndft_mmse = conj(ch_est_tf_ndft) ./ (noise_var+abs(ch_est_tf_ndft).^2);

% assignin('base', 'ch_est_dd_ndft_shift', ch_est_dd_ndft_shift);
% assignin('base', 'ch_est_dd_ndft', ch_est_dd_ndft);

% compensate channel in tf domain
rx_sym_tf_ndft_eq = rx_sym_tf_ndft .* ch_est_tf_ndft_mmse;

% 2d inverse sfft
rx_sym_dd_ndft_eq = sqrt(num.ndft/num.num_ofdmsym_per_subframe)*fft(ifft(rx_sym_tf_ndft_eq, [], 1), [], 2);

% % dump variables
% assignin('base', 'rx_sym_dd_ndft', rx_sym_dd_ndft);
% assignin('base', 'pilot_sym', pilot_sym);
% assignin('base', 'ch_est_dd_ndft', ch_est_dd_ndft);
% assignin('base', 'ch_est_dd_ndft_shift', ch_est_dd_ndft_shift);
% assignin('base', 'ch_est_tf_ndft', ch_est_tf_ndft);
% assignin('base', 'ch_est_tf_ndft_mmse', ch_est_tf_ndft_mmse);

end
