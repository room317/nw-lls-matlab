% len_buff: 2-dimensional

% map plan example (d: data, p: pilot, -: pilot, +: guard, .: buffer)
%   (11)data       :  d d d d d d d d d d d d d d d d
%   (10)data       :  d d d d d d d d d d d d d d d d
%    (9)data&pilot :  d d d d . . . . . . . . d d d d
%    (8)data&pilot :  d d d d . . . . . . . . d d d d
%    (7)data&pilot :  d d d d . . + + - - . . d d d d
%    (6)data&pilot :  d d d d . . + + p - . . d d d d
%    (5)data&guard :  d d d d . . + + + + . . d d d d
%    (4)data&guard :  d d d d . . + + + + . . d d d d
%    (3)data&guard :  d d d d . . . . . . . . d d d d
%    (2)data&guard :  d d d d . . . . . . . . d d d d
%    (1)data       :  d d d d d d d d d d d d d d d d
%    (0)data       :  d d d d d d d d d d d d d d d d

function ch_est_rbs_dd = otfs_ch_est_dd_r1(rx_sym_rbs_dd, num, map_plan)

% demap pilot symbols
[~, rx_sym_pilot_subfrm] = otfs_sym_demap_r2(rx_sym_rbs_dd, num, map_plan);

% remove buffer
rx_sym_pilot_buff = rx_sym_pilot_subfrm(num.num_delay_buff+1:end-num.num_delay_buff, num.num_doppler_buff+1:end-num.num_doppler_buff);

% estimate channels
idx_ch_est_delay = ceil(num.num_delay_data_usr/2)+num.num_delay_buff+1: ...
    ceil(num.num_delay_data_usr/2)+num.num_delay_pilot_usr+num.num_delay_guard_usr-num.num_delay_buff;
idx_ch_est_doppler = ceil(num.num_doppler_data_usr/2)+num.num_doppler_buff+1: ...
    ceil(num.num_doppler_data_usr/2)+num.num_doppler_pilot_usr+num.num_doppler_guard_usr-num.num_doppler_buff;
% scale_ch_est_rbs = sqrt((num.num_delay_usr*num.num_doppler_usr)/ ...
%     ((num.num_delay_pilot_usr+num.num_delay_guard_usr-(2*num.num_delay_buff))*(num.num_doppler_pilot_usr+num.num_doppler_guard_usr-(2*num.num_doppler_buff))));
scale_ch_est_rbs = sqrt((num.num_delay_usr*num.num_doppler_usr)/ ...
    ((num.num_delay_pilot_usr+num.num_delay_guard_usr)*(num.num_doppler_pilot_usr+num.num_doppler_guard_usr)));
ch_est_rbs_raw = zeros(num.num_delay_usr, num.num_doppler_usr);
ch_est_rbs_raw(idx_ch_est_delay, idx_ch_est_doppler) = rx_sym_pilot_buff*scale_ch_est_rbs;

% interpolate channels (linear)
if num.num_delay_data_usr > 0 || num.num_delay_buff > 0
    idx_ch_est_delay_interp = [1, idx_ch_est_delay, num.num_delay_usr];
else
    idx_ch_est_delay_interp = idx_ch_est_delay;
end
if num.num_doppler_data_usr > 0 || num.num_doppler_buff > 0
    idx_ch_est_doppler_interp = [1; idx_ch_est_doppler'; num.num_doppler_usr];
else
    idx_ch_est_doppler_interp = idx_ch_est_doppler';
end

ch_est_rbs_interpl = interp2(idx_ch_est_doppler_interp, idx_ch_est_delay_interp, ...
    ch_est_rbs_raw(idx_ch_est_delay_interp, idx_ch_est_doppler_interp), (1:num.num_doppler_usr)', 1:num.num_delay_usr);

% output
ch_est_rbs_dd = circshift(ch_est_rbs_interpl, [-ceil(num.num_delay_usr/2), -ceil(num.num_doppler_usr/2)]);
% ch_est_rbs_dd = circshift(ch_est_rbs_raw, [-ceil(num.num_delay_usr/2), -ceil(num.num_doppler_usr/2)]);

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
