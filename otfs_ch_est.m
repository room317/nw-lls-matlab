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

function [ch_est_rbs_dd, ch_est_rbs_tf] = otfs_ch_est(tx_sym_rbs_tf, rx_sym_rbs_dd, ch_out_time_serial, num, chest_option, fading_ch, ch_path_gain, map_plan)

if strcmp(chest_option, 'dd_tone')
    
    % demap pilot symbols
    [~, rx_sym_pilot_subfrm] = otfs_sym_demap_r2(rx_sym_rbs_dd, num, map_plan);
    
    % remove buffer
    rx_sym_pilot_buff = rx_sym_pilot_subfrm(num.num_delay_buff+1:end-num.num_delay_buff, num.num_doppler_buff+1:end-num.num_doppler_buff);
    
    % estimate channels
    idx_ch_est_delay = ceil(num.num_delay_data_usr/2)+num.num_delay_buff+1: ...
        ceil(num.num_delay_data_usr/2)+num.num_delay_pilot_usr+num.num_delay_guard_usr-num.num_delay_buff;
    idx_ch_est_doppler = ceil(num.num_doppler_data_usr/2)+num.num_doppler_buff+1: ...
        ceil(num.num_doppler_data_usr/2)+num.num_doppler_pilot_usr+num.num_doppler_guard_usr-num.num_doppler_buff;
    scale_ch_est_rbs = sqrt((num.num_delay_usr*num.num_doppler_usr)/((num.num_delay_pilot_usr+num.num_delay_guard_usr)*(num.num_doppler_pilot_usr+num.num_doppler_guard_usr)));
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
%     ch_est_rbs_dd = circshift(ch_est_rbs_interpl, [-ceil(num.num_delay_usr/2), -ceil(num.num_doppler_usr/2)]);
    ch_est_rbs_dd = circshift(ch_est_rbs_raw, [-ceil(num.num_delay_usr/2), -ceil(num.num_doppler_usr/2)]);
    ch_est_rbs_tf = [];
    
elseif strcmp(chest_option, 'real')
    
    % get channel info
    ch_info = info(fading_ch);
    ch_filter_coeff = ch_info.ChannelFilterCoefficients;    % ch_filter_coeff: channel filter coefficient per path
    num_ch_path = length(fading_ch.PathDelays);                % num_path: number of path
    
    % reproduce real channel (time saving)
    len_ch_filter = size(ch_filter_coeff, 2);
    ch_coeff_mat = zeros(num.nfft, num.nfft, num_ch_path);
    ch_real_mat_t = zeros(num.nfft, num.nfft, num.num_doppler_usr);
    ch_real_halfmap_mat_tf = zeros(num.nfft, num.nfft, num.num_doppler_usr);
    
    idx_doppler_usr = 1:num.num_doppler_usr;
    ch_path_gain_reshape = reshape(ch_path_gain, num.nfft+num.num_cp, num.num_ofdmsym_subfrm, num_ch_path);
    ch_path_gain_reshape(1:num.num_cp, :, :) = [];
    ch_path_gain_reshape = ch_path_gain_reshape(:, idx_doppler_usr, :);
    for idx_ch_path = 1:num_ch_path
        ch_coeff_row = circshift([ch_filter_coeff(idx_ch_path, end:-1:1) zeros(1, num.nfft-len_ch_filter)], -len_ch_filter+1);
        ch_coeff_col = [ch_filter_coeff(idx_ch_path, :) zeros(1, num.nfft-len_ch_filter)];
        ch_coeff_mat(:, :, idx_ch_path) = toeplitz(ch_coeff_col, ch_coeff_row);
    end
    for idx_ch_doppler = 1:num.num_doppler_usr
        path_gain_mat = ch_path_gain_reshape(:, idx_ch_doppler, :);
        ch_mat_per_path = path_gain_mat .* ch_coeff_mat;
        ch_real_mat_t(:, :, idx_ch_doppler) = sum(ch_mat_per_path, 3);
        ch_real_halfmap_mat_tf(:, :, idx_ch_doppler) = ifft(fft(ch_real_mat_t(:, :, idx_ch_doppler), [], 1), [], 2);
    end
    ch_real_mat_shift_tf = fftshift(fftshift(ch_real_halfmap_mat_tf, 2), 1);
    ch_real_mat_tf = ...
        ch_real_mat_shift_tf((num.nfft/2)-(num.num_subc_bw/2)+1:(num.nfft/2)+(num.num_subc_bw/2), ...
        (num.nfft/2)-(num.num_subc_bw/2)+1:(num.nfft/2)+(num.num_subc_bw/2), :);
    
    % extract diagonal elements of real channel
    ch_real_rbs = zeros(num.num_delay_usr, num.num_doppler_usr);
    for idx_ch_doppler = 1:num.num_doppler_usr
        ch_real_rbs(:, idx_ch_doppler) = diag(ch_real_mat_tf(1:num.num_delay_usr, 1:num.num_delay_usr, idx_ch_doppler));  % diagonal term only
    end
    
%     figure
%     subplot(1, 2, 1), mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, real(ch_real_rbs))
%     subplot(1, 2, 2), mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, imag(ch_real_rbs))
%     figure
%     subplot(1, 2, 1), mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, real(ch_est_rbs))
%     subplot(1, 2, 2), mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, imag(ch_est_rbs))
%     figure
%     subplot(1, 2, 1), mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, real(ch_est_rbs-ch_real_rbs))
%     subplot(1, 2, 2), mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, imag(ch_est_rbs-ch_real_rbs))
%     fprintf('rmse:%10.4f\n', sqrt(mean((ch_est_rbs(:)-ch_real_rbs(:)).^2)))
%     pause
    
    % output
    ch_est_rbs_tf = ch_real_rbs;
    ch_est_rbs_dd = [];
    
elseif strcmp(chest_option, 'perfect')
    
    % demap channel input
    ch_in_sym_rbs_sfft = tx_sym_rbs_tf;
    
    % demap channel output
    ch_out_ofdmsym = reshape(ch_out_time_serial, num.nfft+num.num_cp, num.num_ofdmsym_subfrm);
    ch_out_sym_nfft = (1/sqrt(num.nfft)) * fft(ch_out_ofdmsym(num.num_cp+1:end, :), [], 1);
    ch_out_sym_nfft_shift = fftshift(ch_out_sym_nfft, 1);
    ch_out_sym_bw = ch_out_sym_nfft_shift(num.nfft/2-num.num_subc_bw/2+1:num.nfft/2+num.num_subc_bw/2, :);
    ch_out_sym_rbs_sfft = ch_out_sym_bw(1:num.num_delay_usr, 1:num.num_doppler_usr);      % temporary
    
    % estimate channels
    ch_perfect_rbs_zf = ch_out_sym_rbs_sfft ./ ch_in_sym_rbs_sfft;  % estimate channel
    ch_sat = 100;   % ch_est saturation (zero-division prevention)
    ch_perfect_rbs_sat = complex(sign(real(ch_perfect_rbs_zf)).*min(abs(real(ch_perfect_rbs_zf)),ch_sat), ...
        sign(imag(ch_perfect_rbs_zf)).*min(abs(imag(ch_perfect_rbs_zf)),ch_sat));
    
    % output
    ch_est_rbs_tf = ch_perfect_rbs_sat;
    ch_est_rbs_dd = [];
    
else
    error('''chest_option'' value must be one of these: {''dd_tone'', ''perfect'', ''real''}')
end

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
