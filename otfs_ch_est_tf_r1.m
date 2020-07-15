function ch_est_rbs_tf = otfs_ch_est_tf_r1(tx_sym_rbs_tf, ch_out_time_serial, num, chest_option, fading_ch, ch_path_gain)

if strcmp(chest_option, 'real')
    
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
        ch_real_halfmap_mat_tf(:, :, idx_ch_doppler) = ifft(fft(ch_real_mat_t(:, :, idx_ch_doppler), [], 1), [], 2);    % circulant matrix svd
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
    
else
    error('Use ''otfs_ch_est_tf_r1'' for other ''chest_option'' values.')
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
