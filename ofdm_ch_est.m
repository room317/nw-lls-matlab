% frequency domain channel estimation
%   - ch_in_time_serial : for perfert or real ch. est.
%   - ch_out_time_serial: for perfert or real ch. est.
%   - tx_sym_pilot_ndft : for pilot ch. est.
%   - rx_sym_ndft       : for pilot ch. est.


function ch_est_rbs = ofdm_ch_est(tx_sym_rbs, rx_sym_rbs, ch_out_time_serial, num, chest_option, fading_ch, ch_path_gain)

if strcmp(chest_option, 'tf_ltedown') || strcmp(chest_option, 'tf_lteup') || strcmp(chest_option, 'tf_nr')
    
    % set parameters
    ofdmsym_avg_win_size = 1;   % ofdm symbol average window size for channel estimation
    
    rx_sym_pilot_rbs = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
    rx_ofdmsym_pilot_interp_1d = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
    rx_ofdmsym_pilot_interp = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
    for i = 1:num.num_ofdmsym_pilot_usr
        % demap pilots
        rx_sym_pilot_rbs(num.idx_subc_pilot_usr(:, i), num.idx_ofdmsym_pilot_usr(i)) = ...
            rx_sym_rbs(num.idx_subc_pilot_usr(:, i), num.idx_ofdmsym_pilot_usr(i)) .* ...
            conj(tx_sym_rbs(num.idx_subc_pilot_usr(:, i), num.idx_ofdmsym_pilot_usr(i)));
        
        % find symbol index within averaging window
        idx_ofdmsym_pilot_avg = ...
            num.idx_ofdmsym_pilot_usr( ...
            num.idx_ofdmsym_pilot_usr >= max(1, num.idx_ofdmsym_pilot_usr(i)-ofdmsym_avg_win_size) & ...
            num.idx_ofdmsym_pilot_usr <= num.idx_ofdmsym_pilot_usr(i));
        num_ofdmsym_pilot_avg = length(idx_ofdmsym_pilot_avg);
        
        % average pilots across symbols
        cnt_ofdmsym_pilot_avg = ones(num.num_subc_usr, 1);
        cnt_ofdmsym_pilot_avg(num.idx_subc_pilot_usr(:, i), 1) = ...
            sum(double(ismember(num.idx_subc_pilot_usr(:, i-num_ofdmsym_pilot_avg+1:i), num.idx_subc_pilot_usr(:, i))), 2);
        mask_subc_avg = zeros(num.num_subc_usr, 1);
        mask_subc_avg(num.idx_subc_pilot_usr(:, i), 1) = ones(num.num_subc_pilot_usr, 1);
        rx_ofdmsym_pilot_avg = sum(rx_sym_pilot_rbs(:, idx_ofdmsym_pilot_avg), 2) .* mask_subc_avg ./ cnt_ofdmsym_pilot_avg;
        
        % generate virtual pilots (ignored)
        rx_ofdmsym_pilot_vp = rx_ofdmsym_pilot_avg;
        
        % interpolate 1st dimension (along subcarriers)
        rx_ofdmsym_pilot_interp_1d(:, num.idx_ofdmsym_pilot_usr(i)) = ...
            interp1(num.idx_subc_pilot_usr(:, i), rx_ofdmsym_pilot_vp(num.idx_subc_pilot_usr(:, i), 1), ...
            1:num.num_subc_usr, 'linear', 'extrap');
        
    end
    
    % interpolate 2nd dimension (along ofdm symbols)
    for i = 1:num.num_subc_usr
        rx_ofdmsym_pilot_interp(i, :) = ...
            interp1(num.idx_ofdmsym_pilot_usr, rx_ofdmsym_pilot_interp_1d(i, num.idx_ofdmsym_pilot_usr), ...
            1:num.num_ofdmsym_usr, 'linear', 'extrap');
    end
    
    % output
    ch_est_rbs = rx_ofdmsym_pilot_interp;
    
elseif strcmp(chest_option, 'real')
    
    % get channel info
    ch_info = info(fading_ch);
    ch_filter_coeff = ch_info.ChannelFilterCoefficients;    % ch_filter_coeff: channel filter coefficient per path
    num_ch_path = length(fading_ch.PathDelays);                % num_path: number of path
    
    % reproduce real channel (time saving)
    len_ch_filter = size(ch_filter_coeff, 2);
    ch_coeff_mat = zeros(num.nfft, num.nfft, num_ch_path);
    ch_real_mat_t = zeros(num.nfft, num.nfft, num.num_ofdmsym_usr);
    ch_real_halfmap_mat_tf = zeros(num.nfft, num.nfft, num.num_ofdmsym_usr);
    
    idx_ofdmsym_usr = 1:num.num_ofdmsym_usr;
    ch_path_gain_reshape = reshape(ch_path_gain, num.nfft+num.num_cp, num.num_ofdmsym_subfrm, num_ch_path);
    ch_path_gain_reshape(1:num.num_cp, :, :) = [];
    ch_path_gain_reshape = ch_path_gain_reshape(:, idx_ofdmsym_usr, :);
    for idx_ch_path = 1:num_ch_path
        ch_coeff_row = circshift([ch_filter_coeff(idx_ch_path, end:-1:1) zeros(1, num.nfft-len_ch_filter)], -len_ch_filter+1);
        ch_coeff_col = [ch_filter_coeff(idx_ch_path, :) zeros(1, num.nfft-len_ch_filter)];
        ch_coeff_mat(:, :, idx_ch_path) = toeplitz(ch_coeff_col, ch_coeff_row);
    end
    for idx_ch_sym = 1:num.num_ofdmsym_usr
        path_gain_mat = ch_path_gain_reshape(:, idx_ch_sym, :);
        ch_mat_per_path = path_gain_mat .* ch_coeff_mat;
        ch_real_mat_t(:, :, idx_ch_sym) = sum(ch_mat_per_path, 3);
        ch_real_halfmap_mat_tf(:, :, idx_ch_sym) = ifft(fft(ch_real_mat_t(:, :, idx_ch_sym), [], 1), [], 2);
    end
    ch_real_mat_shift_tf = fftshift(fftshift(ch_real_halfmap_mat_tf, 2), 1);
    ch_real_mat_tf = ...
        ch_real_mat_shift_tf((num.nfft/2)-(num.num_subc_usr/2)+1:(num.nfft/2)+(num.num_subc_usr/2), ...
        (num.nfft/2)-(num.num_subc_usr/2)+1:(num.nfft/2)+(num.num_subc_usr/2), :);
    
    % extract diagonal elements of real channel
    ch_real_rbs = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
    for idx_ch_sym = 1:num.num_ofdmsym_usr
        ch_real_rbs(:, idx_ch_sym) = diag(ch_real_mat_tf(:, :, idx_ch_sym));  % diagonal term only
    end
    
%     figure
%     subplot(1, 2, 1), mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, real(ch_real_rbs))
%     subplot(1, 2, 2), mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, imag(ch_real_rbs))
%     figure
%     subplot(1, 2, 1), mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, real(ch_est_rbs))
%     subplot(1, 2, 2), mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, imag(ch_est_rbs))
%     figure
%     subplot(1, 2, 1), mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, real(ch_est_rbs-ch_real_rbs))
%     subplot(1, 2, 2), mesh(1:num.num_ofdmsym_usr, 1:num.num_subc_usr, imag(ch_est_rbs-ch_real_rbs))
%     fprintf('rmse:%10.4f\n', sqrt(mean((ch_est_rbs(:)-ch_real_rbs(:)).^2)))
%     pause
    
    % output
    ch_est_rbs = ch_real_rbs;
    
elseif strcmp(chest_option, 'perfect')
    
    % demap channel input
    ch_in_sym_rbs = tx_sym_rbs;
    
    % demap channel output
    ch_out_ofdmsym = reshape(ch_out_time_serial, num.nfft+num.num_cp, num.num_ofdmsym_subfrm);
    ch_out_sym_nfft = (1/sqrt(num.nfft)) * fft(ch_out_ofdmsym(num.num_cp+1:end, :), [], 1);
    ch_out_sym_nfft_shift = fftshift(ch_out_sym_nfft, 1);
    ch_out_sym_bw = ch_out_sym_nfft_shift(num.nfft/2-num.num_subc_usr/2+1:num.nfft/2+num.num_subc_usr/2, :);
    ch_out_sym_rbs = ch_out_sym_bw;
    
    % estimate channels
    ch_est_rbs_zf = ch_out_sym_rbs ./ ch_in_sym_rbs;  % estimate channel
    ch_est_sat = 100;   % ch_est saturation (zero-division prevention)
    ch_est_rbs = complex(sign(real(ch_est_rbs_zf)).*min(abs(real(ch_est_rbs_zf)),ch_est_sat), ...
        sign(imag(ch_est_rbs_zf)).*min(abs(imag(ch_est_rbs_zf)),ch_est_sat));
    
else
    error('''chest_option'' value must be one of these: {''tf_ltedown'', ''tf_lteup'', ''tf_nr'', ''perfect'', ''real''}')
end

end
