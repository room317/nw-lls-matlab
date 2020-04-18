function [rx_sym_tf_ndft_eq, ch_est_ndft] = ofdm_ch_comp_tf_mmse_perfect(rx_sym_tf_ndft, ch_in_time_serial, ch_out_time_serial, num, noise_var, chest_option, ch_real)

% demap data symbol
rx_sym_tf_ndft_data = rx_sym_tf_ndft(:, num.idx_data_ofdmsym);

if strcmp(chest_option, 'real')
    ch_est_ndft = ch_real;
else
    % demap channel input
    ch_in_time = reshape(ch_in_time_serial, num.nfft+num.num_cp, num.num_ofdmsym_per_subframe);
    ch_in_tf_nfft = (1/sqrt(num.nfft)) * fft(ch_in_time(num.num_cp+1:end,num.idx_data_ofdmsym), [], 1);
    ch_in_tf_nfft_shift = fftshift(ch_in_tf_nfft, 1);
    ch_in_tf_ndft = ch_in_tf_nfft_shift(num.nfft/2-num.ndft/2+1:num.nfft/2+num.ndft/2, :);
    
    % demap channel output
    ch_out_time = reshape(ch_out_time_serial, num.nfft+num.num_cp, num.num_ofdmsym_per_subframe);
    ch_out_tf_nfft = (1/sqrt(num.nfft)) * fft(ch_out_time(num.num_cp+1:end,num.idx_data_ofdmsym), [], 1);
    ch_out_tf_nfft_shift = fftshift(ch_out_tf_nfft, 1);
    ch_out_tf_ndft = ch_out_tf_nfft_shift(num.nfft/2-num.ndft/2+1:num.nfft/2+num.ndft/2, :);
    
    % estimate channel
    ch_est_ndft = ch_out_tf_ndft ./ ch_in_tf_ndft;
end

% calculate mmse
ch_est_sat = 100;   % ch_est saturation (zero-division prevention)
ch_est_ndft_clip = complex(sign(real(ch_est_ndft)).*min(abs(real(ch_est_ndft)),ch_est_sat), ...
    sign(imag(ch_est_ndft)).*min(abs(imag(ch_est_ndft)),ch_est_sat));
ch_est_ndft_mmse = conj(ch_est_ndft_clip) ./ (noise_var+abs(ch_est_ndft_clip).^2);

% compensate channel
rx_sym_tf_ndft_eq = rx_sym_tf_ndft_data .* ch_est_ndft_mmse;

end
