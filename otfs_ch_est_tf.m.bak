function ch_est_ndft = otfs_ch_est_tf(ch_in_time_serial, ch_out_time_serial, num, chest_option)

% set symbol index
if strcmp(chest_option, 'perfect') || strcmp(chest_option, 'real')
    idx_sym = 1:num.num_ofdmsym_per_subframe;
elseif strcmp(chest_option, 'tf_pilot')
    idx_sym = num.idx_pilot_ofdmsym;
else
    idx_sym = [];
end

% demap channel input
ch_in_time = reshape(ch_in_time_serial, num.nfft+num.num_cp, num.num_ofdmsym_per_subframe);
ch_in_tf_nfft = (1/sqrt(num.nfft)) * fft(ch_in_time(num.num_cp+1:end, idx_sym), [], 1);
ch_in_tf_nfft_shift = fftshift(ch_in_tf_nfft, 1);
ch_in_tf_subc = ch_in_tf_nfft_shift(num.nfft/2-num.num_subcarrier/2+1:num.nfft/2+num.num_subcarrier/2, :);
ch_in_tf_ndft = ch_in_tf_subc;

% demap channel output
ch_out_time = reshape(ch_out_time_serial, num.nfft+num.num_cp, num.num_ofdmsym_per_subframe);
ch_out_tf_nfft = (1/sqrt(num.nfft)) * fft(ch_out_time(num.num_cp+1:end, idx_sym), [], 1);
ch_out_tf_nfft_shift = fftshift(ch_out_tf_nfft, 1);
ch_out_tf_subc = ch_out_tf_nfft_shift(num.nfft/2-num.num_subcarrier/2+1:num.nfft/2+num.num_subcarrier/2, :);
ch_out_tf_ndft = ch_out_tf_subc;

% estimate channels
ch_est_ndft = ch_out_tf_ndft ./ ch_in_tf_ndft;  % estimate channel
ch_est_sat = 100;   % ch_est saturation (zero-division prevention)
ch_est_ndft_clip = complex(sign(real(ch_est_ndft)).*min(abs(real(ch_est_ndft)),ch_est_sat), ...
    sign(imag(ch_est_ndft)).*min(abs(imag(ch_est_ndft)),ch_est_sat));

% interpolate channels
if strcmp(chest_option, 'perfect') || strcmp(chest_option, 'real')
    ch_est_ndft = ch_est_ndft_clip;
elseif strcmp(chest_option, 'tf_pilot')
    % ch_est_ndft_mmse_interp = dft_interp(ch_est_ndft_mmse, size(rx_sym_tf_ndft, 2), 2, 3);
    % ch_est_ndft_mmse_interp = lin_interp(ch_est_ndft_mmse_pilot, num.idx_pilot_ofdmsym, 1:size(rx_sym_tf_ndft, 2), 2);
    ch_est_ndft_interp = lin_interp2(ch_est_ndft_clip, num.idx_pilot_ofdmsym, 1:num.num_ofdmsym_per_subframe, 2);
    ch_est_ndft = ch_est_ndft_interp;
else
    ch_est_ndft = [];
end

end
