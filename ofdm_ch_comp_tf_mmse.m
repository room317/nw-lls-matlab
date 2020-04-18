function [rx_sym_tf_ndft_eq, ch_est_ndft] = ofdm_ch_comp_tf_mmse(rx_sym_tf_ndft, tx_sym_tf_ndft_pilot, num, noise_var)

% demap data and pilot symbol
rx_sym_tf_ndft_data = rx_sym_tf_ndft(:, num.idx_data_ofdmsym);
rx_sym_tf_ndft_pilot = rx_sym_tf_ndft(:, num.idx_pilot_ofdmsym);

% estimate channel
% ch_est_ndft = ch_out_tf_ndft ./ ch_in_tf_ndft;
ch_est_ndft = rx_sym_tf_ndft_pilot .* tx_sym_tf_ndft_pilot;
ch_est_sat = 100;   % ch_est saturation (zero-division prevention)
ch_est_ndft_clip = complex(sign(real(ch_est_ndft)).*min(abs(real(ch_est_ndft)),ch_est_sat), ...
    sign(imag(ch_est_ndft)).*min(abs(imag(ch_est_ndft)),ch_est_sat));
ch_est_ndft_mmse = conj(ch_est_ndft_clip) ./ (noise_var+abs(ch_est_ndft_clip).^2);

% interpolate channels
% ch_est_ndft_mmse_interp = dft_interp(ch_est_ndft_mmse_pilot, size(rx_sym_tf_ndft, 2), 2, 3);
% ch_est_ndft_mmse_interp = lin_interp(ch_est_ndft_mmse_pilot, num.idx_pilot_ofdmsym, nw_num.idx_data_ofdmsym, 2);
ch_est_ndft_mmse_interp = lin_interp2(ch_est_ndft_mmse, num.idx_pilot_ofdmsym, num.idx_data_ofdmsym, 2);

% compensate channel
rx_sym_tf_ndft_eq = rx_sym_tf_ndft_data .* ch_est_ndft_mmse_interp;

end
