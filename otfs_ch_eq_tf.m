function rx_sym_tf_ndft_eq = otfs_ch_eq_tf(rx_sym_tf_ndft, ch_est_tf_ndft, noise_var, cheq_option)

if strcmp(cheq_option, 'tfeq')
    % equalize channel
    rx_sym_tf_ndft_eq = rx_sym_tf_ndft ./ ch_est_tf_ndft;
elseif strcmp(cheq_option, 'tfeq_mmse')
    % calculate mmse channel
    ch_est_ndft_mmse = conj(ch_est_tf_ndft) ./ (noise_var+abs(ch_est_tf_ndft).^2);
    
    % equalize channel
    rx_sym_tf_ndft_eq = rx_sym_tf_ndft .* ch_est_ndft_mmse;
else
    rx_sym_tf_ndft_eq = [];
end

end
