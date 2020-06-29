function rx_sym_rbs_tf_eq = otfs_ch_eq_tf_r2(rx_sym_rbs_tf, ch_est_rbs_tf, noise_var, cheq_option)

if strcmp(cheq_option, 'tfeq_mmse')
    
    % calculate mmse channel
    ch_est_rbs_tf_mmse = conj(ch_est_rbs_tf) ./ (noise_var+abs(ch_est_rbs_tf).^2);
    
    % equalize channel
    rx_sym_rbs_tf_eq = rx_sym_rbs_tf .* ch_est_rbs_tf_mmse;
    
elseif strcmp(cheq_option, 'tfeq_zf')
    
    % equalize channel
    rx_sym_rbs_tf_eq = rx_sym_rbs_tf ./ ch_est_rbs_tf;
    
else
    error('Use ''otfs_ch_eq_tf_r2'' for other ''cheq_option'' values.')
end

end
