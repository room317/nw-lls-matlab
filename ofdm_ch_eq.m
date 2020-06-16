function rx_sym_rbs_eq = ofdm_ch_eq(rx_sym_rbs, ch_est_rbs, noise_var, cheq_option)

if strcmp(cheq_option, 'tfeq_zf')
    % equalize channel
    rx_sym_rbs_eq = rx_sym_rbs ./ ch_est_rbs;
elseif strcmp(cheq_option, 'tfeq_mmse')
    % calculate mmse channel
    ch_est_ndft_mmse = conj(ch_est_rbs) ./ (noise_var+abs(ch_est_rbs).^2);
    
    % equalize channel
    rx_sym_rbs_eq = rx_sym_rbs .* ch_est_ndft_mmse;
else
    error('''cheq_option'' value must be one of these: {''tfeq_zf'', ''tfeq_mmse''}')
end

end
