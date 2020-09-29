% ofdm_ch_eq_r1 equalizes rx ofdm signal in time-frequency domain.
% ofdm_ch_eq_r1 requires effective channels only when the channel is real.
% rx_sym_rbs_eq = ofdm_ch_eq_r1(rx_sym_rbs, ch_est_rbs, ch_real_eff, num, noise_var, chest_option, cheq_option)
%   - rx_sym_rbs_eq: rx signal after equalization
%   - rx_sym_rbs: rx signal before equalization
%   - ch_est_rbs: time-frequency domain channel (1-tap)
%   - ch_real_eff: real effective channel (full-tap)
%   - num: numerology
%   - noise_var: noise variance
%   - chest_option: channel estimation option (real or not)
%   - cheq_option: channel equalization option (zf or mmse)
%   - test_option: test option (full-tap or 1-tap)

function rx_sym_rbs_eq = ofdm_ch_eq_r1(rx_sym_rbs, ch_est_rbs, ch_real_eff, num, noise_var, chest_option, cheq_option, test_option)

if strcmp(cheq_option, 'tfeq_zf')
    if strcmp(chest_option, 'real') && test_option.fulltap_eq
        % equalize channel (full-tap zf)
        rx_sym_rbs_eq_vec = ch_real_eff\rx_sym_rbs(:);
        rx_sym_rbs_eq = reshape(rx_sym_rbs_eq_vec, num.num_subc_usr, num.num_ofdmsym_usr);
    else
        % equalize channel (one-tap zf)
        rx_sym_rbs_eq = rx_sym_rbs./ch_est_rbs;
    end
elseif strcmp(cheq_option, 'tfeq_mmse')
    if strcmp(chest_option, 'real') && test_option.fulltap_eq
        % calculate full-tap mmse channel
        ch_real_eff_mmse = ch_real_eff'/(ch_real_eff*ch_real_eff'+noise_var*eye(num.num_subc_usr*num.num_ofdmsym_usr));
        
        % equalize channel (full-tap mmse)
        rx_sym_rbs_eq_vec = ch_real_eff_mmse*rx_sym_rbs(:);
        rx_sym_rbs_eq = reshape(rx_sym_rbs_eq_vec, num.num_subc_usr, num.num_ofdmsym_usr);
    else
        % calculate mmse channel
        ch_est_rbs_mmse = conj(ch_est_rbs)./(noise_var+abs(ch_est_rbs).^2);
        
        % equalize channel (one-tap)
        rx_sym_rbs_eq = rx_sym_rbs.*ch_est_rbs_mmse;
    end
else
    error('''cheq_option'' value must be one of these: {''tfeq_zf'', ''tfeq_mmse''}')
end

end
