% otfs_ch_eq_tf_r3 equalizes rx signal in time-frequency domain.
% otfs_ch_eq_tf_r3 requires effective channels only when the channel is
% real.
% rx_sym_rbs_tf_eq = otfs_ch_eq_tf_r3(rx_sym_rbs_tf, ch_est_rbs_tf, ch_est_rbs_dd, ch_real_eff_tf, noise_var, chest_option, cheq_option, test_option)
%   - rx_sym_rbs_tf_eq: rx signal after equalization
%   - rx_sym_rbs_tf: rx signal before equalization
%   - ch_est_rbs_tf: time-frequency domain channel (1-tap)
%   - ch_est_rbs_dd: delay-doppler domain channel (1-tap)
%   - ch_real_eff_tf: real effective channel (full-tap)
%   - num: numerology
%   - noise_var: noise variance
%   - chest_option: channel estimation option (real or not)
%   - cheq_option: channel equalization option (zf or mmse)
%   - test_option: test option (full-tap or 1-tap)

function rx_sym_rbs_tf_eq = otfs_ch_eq_tf_r3(rx_sym_rbs_tf, ch_est_rbs_tf, ch_est_rbs_dd, ch_real_eff_tf, num, noise_var, chest_option, cheq_option, test_option)

% equalize channel
if strcmp(cheq_option, 'tfeq_mmse')
    if strcmp(chest_option, 'real') && test_option.fulltap_eq
        % calculate full-tap mmse channel
        ch_eff_tf_mmse = ch_real_eff_tf'/(ch_real_eff_tf*ch_real_eff_tf'+noise_var*eye(num.num_subc_usr*num.num_ofdmsym_usr));
        
        % equalize channel (full-tap mmse)
        rx_sym_rbs_tf_eq_vec = ch_eff_tf_mmse*rx_sym_rbs_tf(:);
        rx_sym_rbs_tf_eq = reshape(rx_sym_rbs_tf_eq_vec, num.num_subc_usr, num.num_ofdmsym_usr);
    else
        % check time-frequency channel
        if isempty(ch_est_rbs_tf)
            % 2d inverse sfft for channel transformation
            ch_est_rbs_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(ch_est_rbs_dd, [], 2), [], 1);
        end
        
        % calculate one-tap mmse channel
        ch_est_rbs_tf_mmse = conj(ch_est_rbs_tf)./(noise_var+abs(ch_est_rbs_tf).^2);
        
        % equalize channel (one-tap mmse)
        rx_sym_rbs_tf_eq = rx_sym_rbs_tf.*ch_est_rbs_tf_mmse;
    end
elseif strcmp(cheq_option, 'tfeq_zf')
    if strcmp(chest_option, 'real') && test_option.fulltap_eq
        % equalize channel (full-tap zf)
        rx_sym_rbs_tf_eq_vec = ch_real_eff_tf\rx_sym_rbs_tf(:);
        rx_sym_rbs_tf_eq = reshape(rx_sym_rbs_tf_eq_vec, num.num_subc_usr, num.num_ofdmsym_usr);
    else
        % check time-frequency channel
        if isempty(ch_est_rbs_tf)
            % 2d inverse sfft for channel transformation
            ch_est_rbs_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(ch_est_rbs_dd, [], 2), [], 1);
        end
        
        % equalize channel (one-tap zf)
        rx_sym_rbs_tf_eq = rx_sym_rbs_tf./ch_est_rbs_tf;
    end
else
    error('cheq_option value must be one of these: {tfeq_zf, tfeq_mmse}')
end

end
