% otfs_ch_eq_dd_r4 equalizes rx signal in delay-doppler domain.
% otfs_ch_eq_dd_r4 always requires effective channels.
% [rx_sym_rbs_eq_dd, noise_var_mat_dd] = otfs_ch_eq_dd_r4(rx_sym_rbs_dd, ch_est_rbs_tf, ch_est_rbs_dd, ch_real_eff_dd, num, noise_var, chest_option, cheq_option, test_option)
%   - rx_sym_rbs_eq_dd: rx signal after equalization
%   - rx_sym_rbs_dd: rx signal before equalization
%   - ch_est_rbs_tf: time-frequency domain channel (1-tap)
%   - ch_est_rbs_dd: delay-doppler domain channel (1-tap)
%   - ch_real_eff_dd: real effective channel (full-tap)
%   - num: numerology
%   - noise_var: noise variance
%   - chest_option: channel estimation option (real or not)
%   - cheq_option: channel equalization option (zf or mmse)
%   - test_option: test option (full-tap or 1-tap)

function [rx_sym_rbs_eq_dd, noise_var_mat_dd] = otfs_ch_eq_dd_r4(rx_sym_rbs_dd, ch_est_rbs_tf, ch_est_rbs_dd, ch_real_eff_dd, num, noise_var, chest_option, cheq_option, test_option)

% generate effective channel
if strcmp(chest_option, 'real') && test_option.fulltap_eq
    % get real effective channel
    ch_eff_dd = ch_real_eff_dd;
else
    if isempty(ch_est_rbs_dd)
        % 2d inverse sfft for channel transformation
        ch_est_rbs_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(ch_est_rbs_tf, [], 1), [], 2);
    end
    
    % generate block circular channel matrix
    ch_eff_dd = gen_eff_ch(ch_est_rbs_dd);
end

% equalize channel in dd domain
if strcmp(cheq_option, 'ddeq_zf')
    % one-tap equalization (matrix inversion)
    rx_sym_rbs_eq_vec_dd = ch_eff_dd\rx_sym_rbs_dd(:);
    rx_sym_rbs_eq_dd = reshape(rx_sym_rbs_eq_vec_dd, num.num_delay_usr, num.num_doppler_usr);
    
%     % calculate noise variance
%     noise_var_mat_dd = noise_var*(abs(reshape(ch_eff_dd\ones(num.num_delay_usr*num.num_doppler_usr, 1), num.num_delay_usr, num.num_doppler_usr)).^2);
elseif strcmp(cheq_option, 'ddeq_mmse')
    % full-tap equalization (matrix inversion)
    ch_eff_mmse_dd = ch_eff_dd'/(ch_eff_dd*ch_eff_dd'+noise_var*eye(num.num_delay_usr*num.num_doppler_usr));
    rx_sym_rbs_eq_vec_dd = ch_eff_mmse_dd*rx_sym_rbs_dd(:);
    rx_sym_rbs_eq_dd = reshape(rx_sym_rbs_eq_vec_dd, num.num_delay_usr, num.num_doppler_usr);
    
%     % calculate noise variance
%     noise_var_mat_dd = noise_var*(abs(reshape(sum(ch_eff_mmse_dd, 2), num.num_delay_usr, num.num_doppler_usr)).^2);
else
    error('cheq_option value must be one of these: {ddeq_zf, ddeq_mmse}')
end

% calculate noise variance (common for both)
noise_var_mat_dd = reshape(noise_var*abs(sum(ch_eff_dd, 2)).^2, num.num_delay_usr, num.num_doppler_usr);

end
