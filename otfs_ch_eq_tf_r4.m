% otfs_ch_eq_tf_r4 equalizes rx signal in time-frequency domain.
% otfs_ch_eq_tf_r4 requires effective channels only when the channel is
% real.
% [rx_sym_rbs_eq_tf, noise_var_mat_dd] = otfs_ch_eq_tf_r4(rx_sym_rbs_tf, ch_est_rbs_tf, ch_est_rbs_dd, ch_real_eff_tf, noise_var, chest_option, cheq_option, test_option)
%   - rx_sym_rbs_eq_tf: rx signal after equalization
%   - rx_sym_rbs_tf: rx signal before equalization
%   - ch_est_rbs_tf: time-frequency domain channel (1-tap)
%   - ch_est_rbs_dd: delay-doppler domain channel (1-tap)
%   - ch_real_eff_tf: real effective channel (full-tap)
%   - num: numerology
%   - noise_var: noise variance
%   - chest_option: channel estimation option (real or not)
%   - cheq_option: channel equalization option (zf or mmse)
%   - test_option: test option (full-tap or 1-tap)

function [rx_sym_rbs_eq_tf, noise_var_mat_dd] = otfs_ch_eq_tf_r4(rx_sym_rbs_tf, ch_est_rbs_tf, ch_est_rbs_dd, ch_real_eff_tf, ch_real_eff_dd, num, noise_var, chest_option, cheq_option, test_option)

% equalize channel
if strcmp(cheq_option, 'tfeq_zf')
    if strcmp(chest_option, 'real') && test_option.fulltap_eq
        % equalize channel (full-tap zf)
        rx_sym_rbs_eq_vec_tf = ch_real_eff_tf\rx_sym_rbs_tf(:);
        rx_sym_rbs_eq_tf = reshape(rx_sym_rbs_eq_vec_tf, num.num_subc_usr, num.num_ofdmsym_usr);
        
%         % calculate noise variance
%         ch_inv_dd = reshape(sum(inv(ch_real_eff_dd), 2), num.num_delay_usr, num.num_doppler_usr);
%         noise_var_mat_dd = noise_var*(abs(ch_inv_dd).^2);
    else
        % check time-frequency channel
        if isempty(ch_est_rbs_tf)
            % 2d inverse sfft for channel transformation
            ch_est_rbs_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(ch_est_rbs_dd, [], 2), [], 1);
        end
        
        % equalize channel (one-tap zf)
        rx_sym_rbs_eq_tf = rx_sym_rbs_tf./ch_est_rbs_tf;
        
%         % calculate noise variance
%         if isempty(ch_est_rbs_dd)
%             % 2d inverse sfft for channel transformation
%             ch_est_rbs_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(ch_est_rbs_tf, [], 1), [], 2);
%         end
% %         ch_inv_dd = reshape(sum(inv(gen_eff_ch(ch_est_rbs_dd)), 2), num.num_delay_usr, num.num_doppler_usr);
% %         ch_inv_dd = sum(ch_est_rbs_dd, 'all')*ones(num.num_delay_usr, num.num_doppler_usr)/sqrt(numel(ch_est_rbs_dd));  % simplified
% %         noise_var_mat_dd = noise_var*(abs(ch_inv_dd).^2);
% %         noise_var_mat_dd = noise_var/var(ch_est_rbs_dd, 1, 'all');
%         
%         
% %         ch_dd = gen_eff_ch(ch_est_rbs_dd);
%         ch_inv_dd = inv(gen_eff_ch(ch_est_rbs_dd));
% %         assignin('base', 'ch_dd', ch_dd)
% %         assignin('base', 'ch_inv_dd', ch_inv_dd)
% %         pause
%         noise_var_mat_dd = noise_var*sum(abs(ch_inv_dd).^2, 2);
    end
elseif strcmp(cheq_option, 'tfeq_mmse')
    if strcmp(chest_option, 'real') && test_option.fulltap_eq
        % calculate full-tap mmse channel
        ch_eff_mmse_tf = ch_real_eff_tf'/(ch_real_eff_tf*ch_real_eff_tf'+noise_var*eye(num.num_subc_usr*num.num_ofdmsym_usr));
        
        % equalize channel (full-tap mmse)
        rx_sym_rbs_eq_vec_tf = ch_eff_mmse_tf*rx_sym_rbs_tf(:);
        rx_sym_rbs_eq_tf = reshape(rx_sym_rbs_eq_vec_tf, num.num_subc_usr, num.num_ofdmsym_usr);
        
%         % calculate noise variance
%         ch_eff_mmse_dd = sfft_mtx*ch_eff_mmse_tf*isfft_mtx;
%         ch_inv_dd = reshape(sum(ch_eff_mmse_dd, 2), num.num_delay_usr, num.num_doppler_usr);
%         noise_var_mat_dd = noise_var*(abs(ch_inv_dd).^2);
    else
        % check time-frequency channel
        if isempty(ch_est_rbs_tf)
            % 2d inverse sfft for channel transformation
            ch_est_rbs_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(ch_est_rbs_dd, [], 2), [], 1);
        end
        
        % calculate one-tap mmse channel
        ch_est_rbs_mmse_tf = conj(ch_est_rbs_tf)./(noise_var+abs(ch_est_rbs_tf).^2);
        
        % equalize channel (one-tap mmse)
        rx_sym_rbs_eq_tf = rx_sym_rbs_tf.*ch_est_rbs_mmse_tf;
        
%         % calculate noise variance
%         ch_est_rbs_mmse_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(ch_est_rbs_mmse_tf, [], 1), [], 2);
%         ch_inv_dd = reshape(sum(gen_eff_ch(ch_est_rbs_mmse_dd), 2), num.num_delay_usr, num.num_doppler_usr);
%         noise_var_mat_dd = noise_var*(abs(ch_inv_dd).^2);
    end
else
    error('cheq_option value must be one of these: {tfeq_zf, tfeq_mmse}')
end

% calculate noise variance (common for both)
if strcmp(chest_option, 'real') && test_option.fulltap_eq
%     noise_var_mat_dd = noise_var*reshape(abs(sum(ch_eff_dd, 2)).^2, num.num_subc_usr, []);
    noise_var_mat_dd = noise_var*reshape(sum(abs(inv(ch_real_eff_dd)).^2, 2), num.num_subc_usr, []);
else
    noise_var_mat_dd = noise_var*sum(abs(1./ch_est_rbs_tf).^2, 'all')/(num.num_doppler_usr*num.num_delay_usr)*ones(num.num_delay_usr, num.num_doppler_usr);
end

% assignin('base', 'rx_sym_rbs_tf', rx_sym_rbs_tf)
% assignin('base', 'ch_real_eff_dd', ch_real_eff_dd)
% assignin('base', 'ch_real_eff_tf', ch_real_eff_tf)
% assignin('base', 'ch_est_rbs_dd', ch_est_rbs_dd)
% assignin('base', 'ch_est_rbs_tf', ch_est_rbs_tf)
% assignin('base', 'noise_var', noise_var)
% assignin('base', 'noise_var_mat_dd', noise_var_mat_dd)
% pause

end
