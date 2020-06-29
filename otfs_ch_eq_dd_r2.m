function rx_sym_rbs_dd_eq = otfs_ch_eq_dd_r2(rx_sym_rbs_dd, ch_est_rbs_dd, num)

% generate block circular channel matrix
%   ex) reshape(permute(reshape(reshape(permute(a(:, :, [1 3 2 4]), [1 3 2]), 8, 3), 4, 2, 3), [1 3 2]), 4, 6)
sub_ch = zeros(num.num_delay_usr, num.num_delay_usr, num.num_doppler_usr);
for i = 1:num.num_doppler_usr
    sub_ch(:, :, i) = toeplitz(ch_est_rbs_dd(:, i), circshift(ch_est_rbs_dd(end:-1:1, i), 1));
end
idx_ch_order = reshape(toeplitz(1:num.num_doppler_usr, circshift(num.num_doppler_usr:-1:1, 1)), [], 1);
new_ch_reorder1 = reshape(permute(sub_ch(:, :, idx_ch_order), [1 3 2]), num.num_delay_usr*num.num_doppler_usr*num.num_doppler_usr, num.num_delay_usr);
new_ch_reorder2 = reshape(permute(reshape(new_ch_reorder1, num.num_delay_usr*num.num_doppler_usr, num.num_doppler_usr, num.num_delay_usr), [1 3 2]), num.num_delay_usr*num.num_doppler_usr, num.num_delay_usr*num.num_doppler_usr);
new_ch = new_ch_reorder2/sqrt(num.num_subc_usr*num.num_ofdmsym_usr);

% equalize channel in dd domain
rx_sym_rbs_dd_eq_vec = new_ch \ rx_sym_rbs_dd(:);
rx_sym_rbs_dd_eq = reshape(rx_sym_rbs_dd_eq_vec, num.num_delay_usr, num.num_doppler_usr);

% % dump variables
% assignin('base', 'rx_sym_rbs_dd', rx_sym_rbs_dd);
% assignin('base', 'ch_est_rbs_dd', ch_est_rbs_dd);
% assignin('base', 'new_ch', new_ch);
% assignin('base', 'rx_sym_rbs_dd_eq', rx_sym_rbs_dd_eq);

end
