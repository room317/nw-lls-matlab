function rx_sym_rbs_dd_eq = otfs_ch_eq_dd_r2(rx_sym_rbs_dd, ch_est_rbs_dd, num)

% generate block circular channel matrix
eff_ch = gen_eff_ch(ch_est_rbs_dd);

% equalize channel in dd domain
rx_sym_rbs_dd_eq_vec = eff_ch \ rx_sym_rbs_dd(:);
rx_sym_rbs_dd_eq = reshape(rx_sym_rbs_dd_eq_vec, num.num_delay_usr, num.num_doppler_usr);

% % dump variables
% assignin('base', 'rx_sym_rbs_dd', rx_sym_rbs_dd);
% assignin('base', 'ch_est_rbs_dd', ch_est_rbs_dd);
% assignin('base', 'new_ch', new_ch);
% assignin('base', 'rx_sym_rbs_dd_eq', rx_sym_rbs_dd_eq);

end
