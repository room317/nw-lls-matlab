% map otfs data and pilot qam symbols

function rx_sym_dd_ndft_data = otfs_sym_demap(rx_sym_dd_ndft, idx_pilot_sym, num)

% demap data symbols
data_sym_left = rx_sym_dd_ndft(:, 1:idx_pilot_sym(2));
data_sym_up = rx_sym_dd_ndft(1:idx_pilot_sym(1), idx_pilot_sym(2)+1:idx_pilot_sym(2)+num.len_pilot_otfssym);
data_sym_down = rx_sym_dd_ndft(idx_pilot_sym(1)+num.len_pilot_subc+1:end, idx_pilot_sym(2)+1:idx_pilot_sym(2)+num.len_pilot_otfssym);
data_sym_mid = [data_sym_up; data_sym_down];
data_sym_right = rx_sym_dd_ndft(:, idx_pilot_sym(2)+num.len_pilot_otfssym+1:end);

% serialize data symbols
rx_sym_dd_ndft_data = [data_sym_left(:); data_sym_mid(:); data_sym_right(:)];

% % dump variables
% assignin('base', 'idx_pilot_sym', idx_pilot_sym);
% assignin('base', 'rx_sym_dd_ndft', rx_sym_dd_ndft);
% assignin('base', 'data_sym_left', data_sym_left);
% assignin('base', 'data_sym_up', data_sym_up);
% assignin('base', 'data_sym_down', data_sym_down);
% assignin('base', 'data_sym_mid', data_sym_mid);
% assignin('base', 'data_sym_right', data_sym_right);
% assignin('base', 'rx_sym_dd_ndft_data', rx_sym_dd_ndft_data);
% pause

end
