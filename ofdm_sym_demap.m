% demap data qam symbols from user resource blocks
%  - tx_sym_data_slot: slot data per user

function rx_sym_data_slot = ofdm_sym_demap(rx_sym_rbs, num)

% demap data1 from user resource block
rx_sym_data1_slot = rx_sym_rbs(:, num.idx_ofdmsym_data_usr);

% demap data2 from user resource block
rx_sym_data2_slot = zeros(num.num_subc_data_usr, num.num_ofdmsym_pilot_usr);
for i = 1:num.num_ofdmsym_pilot_usr
    rx_sym_data2_slot(:, i) = rx_sym_rbs(num.idx_subc_data_usr(:, i), num.idx_ofdmsym_pilot_usr(i));
end

rx_sym_data_slot = [rx_sym_data1_slot(:); rx_sym_data2_slot(:)];

end
