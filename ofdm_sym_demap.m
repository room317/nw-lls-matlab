% demap data qam symbols from user resource blocks
%  - tx_sym_data_usrfrm: user frame data per user (vector)
%  - tx_sym_pilot_usrfrm: user frame pilot per user (matrix)

function [rx_sym_data_usrfrm, rx_sym_pilot_usrfrm] = ofdm_sym_demap(rx_sym_rbs, num)

% demap data1 from user resource block
rx_sym_data1_usrfrm = rx_sym_rbs(:, num.idx_ofdmsym_data_usr);

% demap data2 from user resource block
rx_sym_data2_usrfrm = zeros(num.num_subc_data_usr, num.num_ofdmsym_pilot_usr);
for i = 1:num.num_ofdmsym_pilot_usr
    rx_sym_data2_usrfrm(:, i) = rx_sym_rbs(num.idx_subc_data_usr(:, i), num.idx_ofdmsym_pilot_usr(i));
end

rx_sym_data_usrfrm = [rx_sym_data1_usrfrm(:); rx_sym_data2_usrfrm(:)];

% demap pilot from user resource block
rx_sym_pilot_usrfrm = zeros(num.num_subc_pilot_usr, num.num_ofdmsym_pilot_usr);
for i = 1:num.num_ofdmsym_pilot_usr
    rx_sym_pilot_usrfrm(:, i) =  rx_sym_rbs(num.idx_subc_pilot_usr(:, i), num.idx_ofdmsym_pilot_usr(i));
end

end
