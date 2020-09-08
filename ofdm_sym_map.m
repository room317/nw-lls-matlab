% map qam symbols to user resource blocks
%  - tx_sym_data_usrfrm: user frame data per user
%  - tx_sym_pilot_usrfrm: for error variance calculation

function [tx_sym_rbs, tx_sym_pilot_usrfrm] = ofdm_sym_map(tx_sym_data_usrfrm, num)

% generate pilot symbols (random sequence for now)
tx_sym_pilot_usrfrm = 1/sqrt(2)* ...
    complex(randi([0 1], num.num_subc_pilot_usr, num.num_ofdmsym_pilot_usr)*2-1, ...
    randi([0 1], num.num_subc_pilot_usr, num.num_ofdmsym_pilot_usr)*2-1);

% reshape data
tx_sym_data1_usrfrm = reshape(tx_sym_data_usrfrm(1:num.num_subc_usr*num.num_ofdmsym_data_usr), num.num_subc_usr, num.num_ofdmsym_data_usr);
tx_sym_data2_usrfrm = reshape(tx_sym_data_usrfrm(num.num_subc_usr*num.num_ofdmsym_data_usr+1:end), num.num_subc_data_usr, num.num_ofdmsym_pilot_usr);

% map data1
tx_sym_rbs = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
tx_sym_rbs(:, num.idx_ofdmsym_data_usr) = tx_sym_data1_usrfrm;

% map data2 and pilot to user resource block
for i = 1:num.num_ofdmsym_pilot_usr
    tx_sym_rbs(num.idx_subc_data_usr(:, i), num.idx_ofdmsym_pilot_usr(i)) = tx_sym_data2_usrfrm(:, i);
    tx_sym_rbs(num.idx_subc_pilot_usr(:, i), num.idx_ofdmsym_pilot_usr(i)) = tx_sym_pilot_usrfrm(:, i);
end

end
