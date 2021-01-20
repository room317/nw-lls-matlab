% gen_rx_bit generates ofdm/otfs rx bits per user.
% [rx_bit, rx_sym_serial, rx_crc_error] = gen_rx_bit_usr_r1(rx_sym, turbo_dec, error_var, sim, cc, rm)
%   - rx_bit: received bits
%   - rx_crc_error: rx crc error
%   - rx_sym: received qam symbols per user
%   - error_var: qam symbol error variance for llr calculation
%   - rx_crc: crc detector object
%   - turbo_dec: turbo decoder object
%   - sim: simulation parameters
%   - cc: channel coding parameters
%   - rm: rate matching parameters

function [rx_bit, rx_crc_error] = gen_rx_bit_usr_r1(rx_sym, error_var, rx_crc, turbo_dec, sim, cc, rm)

% initialize
rx_bit_crc_removed_cb = zeros(cc.K-cc.L, cc.C);
rx_crc_error_cb = zeros(1, cc.C);

% simulate per-cb rx process
for idx_cb = 1:cc.C
    % serialize qam symbols
    rx_sym_cb = rx_sym(:, sum(rm.num_usrfrm_cb(1:idx_cb-1))+1:sum(rm.num_usrfrm_cb(1:idx_cb)));
    rx_sym_serial_cb = reshape(rx_sym_cb, [], 1);
    
    % serialize noise variance
    error_var_cb = reshape(error_var(:, sum(rm.num_usrfrm_cb(1:idx_cb-1))+1:sum(rm.num_usrfrm_cb(1:idx_cb))), [], 1);
%     error_var_cb = abs(rx_sym_serial_cb-tx_sym_serial).^2;
    
    % demap qam symbols
    % noise_var_tmp = var(tx_sym_serial-rx_sym_serial);
%     rx_bit_demod_cb = (-1) * qamdemod(rx_sym_serial_cb, 2^rm.Qm, 'UnitAveragePower', true, ...
%         'OutputType', 'approxllr', 'NoiseVariance', error_var_cb);
    
    % test
    rx_bit_demod_cb = (-1) * qamdemod(rx_sym_serial_cb.', 2^rm.Qm, 'UnitAveragePower', true, ...
        'OutputType', 'approxllr', 'NoiseVariance', error_var_cb.');
    rx_bit_demod_cb = rx_bit_demod_cb(:);
    
    % rate match
    % output fixed (rm.Kw)
    rx_bit_ratematch_cb = rx_ratematch_r2(rx_bit_demod_cb, rm);
    
    % turbo decode data
    rx_bit_dec_cb = turbo_dec(rx_bit_ratematch_cb);
    
    % detect crc
    [rx_bit_crc_removed_cb(:, idx_cb), rx_crc_error_cb(1, idx_cb)] = rx_crc(rx_bit_dec_cb);
end

% buffer per code block
rx_bit_buff = reshape(rx_bit_crc_removed_cb, [], 1);

% remove padded bits
rx_bit = rx_bit_buff(1:sim.len_tb_bit);

% rearrange crc
rx_crc_error = sum(double(rx_crc_error_cb))>0;

end
