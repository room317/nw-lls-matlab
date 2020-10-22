% gen_rx_bit generates ofdm/otfs rx bits per user.
% [rx_bit, rx_sym_serial, rx_crc_error] = gen_rx_bit(rx_sym, turbo_dec, error_var, sim, cc, rm)
%   - rx_bit: received bits
%   - rx_sym_serial: serialized received qam symbols
%   - rx_crc_error: rx crc error
%   - rx_sym: received qam symbols per user
%   - turbo_dec: turbo decoder object
%   - error_var: qam symbol error variance for llr calculation
%   - sim: simulation parameters
%   - cc: channel coding parameters
%   - rm: rate matching parameters

function [rx_bit, rx_sym_serial, rx_crc_error] = gen_rx_bit(rx_sym, rx_crc, turbo_dec, error_var, sim, cc, rm)

% serialize qam symbols
rx_sym_serial = rx_sym(:);

% demap the qam symbols
% noise_var_tmp = var(tx_sym_serial-rx_sym_serial);
rx_bit_demod = (-1) * qamdemod(rx_sym_serial, 2^rm.Qm, 'UnitAveragePower', true, ...
    'OutputType', 'approxllr', 'NoiseVariance', error_var);

% buffer per codeword
rx_bit_buff = reshape(rx_bit_demod, [], cc.C);  % rm.E*cc.C

% rate match
% rx_bit_ratematch = rx_ratematch(rx_bit_demod, rm);
rx_bit_ratematch = rx_ratematch_r2(rx_bit_buff, rm);

% turbo decode data
rx_bit_dec = zeros(cc.K, cc.C);
for idx_cw = 1:cc.C
    rx_bit_dec(:, idx_cw) = turbo_dec(rx_bit_ratematch(:, idx_cw));
end

% detect crc
[rx_bit_crc_removed, rx_crc_error] = rx_crc(rx_bit_dec);

% remove padded bits
rx_bit = rx_bit_crc_removed(1 : sim.len_tb_bit);

end
