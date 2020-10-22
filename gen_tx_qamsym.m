% generate ofdm/otfs qam symbols per user
% [tx_bit, tx_sym_serial, tx_sym] = gen_tx_qamsym(sim, cc, rm, num, turbo_enc)
%   - tx_bit, tx_sym_serial, tx_sym] = gen_tx_qamsym(sim, cc, rm, num, turbo_enc)
%   - tx_sym_serial, tx_sym] = gen_tx_qamsym(sim, cc, rm, num, turbo_enc)
%   - tx_sym] = gen_tx_qamsym(sim, cc, rm, num, turbo_enc)
%   - sim: simulation parameters
%   - cc: channel coding parameters
%   - rm: rate matching parameters
%   - num: numerology parameters
%   - turbo_enc: turbo encoder object

function [tx_bit, tx_sym_serial, tx_sym] = gen_tx_qamsym(tx_crc, turbo_enc, sim, cc, rm, num)

% generate bit stream
tx_bit = randi([0 1], sim.len_tb_bit, 1);

% pad bits
tx_bit_pad = [tx_bit; zeros(cc.F, 1)];

% buffer per codeword
tx_bit_buff = reshape(tx_bit_pad, (sim.len_tb_bit+cc.F)/cc.C, cc.C);

% generate crc
tx_bit_crc = tx_crc(tx_bit_buff);

% turbo encode data
tx_bit_enc = zeros(rm.D*3, cc.C);
for idx_cw = 1 : cc.C
    tx_bit_enc(:, idx_cw) = turbo_enc(tx_bit_crc(:, idx_cw));
end

% rate match
tx_bit_ratematch = tx_ratematch_r2(tx_bit_enc, rm);

% modulate bit stream
tx_sym_serial = qammod(tx_bit_ratematch(:), 2^rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);

% simulate per user frame
tx_sym = reshape(tx_sym_serial, num.num_qamsym_usr, []);

end
