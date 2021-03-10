% single ofdm symbol transmitter
% [tx_bit, tx_sym, tx_sym_rbs] = gen_tx_qamsym_usr_r1(sim, cc, rm, num, tx_crc, turbo_enc)
%   - sim, cc, rm, num: simulation parameters
%   - tx_crc: crc generator
%   - turbo_enc: turbo encoder
%   - tx_bit: tx bits for a single user
%   - tx_sym: tx qam symbols per code block for a single user
%   - tx_sym_rbs: tx qam symbols per rbs for a single user
%   - tx_sym_bw: tx qam symbols per full bandwidth for a single user
% created: 2019.10.15
% modified:
%   - channel coding added: 2019.10.21
%   - rate matching added: 2019.10.22
%   - structure updated: 2019.10.23
%   - rate matching bug fixed: 2019.11.08
%   - slot buffer fixed: 2019.12.10
%   - real channel estimation with single tone pilot added: 2020.02.09
%   - slot-based to user-frame-based simulation (user frame = multiple slots)

function [tx_bit, tx_sym] = gen_tx_qamsym_usr_r1(sim, cc, rm, num, tx_crc, turbo_enc)

% generate bit stream
tx_bit = randi([0 1], sim.len_tb_bit, 1);

% pad bits
tx_bit_pad = [tx_bit; zeros(cc.F, 1)];

% buffer per code block
tx_bit_buff = reshape(tx_bit_pad, (sim.len_tb_bit+cc.F)/cc.C, cc.C);

% initialize
tx_sym = zeros(num.num_qamsym_usr, sum(rm.num_usrfrm_cb));

% simulate per-cb tx process
for idx_cb = 1:cc.C
    % generate crc
    tx_bit_crc_cb = tx_crc(tx_bit_buff(:, idx_cb));
    
    % turbo encode data
    tx_bit_enc_cb = turbo_enc(tx_bit_crc_cb);
    
    % rate match
    tx_bit_ratematch_cb = tx_ratematch_r2(tx_bit_enc_cb, idx_cb, rm);
    
%     % temporary for otfs test
%     tx_bit_enc_cb = repmat([tx_bit_crc_cb; tx_bit_crc_cb(1:4)], 3, 1);
%     tx_bit_ratematch_cb = tx_bit_enc_cb(1:rm.E_true(idx_cb), :);
    
    % modulate bit stream
    % accumulate output per code block serially
    tx_sym_serial_cb = qammod(tx_bit_ratematch_cb(:), 2^rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);
    
    % simulate per user frame
    tx_sym_cb = reshape(tx_sym_serial_cb, num.num_qamsym_usr, []);
    tx_sym(:, sum(rm.num_usrfrm_cb(1:idx_cb-1))+1:sum(rm.num_usrfrm_cb(1:idx_cb))) = tx_sym_cb;
    
%     if idx_cb == 1
%         assignin('base', 'tx_bit_crc_cb', tx_bit_crc_cb)
%         assignin('base', 'tx_bit_enc_cb', tx_bit_enc_cb)
%         assignin('base', 'tx_bit_ratematch_cb', tx_bit_ratematch_cb)
%         assignin('base', 'tx_sym_serial_cb', tx_sym_serial_cb)
%         assignin('base', 'tx_sym_cb', tx_sym_cb)
%     end
end

% % initialize
% tx_sym_rbs = zeros(num.num_subc_usr, num.num_ofdmsym_usr, sum(rm.num_usrfrm_cb));
% 
% % simulate per-user-frame tx process
% for idx_usrfrm = 1:sum(rm.num_usrfrm_cb)
%     % extract user frame data per user
%     tx_sym_data_usrfrm = squeeze(tx_sym(:, idx_usrfrm));
%     
%     % map qam symbols to user physical resource blocks
%     [tx_sym_rbs_usrfrm, ~] = ofdm_sym_map(tx_sym_data_usrfrm, num);
%     tx_sym_rbs(:, :, idx_usrfrm) = tx_sym_rbs_usrfrm;
% end

% % dump
% assignin('base', 'tx_bit', tx_bit);
% assignin('base', 'tx_sym', tx_sym);
% assignin('base', 'tx_sym_rbs', tx_sym_rbs);

end
