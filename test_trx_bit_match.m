
% create crc generator and detector objects
tx_crc = comm.CRCGenerator(nw_cc.gCRC24A);
rx_crc = comm.CRCDetector(nw_cc.gCRC24A);

% create turbo encoder/decoder objects
turbo_enc = comm.TurboEncoder('TrellisStructure', nw_cc.tc_trellis, ...
    'InterleaverIndices', nw_cc.PI+1);
turbo_dec = comm.TurboDecoder('TrellisStructure', nw_cc.tc_trellis, ...
    'InterleaverIndices', nw_cc.PI+1, 'NumIterations', nw_cc.num_iter_max);

% evaluate the parameters
carrier_freq_mhz = 30000;
bw_mhz = 200;
scs_khz = 120;
num_slot = 1;
len_tb_bit = 4053;
mcs = 17;
num_sim = 1;
chest_option = 'real';
waveform = 'ofdm';
usr_option = 'mu';
sim_option.override = false;

% set equalization option
nw_num = nw_num_prm_r1(carrier_freq_mhz, scs_khz, bw_mhz, num_slot, waveform, chest_option, usr_option, sim_option);
nw_cc = nw_cc_prm(len_tb_bit);
nw_rm = nw_rm_prm(len_tb_bit, mcs, nw_num, nw_cc, sim_option, test_option);
nw_sim = nw_sim_prm(len_tb_bit, num_sim, nw_num, nw_cc, nw_rm);

% generate bit stream
tx_bit = randi([0 1], nw_sim.len_tb_bit, 1);

% pad bits
tx_bit_pad = [tx_bit; zeros(nw_cc.F, 1)];

% buffer per code block
tx_bit_buff = reshape(tx_bit_pad, (nw_sim.len_tb_bit+nw_cc.F)/nw_cc.C, nw_cc.C);

% initialize
tx_sym = zeros(nw_num.num_qamsym_usr, sum(nw_rm.num_usrfrm_cb));

% simulate per-cb tx process
for idx_cb = 1:nw_cc.C
    % generate crc
    tx_bit_crc_cb = tx_crc(tx_bit_buff(:, idx_cb));
    
    % turbo encode data
    tx_bit_enc_cb = turbo_enc(tx_bit_crc_cb);
%     tx_bit_enc_cb = lteTurboEncode(tx_bit_crc_cb);
    
    % rate match
    tx_bit_ratematch_cb = tx_ratematch_r2(tx_bit_enc_cb, idx_cb, nw_rm);
    
    % modulate bit stream
    % accumulate output per code block serially
    tx_sym_serial_cb = qammod(tx_bit_ratematch_cb(:), 2^nw_rm.Qm, 'InputType', 'bit', 'UnitAveragePower', true);
    
    % simulate per user frame
    tx_sym_cb = reshape(tx_sym_serial_cb, nw_num.num_qamsym_usr, []);
    tx_sym(:, sum(nw_rm.num_usrfrm_cb(1:idx_cb-1))+1:sum(nw_rm.num_usrfrm_cb(1:idx_cb))) = tx_sym_cb;
end

rx_sym = tx_sym;
error_var = 1;

% initialize
rx_bit_crc_removed_cb = zeros(nw_cc.K-nw_cc.L, nw_cc.C);
rx_crc_error_cb = zeros(1, nw_cc.C);

% simulate per-cb rx process
for idx_cb = 1:nw_cc.C
    % serialize qam symbols
    rx_sym_cb = rx_sym(:, sum(nw_rm.num_usrfrm_cb(1:idx_cb-1))+1:sum(nw_rm.num_usrfrm_cb(1:idx_cb)));
    rx_sym_serial_cb = reshape(rx_sym_cb, [], 1);
    
    % serialize noise variance
    error_var_cb = reshape(error_var(:, sum(nw_rm.num_usrfrm_cb(1:idx_cb-1))+1:sum(nw_rm.num_usrfrm_cb(1:idx_cb))), [], 1);
    
    % demap qam symbols
    rx_bit_demod_cb = (-1) * qamdemod(rx_sym_serial_cb.', 2^nw_rm.Qm, 'UnitAveragePower', true, ...
        'OutputType', 'approxllr', 'NoiseVariance', error_var_cb.');
    rx_bit_demod_cb = rx_bit_demod_cb(:);
    
    % rate match
    % output fixed (nw_rm.Kw)
    rx_bit_ratematch_cb = rx_ratematch_r2(rx_bit_demod_cb, nw_rm);
    
    % turbo decode data
    rx_bit_dec_cb = turbo_dec(rx_bit_ratematch_cb);
    
    % detect crc
    [rx_bit_crc_removed_cb(:, idx_cb), rx_crc_error_cb(1, idx_cb)] = rx_crc(rx_bit_dec_cb);
end

% buffer per code block
rx_bit_buff = reshape(rx_bit_crc_removed_cb, [], 1);

% remove padded bits
rx_bit = rx_bit_buff(1:nw_sim.len_tb_bit);

% rearrange crc
rx_crc_error = sum(double(rx_crc_error_cb))>0;



figure, plot(real(tx_sym(:)), '-bo'), hold on, plot(real(rx_sym(:)), ':rx'), hold off, grid minor
figure, plot(real(tx_sym_cb(:)), '-bo'), hold on, plot(real(rx_sym_cb(:)), ':rx'), hold off, grid minor
figure, plot(real(tx_sym_serial_cb(:)), '-bo'), hold on, plot(real(rx_sym_serial_cb(:)), ':rx'), hold off, grid minor
figure, plot(tx_bit_ratematch_cb(:), '-bo'), hold on, plot(rx_bit_demod_cb(:), ':rx'), hold off, grid minor
figure, plot(tx_bit_enc_cb(:), '-bo'), hold on, plot(rx_bit_ratematch_cb(:), ':rx'), hold off, grid minor
figure, plot(tx_bit_crc_cb(:), '-bo'), hold on, plot(rx_bit_dec_cb(:), ':rx'), hold off, grid minor
figure, plot(tx_bit_buff(:), '-bo'), hold on, plot(rx_bit_crc_removed_cb(:), ':rx'), hold off, grid minor
figure, plot(tx_bit(:), '-bo'), hold on, plot(rx_bit(:), ':rx'), hold off, grid minor
