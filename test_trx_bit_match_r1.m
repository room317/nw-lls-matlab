% lte toolbox vs. manual coding

% evaluate the parameters
carrier_freq_mhz = 30000;
bw_mhz = 200;
scs_khz = 120;
num_slot = 1;
len_tb_bit = 20000; % 4053;
mcs = 17;
num_sim = 1;
chest_option = 'real';
waveform = 'ofdm';
usr_option = 'mu';
sim_option.override = false;

rv_idx = 0;

test_option.license = true;
test_option.rm_version = 2;

% set equalization option
nw_cc = nw_cc_prm_r1(len_tb_bit);
nw_num = nw_num_prm_r1(carrier_freq_mhz, scs_khz, bw_mhz, num_slot, waveform, chest_option, usr_option, sim_option);
nw_rm = nw_rm_prm_r1(mcs, nw_cc, nw_num, sim_option);

% create crc generator and detector objects
tx_crc_a = comm.CRCGenerator(nw_cc.gCRC24A);
tx_crc_b = comm.CRCGenerator(nw_cc.gCRC24B);
rx_crc_a = comm.CRCDetector(nw_cc.gCRC24A);
rx_crc_b = comm.CRCDetector(nw_cc.gCRC24B);

% % create turbo encoder/decoder objects
% turbo_enc = comm.TurboEncoder('TrellisStructure', nw_cc.tc_trellis, ...
%     'InterleaverIndices', nw_cc.PI+1);
% turbo_dec = comm.TurboDecoder('TrellisStructure', nw_cc.tc_trellis, ...
%     'InterleaverIndices', nw_cc.PI+1, 'NumIterations', nw_cc.num_iter_max);


% tx--------------

% generate bit stream
tx_bit = randi([0 1], nw_cc.A, 1);

% channel coding, multiplexing and interleaving
if test_option.license
    % crc calculation
    tx_bit_crc = lteCRCEncode(tx_bit, '24A');
    
    % code block segmentation and code block crc attachment
    tx_bit_seg = lteCodeBlockSegment(tx_bit_crc);
    
    % channel coding
    tx_bit_enc = lteTurboEncode(tx_bit_seg);
    
    % rate matching and code block concatenation
    tx_bit_ratematch = lteRateMatchTurbo(tx_bit_enc, nw_rm.G, rv_idx);
    
    % symbol modulate the codeword
    switch nw_rm.Q_m
        case 2
            qam_mod = 'QPSK';
        case 4
            qam_mod = '16QAM';
        case 6
            qam_mod = '64QAM';
        otherwise
            error('Wrong modulation order.')
    end
    tx_sym_serial = lteSymbolModulate(tx_bit_ratematch, qam_mod);
else
    % crc calculation
    tx_bit_crc = tx_crc_a(tx_bit);
    
    % code block segmentation and code block crc attachment
    tx_bit_seg = cell(1, nw_cc.C);
    if nw_cc.C > 1
        for r = 1:nw_cc.C
            if r == 1
                % insert filler bits and calculate crc
                tx_bit_seg{r} = tx_crc_b([zeros(nw_cc.F, 1); tx_bit_crc(1:nw_cc.K_r(r)-nw_cc.L_-nw_cc.F)]);
            else
                % calculate crc
                tx_bit_seg{r} = tx_crc_b(tx_bit_crc(sum(nw_cc.K_r(1:r-1))-(r-1)*nw_cc.L_-nw_cc.F+1:sum(nw_cc.K_r(1:r))-r*nw_cc.L_-nw_cc.F));
            end
        end
    else
        tx_bit_seg{1} = [zeros(nw_cc.F, 1); tx_bit_crc(1:nw_cc.K_r(1)-nw_cc.L_-nw_cc.F)];
    end
    
    % channel coding
    tx_bit_enc = cell(1, nw_cc.C);
    for r = 1:nw_cc.C
        % create turbo encoder/decoder objects
        turbo_enc = comm.TurboEncoder('TrellisStructure', nw_cc.tc_trellis, ...
            'InterleaverIndices', nw_cc.PI_r{r});
        
        % turbo encode data
        tx_bit_enc{r} = turbo_enc(tx_bit_seg{r});
    end
    
    % rate matching and code block concatenation
    tx_bit_ratematch = tx_ratematch_r3(tx_bit_enc, nw_cc, nw_rm, rv_idx, test_option);
    
    % modulate bit stream
    switch nw_rm.Q_m
        case 2
            lteSymMap = 'gray';
        case 4
            lteSymMap = [11 10 14 15 ... 
                          9  8 12 13 ...
                          1  0  4  5 ...
                          3  2  6  7];
        case 6
            lteSymMap = 'gray';
        otherwise
            error('Wrong modulation order.')
    end
    tx_sym_serial = qammod(tx_bit_ratematch(:), 2^nw_rm.Q_m, lteSymMap, 'InputType', 'bit', 'UnitAveragePower', true);
end

% simulate per resource block
tx_sym = reshape(tx_sym_serial, nw_num.num_qamsym_usr, []);


% ch--------------

rx_sym = tx_sym;
error_var = 1;


% rx--------------

% serialize rx symbols
rx_sym_serial = reshape(rx_sym, [], 1);

% channel coding, multiplexing and interleaving
if test_option.license
    % symbol modulate the codeword
    rx_bit_demod = lteSymbolDemodulate(rx_sym_serial, qam_mod, 'Soft');
    
    % rate recovery
    rx_bit_raterecover = lteRateRecoverTurbo(rx_bit_demod, len_tb_bit, rv_idx);
    
    % turbo decoder
    rx_bit_dec = lteTurboDecode(rx_bit_raterecover, nw_cc.num_iter_max);
    
    % code block segmentation and code block crc attachment
    [rx_bit_cbrecover, segError] = lteCodeBlockDesegment(rx_bit_dec, len_tb_bit+24);
    
    % crc decoder
    [rx_bit, crcError] = lteCRCDecode(rx_bit_cbrecover, '24A');
else
    % demap qam symbols
    rx_bit_demod = (-1) * qamdemod(rx_sym_serial.', 2^nw_rm.Q_m, lteSymMap, 'UnitAveragePower', true, ...
        'OutputType', 'approxllr', 'NoiseVariance', error_var.');
    rx_bit_demod = rx_bit_demod(:);
    
    % rate recovery
    rx_bit_raterecover = lteRateRecoverTurbo(rx_bit_demod, len_tb_bit, rv_idx);
    
    % turbo decoder
    rx_bit_dec = lteTurboDecode(rx_bit_raterecover, nw_cc.num_iter_max);
    
    % code block segmentation and code block crc attachment
    [rx_bit_cbrecover, segError] = lteCodeBlockDesegment(rx_bit_dec, len_tb_bit+24);
    
    % crc decoder
    [rx_bit, crcError] = lteCRCDecode(rx_bit_cbrecover, '24A');
end
