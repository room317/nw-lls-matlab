% gen_tx_qamsym generates ofdm/otfs tx qam symbols per user. (cell-base signal processing)
% [tx_bit, tx_sym] = gen_tx_qamsym_usr_r2(rv_idx, cc, rm, num, tx_crc_a, tx_crc_b, test_option)
%   - tx_bit: tx bits for a single user
%   - tx_sym: tx qam symbols per code block for a single user
%   - rv_idx: redundancy version
%   - cc, rm, num: simulation parameters
%   - tx_crc_a, tx_crc_b: crc generator
%   - test_option: license check
% created: 2019.10.15
% modified:
%   - channel coding added: 2019.10.21
%   - rate matching added: 2019.10.22
%   - structure updated: 2019.10.23
%   - rate matching bug fixed: 2019.11.08
%   - slot buffer fixed: 2019.12.10
%   - real channel estimation with single tone pilot added: 2020.02.09
%   - slot-based to user-frame-based simulation (user frame = multiple slots)
%   - cell-based processing: 2021.10.06 (toolbox cross-checked)

function [tx_bit, tx_sym] = gen_tx_qamsym_usr_r2(rv_idx, cc, rm, num, tx_crc_a, tx_crc_b, test_option)

% generate bit stream
tx_bit = randi([0 1], cc.A, 1);

% channel coding, multiplexing and interleaving
if test_option.license
    % crc calculation
    tx_bit_crc = lteCRCEncode(tx_bit, '24A');
    
    % code block segmentation and code block crc attachment
    tx_bit_seg = lteCodeBlockSegment(tx_bit_crc);
    
    % channel coding
    tx_bit_enc = lteTurboEncode(tx_bit_seg);
    
    % rate matching and code block concatenation
    tx_bit_ratematch = lteRateMatchTurbo(tx_bit_enc, rm.G, rv_idx);
    
    if test_option.qam_modem_toolbox
        % symbol modulate the codeword
        tx_sym_serial = lteSymbolModulate(tx_bit_ratematch, rm.Q_mod);
    else
        % modulate bit stream
        tx_sym_serial = qammod(tx_bit_ratematch(:), 2^rm.Q_m, rm.lteSymMap, 'InputType', 'bit', 'UnitAveragePower', true);
    end
else
    % crc calculation
    tx_bit_crc = tx_crc_a(tx_bit);
    
    % code block segmentation and code block crc attachment
    tx_bit_seg = cell(1, cc.C);
    if cc.C > 1
        for r = 1:cc.C
            if r == 1
                % insert filler bits and calculate crc
                tx_bit_seg{r} = tx_crc_b([zeros(cc.F, 1); tx_bit_crc(1:cc.K_r(r)-cc.L_-cc.F)]);
            else
                % calculate crc
                tx_bit_seg{r} = tx_crc_b(tx_bit_crc(sum(cc.K_r(1:r-1))-(r-1)*cc.L_-cc.F+1:sum(cc.K_r(1:r))-r*cc.L_-cc.F));
            end
        end
    else
        tx_bit_seg{1} = [zeros(cc.F, 1); tx_bit_crc(1:cc.K_r(1)-cc.L_-cc.F)];
    end
    
    % channel coding
    tx_bit_enc = cell(1, cc.C);
    for r = 1:cc.C
        % create turbo encoder/decoder objects
        turbo_enc = comm.TurboEncoder('TrellisStructure', cc.tc_trellis, ...
            'InterleaverIndices', cc.PI_r{r});
    
        % turbo encode data
        tx_bit_enc{r} = turbo_enc(tx_bit_seg{r});
    end
    
    % rate matching and code block concatenation
    tx_bit_ratematch = tx_ratematch_r3(tx_bit_enc, cc, rm, rv_idx, test_option);

    % modulate bit stream
    tx_sym_serial = qammod(tx_bit_ratematch(:), 2^rm.Q_m, rm.lteSymMap, 'InputType', 'bit', 'UnitAveragePower', true);
end

% simulate per user frame
tx_sym = reshape(tx_sym_serial, num.num_qamsym_usr, []);

% % dump
% assignin('base', 'tx_bit', tx_bit);
% assignin('base', 'tx_bit_crc', tx_bit_crc);
% assignin('base', 'tx_bit_seg', tx_bit_seg);
% assignin('base', 'tx_bit_enc', tx_bit_enc);
% assignin('base', 'tx_bit_ratematch', tx_bit_ratematch);
% assignin('base', 'tx_sym_serial', tx_sym_serial);
% assignin('base', 'tx_sym', tx_sym);

end
