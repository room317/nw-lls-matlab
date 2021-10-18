% gen_rx_bit generates ofdm/otfs rx bits per user. (cell-base signal processing)
% [rx_bit, rx_crc_error] = gen_rx_bit_usr_r2(rx_sym, error_var, rv_idx, cc, rm, test_option)
%   - rx_bit: received bits
%   - rx_crc_error: rx crc error
%   - rx_sym: received qam symbols per user
%   - error_var: qam symbol error variance for llr calculation
%   - rv_idx: redundancy version
%   - cc, rm, num: simulation parameters
%   - test_option: license check

function [rx_bit, rx_crc_error] = gen_rx_bit_usr_r2(rx_sym, error_var, rv_idx, rx_crc_a, cc, rm, test_option)

% serialize qam symbols
rx_sym_serial = reshape(rx_sym, [], 1);

% channel coding, multiplexing and interleaving
if test_option.license
    if test_option.qam_modem_toolbox
        % symbol modulate the codeword
        rx_bit_demod = lteSymbolDemodulate(rx_sym_serial, rm.Q_mod, 'Soft');
    else
        % serialize noise variance
        error_var_serial = reshape(error_var, [], 1);
        
        % demap qam symbols
        rx_bit_demod = (-1) * qamdemod(rx_sym_serial.', 2^rm.Q_m, rm.lteSymMap, 'UnitAveragePower', true, ...
            'OutputType', 'approxllr', 'NoiseVariance', error_var_serial.');
        rx_bit_demod = rx_bit_demod(:);
    end
    
    % rate recovery
    rx_bit_raterecover = lteRateRecoverTurbo(rx_bit_demod, cc.A, rv_idx);
    
    % turbo decoder
    rx_bit_dec = lteTurboDecode(rx_bit_raterecover, cc.num_iter_max);
    
    % code block segmentation and code block crc attachment
    [rx_bit_cbrecover, ~] = lteCodeBlockDesegment(rx_bit_dec, cc.A+24);
    
    % crc decoder
    [rx_bit, rx_crc_error] = lteCRCDecode(rx_bit_cbrecover, '24A');
else
    % serialize noise variance
    error_var_serial = reshape(error_var, [], 1);
    
    % demap qam symbols
    rx_bit_demod = (-1) * qamdemod(rx_sym_serial.', 2^rm.Q_m, rm.lteSymMap, 'UnitAveragePower', true, ...
        'OutputType', 'approxllr', 'NoiseVariance', error_var_serial.');
    rx_bit_demod = rx_bit_demod(:);
    
    % rate recovery
    rx_bit_raterecover = lteRateRecoverTurbo(rx_bit_demod, cc.A, rv_idx);
    
    % turbo decoder
    rx_bit_dec = lteTurboDecode(rx_bit_raterecover, cc.num_iter_max);
    
    % code block segmentation and code block crc attachment
    [rx_bit_cbrecover, ~] = lteCodeBlockDesegment(rx_bit_dec, cc.A+24);
    
    % detect crc
    [rx_bit, rx_crc_error] = rx_crc_a(double(rx_bit_cbrecover));
end

% rearrange crc
rx_crc_error = double(rx_crc_error)>0;

% % dump
% assignin('base', 'rx_sym_serial', rx_sym_serial);
% assignin('base', 'rx_bit_demod', rx_bit_demod);
% assignin('base', 'rx_bit_raterecover', rx_bit_raterecover);
% assignin('base', 'rx_bit_dec', rx_bit_dec);
% assignin('base', 'rx_bit_cbrecover', rx_bit_cbrecover);
% assignin('base', 'rx_bit', rx_bit);
% assignin('base', 'rx_crc_error', rx_crc_error);

end
