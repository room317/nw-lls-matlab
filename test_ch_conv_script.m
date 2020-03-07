% test_ch_conv script
% created: 2020.03.04
% modified:
%   -

% set test mode
test_pilot = true;
test_signal = false;
test_gauss_ch = false;

% set parameters
NumSim = 1000;
SNRStart = 10;
SNRCnt = 10;
SNRStep = 1;
ErrorStop = 400;
ErrorBreak = 0.001;

% start simulation
test_sim_result = zeros(SNRCnt, 4);
for snr_idx = 1 : SNRCnt
    snr_db = SNRStart + SNRStep * (snr_idx - 1);
    
    % initialize counter
    total_error_tfeq = 0;
    total_error_ddeq = 0;
    sim_cnt = 0;
    
    for sim_idx = 1 : NumSim
        % simulate single packet
        [bit_error_tfeq, bit_error_ddeq] = test_ch_conv_r1(snr_db, test_pilot, test_signal, test_gauss_ch, false);
        
        % count packet error
        total_error_tfeq = total_error_tfeq + double(bit_error_tfeq > 0);
        total_error_ddeq = total_error_ddeq + double(bit_error_ddeq > 0);
        sim_cnt = sim_cnt + 1;
        if min(total_error_tfeq, total_error_ddeq) > ErrorStop
            break
        end
    end
    
    % compute bit error rate
    test_sim_result(snr_idx, 1) = snr_db;                % snr
    test_sim_result(snr_idx, 2) = sim_cnt;                % total number of packets
    test_sim_result(snr_idx, 3) = total_error_tfeq / sim_cnt;  % per of tfeq
    test_sim_result(snr_idx, 4) = total_error_ddeq / sim_cnt;  % per of ddeq
    
    % save simulation result
    fprintf('    SNR: %4.1f     NUM PKTS: %7d   PER TFEQ: %9.6f   PER DDEQ: %9.6f\n', test_sim_result(snr_idx, :));
    
    if max(test_sim_result(snr_idx, 3:4)) < ErrorBreak
        disp('BREAK');
        break;
    end
end

% plot
assignin('base', 'test_sim_result', test_sim_result);
figure
semilogy(test_sim_result(:, 1), test_sim_result(:, 3), '-bo'), hold on
semilogy(test_sim_result(:, 1), test_sim_result(:, 4), ':rx'), hold off
xlabel('SNR (dB)'); ylabel('PER');
title('QAM over Fading Channel');
axis([SNRStart snr_db 1e-4 1])
grid minor
legend('TFEQ', 'DDEQ')
