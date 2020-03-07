% test_ch_conv script
% created: 2020.03.06
% modified:
%   -

% set test mode
test_pilot = true;
test_cp_position = 'postfix';
test_synch = 7;

% set parameters
NumSim = 10000;
SNRStart = 20;
SNRCnt = 10;
SNRStep = 2;
ErrorStop = 400;
ErrorBreak = 0.001;

% start simulation
test_sim_result = zeros(SNRCnt, 6);
for snr_idx = 1 : SNRCnt
    snr_db = SNRStart + SNRStep * (snr_idx - 1);
    
    % initialize counter
    total_qam_error_tfeq = 0;
    total_qam_error_ddeq = 0;
    total_pkt_error_tfeq = 0;
    total_pkt_error_ddeq = 0;
    sim_cnt = 0;
    
    for sim_idx = 1 : NumSim
        % simulate single packet
        [qam_error_tfeq, qam_error_ddeq, num_qam_per_pkt] = test_dd_pilot(snr_db, test_pilot, test_cp_position, test_synch, false);
        
        % count packet error
        total_qam_error_tfeq = total_qam_error_tfeq + qam_error_tfeq;
        total_qam_error_ddeq = total_qam_error_ddeq + qam_error_ddeq;
        total_pkt_error_tfeq = total_pkt_error_tfeq + double(qam_error_tfeq > 0);
        total_pkt_error_ddeq = total_pkt_error_ddeq + double(qam_error_ddeq > 0);
        sim_cnt = sim_cnt + 1;
        if min(total_pkt_error_tfeq, total_pkt_error_ddeq) > ErrorStop
            break
        end
    end
    
    % compute bit error rate
    test_sim_result(snr_idx, 1) = snr_db;                % snr
    test_sim_result(snr_idx, 2) = sim_cnt;                % total number of packets
    test_sim_result(snr_idx, 3) = total_qam_error_tfeq / sim_cnt / num_qam_per_pkt;  % per of tfeq
    test_sim_result(snr_idx, 4) = total_qam_error_ddeq / sim_cnt / num_qam_per_pkt;  % per of ddeq
    test_sim_result(snr_idx, 5) = total_pkt_error_tfeq / sim_cnt;  % per of tfeq
    test_sim_result(snr_idx, 6) = total_pkt_error_ddeq / sim_cnt;  % per of ddeq
    
    % save simulation result
    fprintf('  SNR: %4.1f | NUM PKTS: %5d | SER TFEQ: %6.4f | SER DDEQ: %6.4f | PER TFEQ: %6.4f | PER DDEQ: %6.4f\n', test_sim_result(snr_idx, :));
    
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
xlabel('SNR (dB)'); ylabel('SER');
title('QAM over Fading Channel');
axis([SNRStart snr_db 1e-3 1])
grid minor
legend('TFEQ', 'DDEQ')

figure
semilogy(test_sim_result(:, 1), test_sim_result(:, 5), '-bo'), hold on
semilogy(test_sim_result(:, 1), test_sim_result(:, 6), ':rx'), hold off
xlabel('SNR (dB)'); ylabel('PER');
title('QAM over Fading Channel');
axis([SNRStart snr_db 1e-3 1])
grid minor
legend('TFEQ', 'DDEQ')

