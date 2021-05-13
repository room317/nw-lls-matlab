% test_dd_pilot script
%   - v0: test qam symbol error rate and packet error rate with respect to snr
% created: 2020.03.06
% modified:
%   -

% set test mode
test_scope = false;
test_seed = -1;
test_synch = 0;
test_chest = 'dd_golay';
test_cheq = 'tf_mmse';
test_ch_rmse = false;
test_papr = false;
test_usr = 'su';

% set numerical parameters
scs_khz = 60;
bw_mhz = 10;
num_slot = 4;

% set sim parameters
NumSim = 10;
SNRStart = 10;
SNRCnt = 10;
SNRStep = 2;
ErrorStop = 400;
ErrorBreak = 0.001;
% NumSim = 10000;
% METRICStart = 0;
% METROCCnt = 20;
% METRICStep = 25;
% ErrorStop = 400;
% ErrorBreak = 0.001;

% start simulation
test_sim_result = zeros(SNRCnt, 5);
for snr_idx = 1 : SNRCnt
    test_snr = SNRStart + SNRStep * (snr_idx - 1);
    
    % initialize counter
    total_qam_error = 0;
    total_pkt_error = 0;
%     total_ch_est_mse = 0;
    sim_cnt = 0;
    
    for sim_idx = 1 : NumSim
        % simulate single packet
        [qam_error, num_qam_per_pkt, ~, ~] = test_dd_pilot_r6( ...
            test_snr, scs_khz, bw_mhz, num_slot, ...
            test_synch, test_scope, test_seed, test_ch_rmse, test_papr, ...
            test_chest, test_cheq, test_usr);
        
        % count packet error
        total_qam_error = total_qam_error + qam_error;
        total_pkt_error = total_pkt_error + double(qam_error > 0);
%         total_ch_est_mse = total_ch_est_mse + ch_est_rmse.^2;
        sim_cnt = sim_cnt + 1;
        if total_pkt_error > ErrorStop
            break
        end
    end
    
    % compute bit error rate
    test_sim_result(snr_idx, 1) = test_snr;                % snr
    test_sim_result(snr_idx, 2) = sim_cnt;                % total number of packets
    test_sim_result(snr_idx, 3) = total_qam_error / sim_cnt / num_qam_per_pkt;  % ser
    test_sim_result(snr_idx, 4) = total_pkt_error / sim_cnt;  % per
%     test_sim_result(metric_idx, 5) = sqrt(total_ch_est_mse / sim_cnt);  % channel rmse
    
    % save simulation result
    fprintf('  SNR: %4.1f | NUM PKTS: %5d | SER: %6.4f | PER: %6.4f | CH RMSE: %6.4f\n', test_sim_result(snr_idx, :));
    
    if test_sim_result(snr_idx, 3) < ErrorBreak
        disp('BREAK');
        break;
    end
end

% plot
assignin('base', 'test_sim_result', test_sim_result);

figure
semilogy(test_sim_result(:, 1), test_sim_result(:, 3), '-bo'), hold on
semilogy(test_sim_result(:, 1), test_sim_result(:, 4), ':rx'), hold off
xlabel('SNR'), ylabel('SER/PER'), title('QAM over Fading Channel')
axis([SNRStart test_snr 1e-2 1])
grid minor
legend('SER', 'PER')

% figure
% semilogy(test_sim_result(:, 1), test_sim_result(:, 5), '-bo')
% xlabel('METRIC'), ylabel('CH. RMSE'), title('QAM over Fading Channel')
% grid minor
