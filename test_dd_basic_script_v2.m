% test_dd_basic script
%   - test channel estimation error with respect to snr
% created: 2020.12.03
% modified:
%   -

% set test mode
test_synch = 0;
test_scope = false;
test_seed = -1;
test_partial = false;

% set simulation parameters
NumSim = 1000;
ListSNR = 0:2:20;

% set numerical parameters
scs_khz = 60;
bw_mhz = 5;
num_slot = 4;
test_pilot = 'impulse';
test_pilot_pos = {[0, 0]}; % {[0, 0], [0, 7]};
test_seq_len = 32;

% print simulation result
fprintf('  Test in progress: %6.2f %%', 0);

% start simulation
test_sim_result = zeros(2, length(ListSNR));
for idx_snr = 1:length(ListSNR)
    % set snr
    snr_db = ListSNR(idx_snr);
    
    % initialize counter
    total_ch_mse_dd = 0;
    total_ch_mse_tf = 0;
    
    for idx_sim = 1:NumSim
        % simulate single packet
        [ch_rmse_dd, ch_rmse_tf] = test_dd_basic_r3(snr_db, scs_khz, bw_mhz, num_slot, test_synch, test_pilot, test_pilot_pos, test_seq_len, test_scope, test_seed, test_partial);
        
        % sum channel mse
        total_ch_mse_dd = total_ch_mse_dd+ch_rmse_dd.^2;
        total_ch_mse_tf = total_ch_mse_tf+ch_rmse_tf.^2;
    end
    
    % compute bit error rate
    test_sim_result(1, idx_snr) = sqrt(total_ch_mse_dd/NumSim);
    test_sim_result(2, idx_snr) = sqrt(total_ch_mse_tf/NumSim);
    
    % print simulation result
    fprintf('\b\b\b\b\b\b\b\b');
    fprintf('%6.2f %%', idx_snr/length(ListSNR)*100);
end

% print simulation result
fprintf('\b\b\b\b\b\b\b\b');
fprintf('%6.2f %%\n', 100);

figure
semilogy(ListSNR, test_sim_result(1, :), '-ko'), hold on
semilogy(ListSNR, test_sim_result(2, :), '-kx'), hold off
grid minor, xlabel('SNR (dB)'), ylabel('Channel MSE'), title('Channel MSE'), legend('dd', 'tf')
