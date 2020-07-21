% test_dd_basic script
%   - test channel estimation error with respect to the position of the
%     pilot tone
% created: 2020.07.19
% modified:
%   -

% set test mode
test_synch = 0;
test_scope = false;
test_seed = -1;

% set parameters
NumSim = 100;
NumSubc = 600;
NumSym = 14;
    
% print simulation result
fprintf('  Test in progress: %6.2f %%', 0);

% start simulation
test_sim_result = zeros(1, NumSym);
for idx_num_rx_sym = 1 : NumSym
    % initialize counter
    total_ch_est_mse = 0;
    
    for idx_sim = 1 : NumSim
        % simulate single packet
%             ch_est_rmse = test_dd_basic(test_synch, [idx_subc, idx_sym], test_scope, test_seed);
%             ch_est_rmse = test_dd_basic_r1(test_synch, [idx_subc, idx_sym], test_scope, test_seed, 1);
        [ch_est_rmse, ~] = test_dd_basic_r2(test_synch, [NumSubc/2, NumSym/2], test_scope, test_seed, 1, idx_num_rx_sym);
        
        % sum channel mse
        total_ch_est_mse = total_ch_est_mse + ch_est_rmse.^2;
    end
    
    % compute bit error rate
    test_sim_result(1, idx_num_rx_sym) = sqrt(total_ch_est_mse/NumSim);
    
    % print simulation result
    fprintf('\b\b\b\b\b\b\b\b');
    fprintf('%6.2f %%', idx_num_rx_sym/NumSym*100);
end

% print simulation result
fprintf('\b\b\b\b\b\b\b\b');
fprintf('%6.2f %%\n', 100);

figure
semilogy(1:NumSym, test_sim_result)
xlabel('Number of symbols received'), ylabel('Channel MSE'), title('Channel MSE')
