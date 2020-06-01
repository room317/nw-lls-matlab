% test_dd_basic script
% created: 2020.05.01
% modified:
%   -

% set test mode
test_synch = 0;
test_scope = false;
test_seed = -1;

% set parameters
NumSim = 20;
NumSubc = 600;
NumSym = 14;
    
% print simulation result
fprintf('  Resource covered: %6.2f %%', 0);

% start simulation
test_sim_result = zeros(NumSubc, NumSym);
for idx_subc = 0 : NumSubc-1
    for idx_sym = 0 : NumSym-1
        % initialize counter
        total_ch_est_mse = 0;
        
        for idx_sim = 1 : NumSim
            % simulate single packet
%             ch_est_rmse = test_dd_basic(test_synch, [idx_subc, idx_sym], test_scope, test_seed);
            ch_est_rmse = test_dd_basic_r1(test_synch, [idx_subc, idx_sym], test_scope, test_seed, 1);

            % sum channel mse
            total_ch_est_mse = total_ch_est_mse + ch_est_rmse.^2;
        end
        
        % compute bit error rate
        test_sim_result(idx_subc+1, idx_sym+1) = sqrt(mean(total_ch_est_mse));
    
        % print simulation result
        fprintf('\b\b\b\b\b\b\b\b');
        fprintf('%6.2f %%', (NumSym*idx_subc+idx_sym)/(NumSubc*NumSym)*100);
    end
end

% print simulation result
fprintf('\b\b\b\b\b\b\b\b');
fprintf('%6.2f %%\n', 100);

% plot
assignin('base', 'test_sim_result', test_sim_result);

figure
mesh(0:NumSym-1, 0:NumSubc-1, test_sim_result)
xlabel('Doppler'), ylabel('delay'), title('Channel MSE on Resource Plain')
