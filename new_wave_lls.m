% new wave simulation script
% created: 2019.10.15
% modified:
% - doppler sim added: 2019.10.16
% - rb granularity added: 2019.10.16
% - channel estimation added: 2019.10.16
% - pilot test added: 2019.10.17
% - walsh-hadamard sequence added: 2019.10.18
% - channel coding added: 2019.10.21
% - rate matching added: 2019.10.22
% - structure updated: 2019.10.23

% set test mode
test_dd_pilot = false;  % allocate whole resource block to pilot
test_real_ch = false;   % real channel with ici into consideration

% set parameters
SNRCnt = 10;
SNRStep = 2;
ErrorStop = 400;
ErrorBreak = 0.001;

% set bandwidth list
bw_list = [1.4 3, 5, 10, 15, 20];

% open files
list_file = './nw_sim_list.lst';
result_file = './nw_sim_result.res';
fp1 = fopen(list_file, 'r');
fp2 = fopen(result_file, 'a');

% start simulation
while 1
    % read file
    file_read = fgets(fp1);
    if ~ischar(file_read)
        break
    end
    if strcmp(file_read(1:2), '//')
        continue
    end
    
    % parse parameters
    nw_parse_prm = nw_sim_parse_prm(file_read);
    
    % evaluate the parameters
    carrier_freq_mhz = nw_parse_prm.fc;
    velocity_kmh = nw_parse_prm.vel;
    delay_spread_rms_us = nw_parse_prm.delay; % us
    idx_fading = nw_parse_prm.ch;
    snr_db_start = nw_parse_prm.snr;
    idx_bandwidth = find(bw_list == nw_parse_prm.bw, 1);
    len_tb_bit = nw_parse_prm.tblen;
    mcs = nw_parse_prm.mcs;
    num_sim = nw_parse_prm.nsim;
    chest_option = nw_parse_prm.chest;
    cheq_option = nw_parse_prm.cheq;
    waveform = nw_parse_prm.wave;
    
    % set equalization option
    nw_num = nw_num_prm(idx_bandwidth, waveform, chest_option);
    nw_cc = nw_cc_prm(len_tb_bit);
    nw_rm = nw_rm_prm(len_tb_bit, mcs, nw_num, nw_cc);
    nw_sim = nw_sim_prm(len_tb_bit, num_sim, nw_num, nw_cc, nw_rm);
    
    % set test option
    test_option.papr = false;               % papr test
    test_option.ch_mse = false;             % channel mmse test
    test_option.sym_err_var = false;        % symbol error variance test
    test_option.ch_edge_interp = false;     % tf-channel edge interpolation test (no use)
    test_option.partial_reception = 14;   % partial reception test for latency enhancement (number of symbols received)
    test_option.otfs_map_plan = 3;          % pilot resource mapping test (refer to 'otfs_sym_map_r2.m')
    test_option.otfs_pilot_impulse_pwr_reduction = false;       % valid only when impulse pilot is used
    test_option.otfs_pilot_spread_seq = zadoffChuSeq(1,37);     % pilot spread seq. (zadoff-chu sequence spreading)
    test_option.otfs_pilot_seq_ones = ones(37, 1);              % pilot spread seq. (all ones)
%     test_option.otfs_pilot_spread_seq = exp(-1i*pi*25*(0:36).*(1:37)/37);
    
    % check test option
    if test_option.otfs_pilot_impulse_pwr_reduction && ~(test_option.otfs_map_plan == 1 || test_option.otfs_map_plan == 2 || test_option.otfs_map_plan == 3)
        error('To use ''otfs_pilot_impulse_pwr_reduction'', ''otfs_map_plan'' must be one of these: {1, 2, 3}')
    end
    
%     % generate walsh-hadamard sequence
%     walsh_seq = 1;
%     for idx = 1 : log2(16)
%         walsh_seq_update = kron([1 1; 1 -1], walsh_seq);
%         walsh_seq = walsh_seq_update;
%     end
%     test_option.walsh_seq = walsh_seq;
    
    fprintf('RUNNING CASE: %s', file_read);
    
    nw_sim_result = zeros(SNRCnt, 3);
    for snr_idx = 1 : SNRCnt
        snr_db = snr_db_start + SNRStep * (snr_idx - 1);
        nw_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, snr_db, delay_spread_rms_us);
        
        % initialize counter
        total_error = 0;
        sim_cnt = 0;
        sum_papr = 0;
        sum_ch_mse = 0;
        sum_sym_err_var = 0;
        
        for sim_idx = 1 : num_sim
            
            % simulate single packet
            if strcmp(waveform, 'ofdm')    % ofdm
                [pkt_error, tx_papr, ch_mse] = ofdm_single_run_r2(nw_sim, nw_cc, nw_rm, nw_num, snr_db, nw_ch, chest_option, cheq_option, test_option);
            else                           % otfs
                [pkt_error, tx_papr, ch_mse, sym_err_var] = otfs_single_run_r4(nw_sim, nw_cc, nw_rm, nw_num, snr_db, nw_ch, chest_option, cheq_option, test_option);
            end
            
            % count packet error
            total_error = total_error + double(pkt_error);
            sim_cnt = sim_cnt + 1;
            if total_error > ErrorStop
                break
            end
            
            % sum papr per packet
            if test_option.papr
                sum_papr = sum_papr + tx_papr;
            end
            
            % sum channel mse
            if test_option.ch_mse
                sum_ch_mse = sum_ch_mse + ch_mse;
            end
            
            % sum channel mse
            if test_option.sym_err_var
                sum_sym_err_var = sum_sym_err_var + sym_err_var;
            end
            
        end
        
        % compute bit error rate
        nw_sim_result(snr_idx, 1) = snr_db;                % snr
        nw_sim_result(snr_idx, 2) = sim_cnt;                % total number of packets
        nw_sim_result(snr_idx, 3) = total_error / sim_cnt;  % per
        
        % save simulation result
        fprintf('    SNR: %4.1f     NUM PKTS: %7d   PER: %9.6f', nw_sim_result(snr_idx, :));
        fprintf(fp2, '%5.1f %10d %10.6f', nw_sim_result(snr_idx, :)');
        
        if test_option.papr
            papr_db = 10 * log10(sum_papr / sim_cnt);       % calculate papr
            fprintf('   PAPR: %6.3f', papr_db);
            fprintf(fp2, ' %10.6f', papr_db);
        end
        
        if test_option.ch_mse
            ch_mmse = sqrt(sum_ch_mse / sim_cnt);           % calculate channel mmse
            fprintf('   CH MMSE: %6.3f', ch_mmse);
            fprintf(fp2, ' %10.6f', ch_mmse);
        end
        
        if test_option.sym_err_var
            mean_sym_err_var = sqrt(sum_sym_err_var / sim_cnt);           % calculate symbol error variance
            fprintf('   SYM ERR VAR: %6.3f   NOISE VAR: %12.9f', mean_sym_err_var, (nw_num.num_subc_usr/nw_num.nfft)*(10^((-0.1)*snr_db)));
            fprintf(fp2, ' %10.6f', mean_sym_err_var);
        end
        
        fprintf('\n');
        fprintf(fp2, '\n');
        
        % set stopping criteria
        if nw_sim_result(snr_idx, 3) < ErrorBreak
            disp('BREAK');
            break;
        end
    end
    
    % plot
    figure
    semilogy(nw_sim_result(:, 1), nw_sim_result(:, 3), '-bo')
    xlabel('SNR (dB)'); ylabel('PER');
    title('QAM over Fading Channel');
    axis([nw_sim_result(1, 1) - 1 snr_db + 1 1e-4 1])
    grid minor
end

fclose(fp1);
fclose(fp2);
