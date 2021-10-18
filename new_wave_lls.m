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

% set parameters
SNRCnt = 10;
SNRStep = 2;
ErrorStop = 400;
ErrorBreak = 0.001;

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
    bw_mhz = nw_parse_prm.bw;
    scs_khz = nw_parse_prm.scs;
    num_slot = nw_parse_prm.slot;
    len_tb_bit = nw_parse_prm.tblen;
    mcs = nw_parse_prm.mcs;
    num_sim = nw_parse_prm.nsim;
    chest_option = nw_parse_prm.chest;
    cheq_option = nw_parse_prm.cheq;
    waveform = nw_parse_prm.wave;
    usr_option = nw_parse_prm.usr;
    
    % check parameters
    if ~(len_tb_bit > 0)
        error('Check TB length should be greater than 0.')
    end
    
    % set test simulation option
    sim_option.override = false;    % override simulation options when true
    sim_option.num_rb = 11;         % number of total rbs
    sim_option.scs_khz = 60;        % subcarrier spacing (khz)
    sim_option.num_slot = 4;        % number of total slots
    sim_option.Qm = 4;              % number of bits per qam symbol
    sim_option.coderate = 666;      % number of information bits per 1024 bits
    
    % set test option
    test_option.license = true;                     % true: with toolbox license, false: no toolbox license
    test_option.fading_only = false;                % no awgn noise when true
    test_option.awgn = false;                       % no fading when true
    test_option.papr = false;                       % papr test
    test_option.ch_mse = false;                     % channel mmse test
    test_option.sym_err_var = false;                % symbol error variance test
    test_option.ch_tf_edge_interp = false;          % tf-channel edge interpolation test (no use)
    test_option.ch_dd_guard_interp = false;         % guard interpolation for channel estimation
    test_option.ch_est_xcorr_prune = false;         % channel estimation pilot resource pruning after sequence pilot xcorr
    test_option.otfs_pilot_pwr_set = false;         % valid only when otfs_pilot_plan is 'impulse'
    test_option.zc_seq_len = 37;                                    % zadoff-chu sequence length (37)
    test_option.zc_seq = zadoffChuSeq(1,test_option.zc_seq_len);    % zadoff-chu sequence
    test_option.golay_seq_len = 32;                                 % golay sequence length (32, 64, 128)
    if test_option.license
        [Ga, Gb] = wlanGolaySequence(test_option.golay_seq_len);    % complementary golay sequence
    else
        Ga = [];    % license problem
        Gb = [];    % license problem
    end
    test_option.golay_seq_a = Ga;           % complementary golay sequence
    test_option.golay_seq_b = Gb;           % complementary golay sequence
    test_option.fulltap_eq = false;         % use full-tap real channel for equalization
    test_option.common_usr_ch = true;       % use common channel per user
    test_option.gpu_flag = false;           % use gpu for real channel generation
    test_option.pilot_only = false;         % use null data
    test_option.ch_clip = 'none';           % clipping real channel ('none': no clipping, 'center': channel for symbol at center, 'random': channel for random symbol)
%     test_option.otfs_pilot_spread_seq = exp(-1i*pi*25*(0:36).*(1:37)/37);
    test_option.rm_version = 2;             % 1: exactly follows standard, 2: faster vectorized version
    test_option.rv_idx = 0;                 % redundancy version numer
    test_option.max_tbs_calc = false;       % skip sfft calc. when max tbs calc.
    test_option.custom_nr_pilot = false;    % true for activating custom nr pilots
    test_option.qam_modem_toolbox = false;  % true to use lte toolbox qam mod/demod
    
    % check test option
    if test_option.otfs_pilot_pwr_set && ~strcmp(chest_option, 'dd_impulse')
        error('To use ''otfs_pilot_pwr_set'', ''chest_option'' must be ''dd_impulse''.')
    end
    
    if ~test_option.license && strncmp(chest_option, 'dd_golay_', 9)
        error('To use ''DDGOLAY..'' option, ''test_option.license'' option must be true.')
    end
    
    % set equalization option
    nw_cc = nw_cc_prm_r1(len_tb_bit);
    nw_num = nw_num_prm_r1(carrier_freq_mhz, scs_khz, bw_mhz, num_slot, waveform, chest_option, usr_option, sim_option, test_option);
    nw_rm = nw_rm_prm_r1(mcs, nw_cc, nw_num, sim_option, test_option);
    nw_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, delay_spread_rms_us);
    
    % set additional test option
    test_option.part_rx = nw_num.num_ofdmsym;       % partial reception test for latency enhancement (number of symbols received)
    
%     % generate walsh-hadamard sequence
%     walsh_seq = 1;
%     for idx = 1 : log2(16)
%         walsh_seq_update = kron([1 1; 1 -1], walsh_seq);
%         walsh_seq = walsh_seq_update;
%     end
%     test_option.walsh_seq = walsh_seq;
    
    % create crc generator and detector objects
    tx_crc_a = comm.CRCGenerator(nw_cc.gCRC24A);
    tx_crc_b = comm.CRCGenerator(nw_cc.gCRC24B);
    rx_crc_a = comm.CRCDetector(nw_cc.gCRC24A);
%     rx_crc_b = comm.CRCDetector(nw_cc.gCRC24B);
    
    % create a rayleigh fading channel object
    if nw_ch.los
        if test_option.license
            fading_ch = nrTDLChannel( ...
                'DelayProfile', 'TDL-E', ...
                'DelaySpread', delay_spread_rms_us*1e-6, ...
                'MaximumDopplerShift', nw_ch.maximum_doppler_shift, ...
                'NumTransmitAntennas', 1, ...
                'NumReceiveAntennas', 1, ...
                'SampleRate', nw_num.sample_rate, ...
                'KFactorScaling', true, ...
                'KFactor', 10*log10(nw_ch.k_factor),...
                'NormalizePathGains', true);
        else
            fading_ch = comm.RicianChannel(...
                'SampleRate', num.sample_rate,...
                'PathDelays', ch.path_delays,...
                'AveragePathGains', ch.average_path_gains,...
                'KFactor', ch.k_factor,...
                'DirectPathDopplerShift', ch.maximum_doppler_shift,...
                'MaximumDopplerShift', ch.maximum_doppler_shift,...
                'DopplerSpectrum', ch.doppler_spectrum,...
                'PathGainsOutputPort', true, ...
                'NormalizePathGains', true);
        end
    else
        if test_option.license
            fading_ch = nrTDLChannel( ...
                'DelayProfile', 'TDL-C', ...
                'DelaySpread', delay_spread_rms_us*1e-6, ...
                'MaximumDopplerShift', nw_ch.maximum_doppler_shift, ...
                'NumTransmitAntennas', 1, ...
                'NumReceiveAntennas', 1, ...
                'SampleRate', nw_num.sample_rate, ...
                'NormalizePathGains', true);

%             fading_ch = nrTDLChannel( ...
%                 'DelayProfile', 'Custom', ...
%                 'AveragePathGains', nw_ch.average_path_gains,...
%                 'PathDelays', nw_ch.path_delays, ...
%                 'MaximumDopplerShift', nw_ch.maximum_doppler_shift, ...
%                 'NumTransmitAntennas', 1, ...
%                 'NumReceiveAntennas', 1, ...
%                 'SampleRate', nw_num.sample_rate, ...
%                 'NormalizePathGains', true);
        else
            fading_ch = comm.RayleighChannel(...
                'SampleRate', nw_num.sample_rate, ...
                'PathDelays', nw_ch.path_delays, ...
                'AveragePathGains', nw_ch.average_path_gains, ...
                'NormalizePathGains', true, ...
                'MaximumDopplerShift', nw_ch.maximum_doppler_shift, ...
                'DopplerSpectrum', nw_ch.doppler_spectrum, ...
                'PathGainsOutputPort', true);
        end
    end
    
    tic
    fprintf('RUNNING CASE: %s', file_read);
    
    nw_sim_result = zeros(SNRCnt, nw_num.num_usr+2);
    for snr_idx = 1:SNRCnt
        snr_db = snr_db_start+SNRStep*(snr_idx-1);
        
        % initialize counter
        total_error = zeros(1, nw_num.num_usr);
        sim_cnt = 0;
        sum_papr = 0;
        sum_ch_mse = zeros(1, nw_num.num_usr);
        sum_sym_err_var = zeros(1, nw_num.num_usr);
        sum_ch_pwr = zeros(nw_num.num_delay_usr, nw_num.num_doppler_usr);
        
        for sim_idx = 1:num_sim
            
            % simulate single packet
            if strcmp(waveform, 'ofdm')    % ofdm
                [pkt_error, tx_papr, ch_mse, sym_err_var] = ...
                    ofdm_dnlink_singlerun_r3( ...
                    test_option.rv_idx, nw_cc, nw_rm, nw_num, snr_db, ...
                    tx_crc_a, tx_crc_b, rx_crc_a, fading_ch, ...
                    chest_option, cheq_option, test_option);
%                 [pkt_error, tx_papr, ch_mse, sym_err_var] = ...
%                     ofdm_dnlink_singlerun_r2( ...
%                     nw_sim, nw_cc, nw_rm, nw_num, snr_db, ...
%                     tx_crc, rx_crc, turbo_enc, turbo_dec, fading_ch, ...
%                     chest_option, cheq_option, test_option);
%                 [pkt_error, tx_papr, ch_mse, sym_err_var] = ...
%                     test_ofdm_otfs_r0( ...
%                     nw_sim, nw_cc, nw_rm, nw_num, snr_db, ...
%                     tx_crc, rx_crc, turbo_enc, turbo_dec, fading_ch, ...
%                     chest_option, cheq_option, test_option);
            else                           % otfs
                [pkt_error, tx_papr, ch_mse, sym_err_var, ch_est_rbs_dd] = ...
                    otfs_dnlink_singlerun_r2( ...
                    test_option.rv_idx, nw_cc, nw_rm, nw_num, snr_db, ...
                    tx_crc_a, tx_crc_b, rx_crc_a, fading_ch, ...
                    chest_option, cheq_option, test_option);
%                 [pkt_error, tx_papr, ch_mse, sym_err_var, ch_est_rbs_dd] = ...
%                     otfs_dnlink_singlerun_r1( ...
%                     nw_sim, nw_cc, nw_rm, nw_num, snr_db, ...
%                     tx_crc, rx_crc, turbo_enc, turbo_dec, fading_ch, ...
%                     chest_option, cheq_option, test_option);
            end
            
            % count packet error
            total_error = total_error+double(pkt_error);
            sim_cnt = sim_cnt+1;
            if min(total_error)>ErrorStop
                break
            end
            
            % sum papr per packet
            if test_option.papr
                sum_papr = sum_papr+tx_papr;
            end
            
            % sum channel mse
            if test_option.ch_mse
                sum_ch_mse = sum_ch_mse+ch_mse;
            end
            
            if test_option.ch_mse && strcmp(waveform, 'otfs')
                sum_ch_pwr = sum_ch_pwr+sum(abs(ch_est_rbs_dd).^2, [3 4]);
            end
            
            % sum symbol errors
            if test_option.sym_err_var
                sum_sym_err_var = sum_sym_err_var+sym_err_var;
            end
            
        end
        
        % compute bit error rate
        nw_sim_result(snr_idx, 1) = snr_db;                % snr
        nw_sim_result(snr_idx, 2) = sim_cnt;                % total number of packets
        nw_sim_result(snr_idx, 3:end) = total_error/sim_cnt;  % per
        
        % save simulation result
        fprintf('    SNR: %4.1f     NUM PKTS: %7d   PER:', nw_sim_result(snr_idx, 1:2))
        fprintf(' %9.6f', nw_sim_result(snr_idx, 3:end))
        fprintf(fp2, '%5.1f %10d', nw_sim_result(snr_idx, 1:2)');
        fprintf(fp2, ' %10.6f', nw_sim_result(snr_idx, 3:end)');
        
        if test_option.papr
            papr_db = 10*log10(sum_papr/sim_cnt);       % calculate papr
            fprintf('   PAPR:')
            fprintf(' %6.3f', papr_db)
            fprintf(fp2, ' %10.6f', papr_db);
        end
        
        if test_option.ch_mse
            ch_rmse = sqrt(sum_ch_mse/sim_cnt);           % calculate channel mmse
            fprintf('   CH RMSE:')
            fprintf(' %6.3f', ch_rmse)
            fprintf(fp2, ' %10.6f', ch_rmse);
        end
        
        if test_option.sym_err_var
            mean_sym_err_var = sqrt(sum_sym_err_var/sim_cnt);           % calculate symbol error variance
            fprintf('   NOISE VAR: %12.9f', (nw_num.num_subc_usr/nw_num.num_fft)*(10^((-0.1)*snr_db)))
            fprintf('   SYM ERR VAR:')
            fprintf(' %6.3f', mean_sym_err_var)
            fprintf(fp2, ' %10.6f', mean_sym_err_var);
        end
        
        fprintf('\n');
        fprintf(fp2, '\n');
        
        % set stopping criteria
        if (max(nw_sim_result(snr_idx, 3:end)) < ErrorBreak) && ~(test_option.papr || test_option.ch_mse || test_option.sym_err_var)
            disp('BREAK');
            break;
        end
    end
    toc
    
    % plot
    for idx_usr = 1:nw_num.num_usr
        figure
        semilogy(nw_sim_result(:, 1), nw_sim_result(:, 2+idx_usr), '-bo')
        xlabel('SNR (dB)'), ylabel('PER'), title('QAM over Fading Channel')
        axis([nw_sim_result(1, 1) - 1 snr_db + 1 1e-4 1])
        grid minor
    end
    
    if test_option.ch_mse && strcmp(waveform, 'otfs')
        for idx_usr = 1:nw_num.num_usr
            figure
            mesh(1:nw_num.num_doppler_usr, 1:nw_num.num_delay_usr, sqrt(fftshift(fftshift(sum_ch_pwr/(sim_cnt*nw_rm.N_RB*nw_num.num_usr), 1), 2)))
            xlabel('Doppler'), ylabel('Delay'), zlabel('Average Channel Amplitude'), title('Channel Estimation')
        end
    end
end

fclose(fp1);
fclose(fp2);
