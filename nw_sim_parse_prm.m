function nw_parse_prm = nw_sim_parse_prm(file_read)

% bandwidth(mhz) (BW1.4, BW3, BW5, BW10, BW15, BW20)
% carrier(mhz) (F5800)
% velocity(km/h) (Vx)
% fading (RMA, RMANLOS, UMI, UMINLOS, UMA, UMANLOS, EVA, ETU)
% mcs (MCS1 ~ MCS15)
% transfer block length (LENx ~)
% number of simulations (SIMx ~)
% ebno(db) (SNRx)

% ex: BW10 FC5800 V30 RMA MCS10 LEN500 SIM100 SNR10

prm_cell = textscan(file_read, '%s');
prm_set = prm_cell{1};

prm_ofdm_idx = find(contains(prm_set, 'OFDM'), 1);
prm_otfs_idx = find(contains(prm_set, 'OTFS'), 1);
prm_bw_idx = find(contains(prm_set, 'BW'), 1);  % cellfun('isempty', strfind(prm_set, 'BW'));
prm_fc_idx = find(contains(prm_set, 'FC'), 1);
prm_vel_idx = find(contains(prm_set, 'V'), 1);
prm_delay_idx = find(contains(prm_set, 'DS'), 1);
prm_rmalos_idx = find(contains(prm_set, 'RMALOS'), 1);
prm_rmanlos_idx = find(contains(prm_set, 'RMANLOS'), 1);
prm_umilos_idx = find(contains(prm_set, 'UMILOS'), 1);
prm_uminlos_idx = find(contains(prm_set, 'UMINLOS'), 1);
prm_umalos_idx = find(contains(prm_set, 'UMALOS'), 1);
prm_umanlos_idx = find(contains(prm_set, 'UMANLOS'), 1);
prm_eva_idx = find(contains(prm_set, 'EVA'), 1);
prm_etu_idx = find(contains(prm_set, 'ETU'), 1);
prm_tdla_idx = find(contains(prm_set, 'TDLA'), 1);
prm_tdlb_idx = find(contains(prm_set, 'TDLB'), 1);
prm_tdlc_idx = find(contains(prm_set, 'TDLC'), 1);
prm_tdld_idx = find(contains(prm_set, 'TDLD'), 1);
prm_tdle_idx = find(contains(prm_set, 'TDLE'), 1);
prm_test_idx = find(contains(prm_set, 'TEST'), 1);
prm_mcs_idx = find(contains(prm_set, 'MCS'), 1);
prm_tblen_idx = find(contains(prm_set, 'LEN'), 1);
prm_nsim_idx = find(contains(prm_set, 'SIM'), 1);
prm_snr_idx = find(contains(prm_set, 'SNR'), 1);

% prm_bw_idx = find(not(cellfun('isempty', strfind(prm_set, 'BW'))), 1);
% prm_fc_idx = find(not(cellfun('isempty', strfind(prm_set, 'F'))), 1);
% prm_vel_idx = find(not(cellfun('isempty', strfind(prm_set, 'V'))), 1);
% prm_rma_idx = find(not(cellfun('isempty', strfind(prm_set, 'RMA'))), 1);
% prm_rmanlos_idx = find(not(cellfun('isempty', strfind(prm_set, 'RMANLOS'))), 1);
% prm_umi_idx = find(not(cellfun('isempty', strfind(prm_set, 'UMI'))), 1);
% prm_uminlos_idx = find(not(cellfun('isempty', strfind(prm_set, 'UMINLOS'))), 1);
% prm_uma_idx = find(not(cellfun('isempty', strfind(prm_set, 'UMA'))), 1);
% prm_umanlos_idx = find(not(cellfun('isempty', strfind(prm_set, 'UMANLOS'))), 1);
% prm_eva_idx = find(not(cellfun('isempty', strfind(prm_set, 'EVA'))), 1);
% prm_etu_idx = find(not(cellfun('isempty', strfind(prm_set, 'ETU'))), 1);
% prm_mcs_idx = find(not(cellfun('isempty', strfind(prm_set, 'MCS'))), 1);
% prm_tblen_idx = find(not(cellfun('isempty', strfind(prm_set, 'LEN'))), 1);
% prm_nsim_idx = find(not(cellfun('isempty', strfind(prm_set, 'SIM'))), 1);
% prm_snr_idx = find(not(cellfun('isempty', strfind(prm_set, 'SNR'))), 1);

% waveform
if sum(double([isempty(prm_ofdm_idx) isempty(prm_otfs_idx)])) ~= 1
    error('Waveform setting: choose one of these waveforms (OFDM, OTFS).')
elseif ~isempty(prm_ofdm_idx)
    prm_wave = 1; % 'OFDM'
else
    prm_wave = 2; % 'OTFS'
end

% bandwidth
if isempty(prm_bw_idx)
    error('Bandwidth setting example: BW1.4, BW3, BW5, BW10, BW15, BW20')
else
    prm_bw_str = prm_set{prm_bw_idx};
    prm_bw = sscanf(prm_bw_str(3 : end), '%f');
end

% carrier frequency
if isempty(prm_fc_idx)
    error('Carrier frequency setting example[MHz]: FC5800')
else
    prm_fc_str = prm_set{prm_fc_idx};
    prm_fc = sscanf(prm_fc_str(3 : end), '%f');
end

% velocity
if isempty(prm_vel_idx)
    error('Velocity setting example[km/h]: V30')
else
    prm_vel_str = prm_set{prm_vel_idx};
    prm_vel = sscanf(prm_vel_str(2 : end), '%f');
end

% velocity
if isempty(prm_delay_idx)
    error('RMS delay spread setting example[us]: DS0.1')
else
    prm_delay_str = prm_set{prm_delay_idx};
    prm_delay = sscanf(prm_delay_str(3 : end), '%f');
end

% fading channel
if sum(double([isempty(prm_rmalos_idx) isempty(prm_rmanlos_idx) isempty(prm_umilos_idx) isempty(prm_uminlos_idx) isempty(prm_umalos_idx) isempty(prm_umanlos_idx) isempty(prm_eva_idx) isempty(prm_etu_idx) isempty(prm_tdla_idx) isempty(prm_tdlb_idx) isempty(prm_tdlc_idx) isempty(prm_tdld_idx) isempty(prm_tdle_idx) isempty(prm_test_idx)])) ~= 13
    error('Fading channel setting: choose one of these models (RMALOS, RMANLOS, UMILOS, UMINLOS, UMALOS, UMANLOS, EVA, ETU, TDLA, TDLB, TDLC, TDLD, TDLE, TEST).')
elseif ~isempty(prm_rmalos_idx)
    prm_ch = 1; % 'RMa LOS'
elseif ~isempty(prm_rmanlos_idx)
    prm_ch = 2; % 'RMa NLOS'
elseif ~isempty(prm_umilos_idx)
    prm_ch = 3; % 'UMi LOS'
elseif ~isempty(prm_uminlos_idx)
    prm_ch = 4; % 'UMi NLOS'
elseif ~isempty(prm_umalos_idx)
    prm_ch = 5; % 'UMa LOS'
elseif ~isempty(prm_umanlos_idx)
    prm_ch = 6; % 'UMa NLOS'
elseif ~isempty(prm_eva_idx)
    prm_ch = 7; % 'EVA'
elseif ~isempty(prm_etu_idx)
    prm_ch = 8; % 'ETU'
elseif ~isempty(prm_tdla_idx)
    prm_ch = 9; % 'TDL-A'
elseif ~isempty(prm_tdlb_idx)
    prm_ch = 10; % 'TDL-B'
elseif ~isempty(prm_tdlc_idx)
    prm_ch = 11; % 'TDL-C'
elseif ~isempty(prm_tdld_idx)
    prm_ch = 12; % 'TDL-D'
elseif ~isempty(prm_tdle_idx)
    prm_ch = 13; % 'TDL-E'
else
    prm_ch = 14; % 'TEST'
end

% mcs
if isempty(prm_mcs_idx)
    error('MCS setting example: MCS1 ~ MCS15')
else
    prm_mcs_str = prm_set{prm_mcs_idx};
    prm_mcs = sscanf(prm_mcs_str(4 : end), '%f');
end

% transfer block length
if isempty(prm_tblen_idx)
    error('TB length setting example: LEN4000')
else
    prm_tblen_str = prm_set{prm_tblen_idx};
    prm_tblen = sscanf(prm_tblen_str(4 : end), '%f');
end

% number of simulations
if isempty(prm_nsim_idx)
    error('Simulation number setting example: SIM1000')
else
    prm_nsim_str = prm_set{prm_nsim_idx};
    prm_nsim = sscanf(prm_nsim_str(4 : end), '%f');
end

% snr
if isempty(prm_snr_idx)
    error('SNR setting example[Eb/No]: SNR10')
else
    prm_snr_str = prm_set{prm_snr_idx};
    prm_snr = sscanf(prm_snr_str(4 : end), '%f');
end

% output
nw_parse_prm.wave = prm_wave;
nw_parse_prm.bw = prm_bw;
nw_parse_prm.fc = prm_fc;
nw_parse_prm.vel = prm_vel;
nw_parse_prm.delay = prm_delay;
nw_parse_prm.ch = prm_ch;
nw_parse_prm.mcs = prm_mcs;
nw_parse_prm.tblen = prm_tblen;
nw_parse_prm.nsim = prm_nsim;
nw_parse_prm.snr = prm_snr;

end
