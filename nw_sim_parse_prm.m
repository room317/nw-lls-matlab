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
prm_scs_idx = find(contains(prm_set, 'SCS'), 1);
prm_slot_idx = find(contains(prm_set, 'SLOT'), 1);
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
prm_chereal_idx = find(contains(prm_set, 'REAL'), 1);
prm_cheperfect_idx = find(contains(prm_set, 'PERFECT'), 1);
prm_cheddimpulse_idx = find(contains(prm_set, 'DDIMPULSE'), 1);
prm_cheddzc_idx = find(contains(prm_set, 'DDZC'), 1);
prm_cheddrandom_idx = find(contains(prm_set, 'DDRANDOM'), 1);
prm_cheddgolayserial_idx = find(contains(prm_set, 'DDGOLAYSER'), 1);
prm_cheddgolayparallel_idx = find(contains(prm_set, 'DDGOLAYPAR'), 1);
prm_cheddgolaydiag_idx = find(contains(prm_set, 'DDGOLAYDIAG'), 1);
prm_chetfltedown_idx = find(contains(prm_set, 'TFLTEDOWN'), 1);
prm_chetflteup_idx = find(contains(prm_set, 'TFLTEUP'), 1);
prm_chetfnr_idx = find(contains(prm_set, 'TFNR'), 1);
prm_tfeqzf_idx = find(contains(prm_set, 'TFEQZF'), 1);
prm_tfeqmmse_idx = find(contains(prm_set, 'TFEQMMSE'), 1);
prm_ddeqzf_idx = find(contains(prm_set, 'DDEQZF'), 1);
prm_ddeqmmse_idx = find(contains(prm_set, 'DDEQMMSE'), 1);
prm_su_idx = find(contains(prm_set, 'SU'), 1);
prm_mu_idx = find(contains(prm_set, 'MU'), 1);

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
waveform_set = [ ...
    ~isempty(prm_ofdm_idx)
    ~isempty(prm_otfs_idx)];
if sum(double(waveform_set)) ~= 1
    error('Waveform setting. Choose one of these: {OFDM, OTFS}')
elseif ~isempty(prm_ofdm_idx)
    prm_wave = 'ofdm'; % 'OFDM'
else
    prm_wave = 'otfs'; % 'OTFS'
end

% bandwidth
if isempty(prm_bw_idx)
    error('Bandwidth setting. Choose one of these: {BW5, BW10, BW15, BW20, BW25, BW30, BW40, BW50, BW60, BW80, BW90, BW100}')
else
    prm_bw_str = prm_set{prm_bw_idx};
    prm_bw = sscanf(prm_bw_str(3 : end), '%f');
end

% subcarrier spacing
if isempty(prm_scs_idx)
    error('Subcarrier spacing setting. Choose one of these: {SCS15, SCS30, SCS60}')
else
    prm_scs_str = prm_set{prm_scs_idx};
    prm_scs = sscanf(prm_scs_str(4 : end), '%f');
end

% number of slots (1 slot = 14 ofdm symbols)
if isempty(prm_slot_idx)
    error('Slot setting example: SLOT1, SLOT2, SLOT3 ...')
else
    prm_slot_str = prm_set{prm_slot_idx};
    prm_slot = sscanf(prm_slot_str(5 : end), '%f');
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
ch_set = [ ...
    ~isempty(prm_rmalos_idx)
    ~isempty(prm_rmanlos_idx)
    ~isempty(prm_umilos_idx)
    ~isempty(prm_uminlos_idx)
    ~isempty(prm_umalos_idx)
    ~isempty(prm_umanlos_idx)
    ~isempty(prm_eva_idx)
    ~isempty(prm_etu_idx)
    ~isempty(prm_tdla_idx)
    ~isempty(prm_tdlb_idx)
    ~isempty(prm_tdlc_idx)
    ~isempty(prm_tdld_idx)
    ~isempty(prm_tdle_idx)
    ~isempty(prm_test_idx)];
if sum(double(ch_set)) ~= 1
    error('Fading channel setting: choose one of these models: RMALOS, RMANLOS, UMILOS, UMINLOS, UMALOS, UMANLOS, EVA, ETU, TDLA, TDLB, TDLC, TDLD, TDLE, TEST.')
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

% channel estimation
chest_set = [ ...
    ~isempty(prm_chereal_idx)
    ~isempty(prm_cheperfect_idx)
    ~isempty(prm_cheddimpulse_idx)
    ~isempty(prm_cheddzc_idx)
    ~isempty(prm_cheddrandom_idx)
    ~isempty(prm_cheddgolayserial_idx)
    ~isempty(prm_cheddgolayparallel_idx)
    ~isempty(prm_cheddgolaydiag_idx)
    ~isempty(prm_chetfltedown_idx)
    ~isempty(prm_chetflteup_idx)
    ~isempty(prm_chetfnr_idx)];
if sum(double(chest_set)) ~= 1
    error('Channel estimation setting: choose one of these: REAL, PERFECT, DDIMPULSE, DDZC, DDRANDOM, DDGOLAYSER, DDGOLAYPAR, DDGOLAYDIAG, TFLTEDOWN, TFLTEUP, TFNR.')
elseif ~isempty(prm_chereal_idx)
    prm_chest = 'real';             % channel from the fading block (interference exists)
elseif ~isempty(prm_cheperfect_idx)
    prm_chest = 'perfect';          % channel from the whole received signals (interferece cancelled)
elseif ~isempty(prm_cheddimpulse_idx)
    prm_chest = 'dd_impulse';           % channel from dd-domain singletone
elseif ~isempty(prm_cheddzc_idx)
    prm_chest = 'dd_zc';                % channel from dd-domain zadoff-chu sequence
elseif ~isempty(prm_cheddrandom_idx)
    prm_chest = 'dd_random';            % channel from dd-domain random(default: zadoff-chu) sequence
elseif ~isempty(prm_cheddgolayserial_idx)
    prm_chest = 'dd_golay_serial';      % channel from dd-domain golay complementary sequence (in serial layout)
elseif ~isempty(prm_cheddgolayparallel_idx)
    prm_chest = 'dd_golay_parallel';    % channel from dd-domain golay complementary sequence (in parallel layout)
elseif ~isempty(prm_cheddgolaydiag_idx)
    prm_chest = 'dd_golay_diag';        % channel from dd-domain golay complementary sequence (in diagonal layout)
elseif ~isempty(prm_chetfltedown_idx)
    prm_chest = 'tf_ltedown';       % channel from tf-domain lte downlink pilots
elseif ~isempty(prm_chetflteup_idx)
    prm_chest = 'tf_lteup';         % channel from tf-domain lte uplink pilots
else
    prm_chest = 'tf_nr';            % channel from tf-domain nr downlink pilots
end

% channel equalization
cheq_set = [ ...
    ~isempty(prm_tfeqzf_idx)
    ~isempty(prm_tfeqmmse_idx)
    ~isempty(prm_ddeqzf_idx)
    ~isempty(prm_ddeqmmse_idx)];
if sum(double(cheq_set)) ~= 1
    error('Channel equalization setting: choose one of these: TFEQZF, TFEQMMSE, DDEQZF, DDEQMMSE.')
elseif ~isempty(prm_tfeqzf_idx)
    prm_cheq = 'tfeq_zf';       % tf-domain equalization
elseif ~isempty(prm_tfeqmmse_idx)
    prm_cheq = 'tfeq_mmse';     % tf-domain equalization
elseif ~isempty(prm_ddeqzf_idx)
    prm_cheq = 'ddeq_zf';       % dd-domain equalization
else
    prm_cheq = 'ddeq_mmse';     % dd-domain equalization
end

% multi-user simulation
usr_set = [ ...
    ~isempty(prm_su_idx)
    ~isempty(prm_mu_idx)];
if sum(double(usr_set)) ~= 1
    error('Number of users: choose one of these: SU, MU.')
elseif ~isempty(prm_su_idx)
    prm_usr = 'su';  % single-user
else
    prm_usr = 'mu';  % multi-user
end

% check errors
if strcmp(prm_wave, 'ofdm') && strncmp(prm_chest, 'dd_', 3)
    error('OFDM cannot use dd-domain pilots.')
end

% output
nw_parse_prm.wave = prm_wave;
nw_parse_prm.bw = prm_bw;
nw_parse_prm.scs = prm_scs;
nw_parse_prm.slot = prm_slot;
nw_parse_prm.fc = prm_fc;
nw_parse_prm.vel = prm_vel;
nw_parse_prm.delay = prm_delay;
nw_parse_prm.ch = prm_ch;
nw_parse_prm.mcs = prm_mcs;
nw_parse_prm.tblen = prm_tblen;
nw_parse_prm.nsim = prm_nsim;
nw_parse_prm.snr = prm_snr;
nw_parse_prm.chest = prm_chest;
nw_parse_prm.cheq = prm_cheq;
nw_parse_prm.usr = prm_usr;

end
