function nw_ch = nw_ch_prm(carrier_freq_mhz, velocity_kmh, idx_fading, ebno_db, delay_spread_rms_us)

% channel parameters
% ref. IST-4-027756 WINNER II D1.1.2 V1.2
% ref. 3gpp tr 38.901
% ref. 3gpp ts 36.101 
% ref. 3gpp ts 36.104 
% carrier_freq_mhz = 4000;
% velocity_kmh = 30;

% D1-Rural Macro (RMa) LOS (Winner II)
% power_profile_rma_los_db = [-15];
% delay_profile_rma_los_us = [0.005];
power_profile_rma_los_db = [-15 -15.5 -16.2 -17.5 -20.5 -18.9 -21.1 -23.6 -26.1 -29.4 -28.3];
delay_profile_rma_los_us = [0.005 0.02 0.02 0.03 0.045 0.065 0.065 0.09 0.125 0.180 0.190];

% D1-Rural Macro (RMa) NLOS (Winner II)
power_profile_rma_nlos_db = [-5.2 -1.8 -3.3 -7.0 -5.3 -7.1 -9.0 -4.2 -12.4 -26.5];
delay_profile_rma_nlos_us = [0.005 0.000 0.005 0.015 0.02 0.025 0.055 0.100 0.170 0.420];

% B1-Urban Micro (UMi) LOS (Winner II)
power_profile_umi_los_db = [0 -12.7 -14.8 -15.8 -13.9 -17.8 -19.6 -31.4];
delay_profile_umi_los_us = [0 0.035 0.055 0.065 0.105 0.115 0.25 0.46];

% B1-Urban Micro (UMi) NLOS (Winner II)
power_profile_umi_nlos_db = [-1.0 -5.2 -6.1 -8.1 -8.6 -11.7 -12.0 -12.9 -19.6 -23.9 -22.1 -25.6 -23.3 -32.3 -31.7 -29.9];
delay_profile_umi_nlos_us = [0 0.095 0.105 0.115 0.230 0.240 0.245 0.285 0.390 0.430 0.460 0.505 0.515 0.595 0.600 0.615];

% C1-Urban macro-cell LOS (Winner II)
power_profile_uma_los_db = [-25.3 -21.6 -26.3 -25.1 -25.4 -22.0 -29.2 -26.5 -23.2 -32.2 -26.5 -32.1 -28.5 -30.5 -32.6];
delay_profile_uma_los_us = [0.005 0.085 0.135 0.135 0.170 0.190 0.275 0.295 0.290 0.410 0.445 0.500 0.620 0.655 0.960];

% C1-Urban macro-cell NLOS (Winner II)
power_profile_uma_nlos_db = [-5.2 -7.5 -10.5 -3.2 -8.3 -14.0 -6.4 -3.1 -4.6 -8.0 -7.2 -3.1 -9.5 -22.4];
delay_profile_uma_nlos_us = [0.005 0.025 0.035 0.035 0.050 0.065 0.065 0.075 0.145 0.160 0.195 0.200 0.205 0.770];

% Extended Vehicular A (EVA) for max. doppler 70 Hz
power_profile_eva_db = [0.0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];
delay_profile_eva_us = [0.0 0.03 0.15 0.31 0.37 0.71 1.09 1.73 2.51];

% Extended Typical Urban (ETU) for max. doppler 300 Hz
power_profile_etu_db = [-1.0 -1.0 -1.0 0.0 0.0 0.0 -3.0 -5.0 -7.0];
delay_profile_etu_us = [0.0 0.05 0.12 0.2 0.23 0.5 1.6 2.3 5.0];

% tdl-a
power_profile_tdl_a_db = [-13.4 0 -2.2 -4 -6 -8.2 -9.9 -10.5 -7.5 -15.9 -6.6 -16.7 -12.4 -15.2 -10.8 -11.3 -12.7 -16.2 -18.3 -18.9 -16.6 -19.9 -29.7];
delay_profile_tdl_a_norm = [0.0000 0.3819 0.4025 0.5868 0.4610 0.5375 0.6708 0.5750 0.7618 1.5375 1.8978 2.2242 2.1718 2.4942 2.5119 3.0582 4.0810 4.4579 4.5695 4.7966 5.0066 5.3043 9.6586];
delay_profile_tdl_a_us = delay_profile_tdl_a_norm * delay_spread_rms_us;

% tdl-b
power_profile_tdl_b_db = [0 -2.2 -4 -3.2 -9.8 -1.2 -3.4 -5.2 -7.6 -3 -8.9 -9 -4.8 -5.7 -7.5 -1.9 -7.6 -12.2 -9.8 -11.4 -14.9 -9.2 -11.3];
delay_profile_tdl_b_norm = [0.0000 0.1072 0.2155 0.2095 0.2870 0.2986 0.3752 0.5055 0.3681 0.3697 0.5700 0.5283 1.1021 1.2756 1.5474 1.7842 2.0169 2.8294 3.0219 3.6187 4.1067 4.2790 4.7834];
delay_profile_tdl_b_us = delay_profile_tdl_b_norm * delay_spread_rms_us;

% tdl-c
power_profile_tdl_c_db = [-4.4 -1.2 -3.5 -5.2 -2.5 0 -2.2 -3.9 -7.4 -7.1 -10.7 -11.1 -5.1 -6.8 -8.7 -13.2 -13.9 -13.9 -15.8 -17.1 -16 -15.7 -21.6 -22.8];
delay_profile_tdl_c_norm = [0 0.2099 0.2219 0.2329 0.2176 0.6366 0.6448 0.6560 0.6584 0.7935 0.8213 0.9336 1.2285 1.3083 2.1704 2.7105 4.2589 4.6003 5.4902 5.6077 6.3065 6.6374 7.0427 8.6523];
delay_profile_tdl_c_us = delay_profile_tdl_c_norm * delay_spread_rms_us;

% tdl-d
% power_profile_tdl_d_los_db = -0.2;
% delay_profile_tdl_d_los_norm = 0;
% delay_profile_tdl_d_los_us = delay_profile_tdl_d_los_norm * delay_spread_rms;
power_profile_tdl_d_db = [-13.5 -18.8 -21 -22.8 -17.9 -20.1 -21.9 -22.9 -27.8 -23.6 -24.8 -30.0 -27.7];
delay_profile_tdl_d_norm = [0 0.035 0.612 1.363 1.405 1.804 2.596 1.775 4.042 7.937 9.424 9.708 12.525];
delay_profile_tdl_d_us = delay_profile_tdl_d_norm * delay_spread_rms_us;
k_factor_tdl_d_db = 13.3;

% tdl-e
% power_profile_tdl_e_los_db = -0.03;
% delay_profile_tdl_e_los_norm = 0;
% delay_profile_tdl_e_los_us = delay_profile_tdl_e_los_norm * delay_spread_rms;
power_profile_tdl_e_db = [-22.03 -15.8 -18.1 -19.8 -22.9 -22.4 -18.6 -20.8 -22.6 -22.3 -25.6 -20.2 -29.8 -29.2];
delay_profile_tdl_e_norm = [0 0.5133 0.5440 0.5630 0.5440 0.7112 1.9092 1.9293 1.9589 2.6426 3.7136 5.4524 12.0034 20.6519];
delay_profile_tdl_e_us = delay_profile_tdl_e_norm * delay_spread_rms_us;
k_factor_tdl_e_db = 22;

% test channel
% power_profile_test_db = [-3.02 0 -2 -6 -7.93 -10.1];
% delay_profile_test_us = [0 0.04 0.1 0.32 0.46 1];
power_profile_test_db = 0;
delay_profile_test_us = 0;

switch idx_fading
    case 1
        power_profile_db = power_profile_rma_los_db;
        delay_profile_us = delay_profile_rma_los_us;
        los = true;
        k_factor_db = 0;
    case 2
        power_profile_db = power_profile_rma_nlos_db;
        delay_profile_us = delay_profile_rma_nlos_us;
        los = false;
    case 3
        power_profile_db = power_profile_umi_los_db;
        delay_profile_us = delay_profile_umi_los_us;
        los = true;
        k_factor_db = 0;
    case 4
        power_profile_db = power_profile_umi_nlos_db;
        delay_profile_us = delay_profile_umi_nlos_us;
        los = false;
    case 5
        power_profile_db = power_profile_uma_los_db;
        delay_profile_us = delay_profile_uma_los_us;
        los = true;
        k_factor_db = 0;
    case 6
        power_profile_db = power_profile_uma_nlos_db;
        delay_profile_us = delay_profile_uma_nlos_us;
        los = false;
    case 7
        power_profile_db = power_profile_eva_db;
        delay_profile_us = delay_profile_eva_us;
        los = false;
    case 8
        power_profile_db = power_profile_etu_db;
        delay_profile_us = delay_profile_etu_us;
        los = false;
    case 9
        power_profile_db = power_profile_tdl_a_db;
        delay_profile_us = delay_profile_tdl_a_us;
        los = false;
    case 10
        power_profile_db = power_profile_tdl_b_db;
        delay_profile_us = delay_profile_tdl_b_us;
        los = false;
    case 11
        power_profile_db = power_profile_tdl_c_db;
        delay_profile_us = delay_profile_tdl_c_us;
        los = false;
    case 12
        power_profile_db = power_profile_tdl_d_db;
        delay_profile_us = delay_profile_tdl_d_us;
        los = true;
        k_factor_db = k_factor_tdl_d_db;
    case 13
        power_profile_db = power_profile_tdl_e_db;
        delay_profile_us = delay_profile_tdl_e_us;
        los = true;
        k_factor_db = k_factor_tdl_e_db;
    otherwise
        power_profile_db = power_profile_test_db;
        delay_profile_us = delay_profile_test_us;
        los = false;
end

% channel parameters
delay_profile = delay_profile_us * 1e-6;
doppler_freq = velocity_kmh * (1000 / 3600) / (3e8 / (carrier_freq_mhz * 1e6));
t_coherence = 0.423 / doppler_freq;                                         % coherence time (sec)
power_profile = 10 .^ (power_profile_db / 10);
t_excess = (power_profile * delay_profile') / sum(power_profile);           % mean excess delay (sec)
t_second = (power_profile * (delay_profile .^ 2)') / sum(power_profile);    % second moment (sec^2)
t_rms_delay = sqrt(t_second - t_excess ^ 2);                                % rms delay spread (sec)
coherence_bw = 1 / (5 * t_rms_delay);                                       % coherence bandwidth (khz)

% output
nw_ch.average_path_gains = power_profile_db;
nw_ch.path_delays = delay_profile;
nw_ch.maximum_doppler_shift = doppler_freq;
nw_ch.t_coherence = t_coherence;
nw_ch.coherence_bw = coherence_bw;
nw_ch.ebno_db = ebno_db;
nw_ch.doppler_spectrum = doppler('Jakes');   % doppler model (create a doppler spectrum structure)
nw_ch.los = los;
if los
    nw_ch.k_factor = 10^(0.1*k_factor_db);
    nw_ch.delay_spread_rms = delay_spread_rms_us * 1e-6;
end

end
