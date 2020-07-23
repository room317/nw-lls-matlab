% simulator: new_wave_lls
% commit: 

% RUNNING CASE: OTFS BW10 FC4000 V120 DS0.1 TDLA MCS8 LEN6120 SIM10000 SNR0 DDTONE TFEQMMSE
% 9 symbol received
per_case01 = [ ...
     0.0    1.000000
     2.0    1.000000
     4.0    1.000000
     6.0    1.000000
     8.0    1.000000
    10.0    1.000000
    12.0    1.000000
    14.0    1.000000
    16.0    1.000000
    18.0    1.000000];

% 10 symbols received
per_case02 = [ ...
     0.0    1.000000
     2.0    1.000000
     4.0    1.000000
     6.0    1.000000
     8.0    1.000000
    10.0    1.000000
    12.0    1.000000
    14.0    1.000000
    16.0    0.997512
    18.0    0.968599];

% 11 symbol received
per_case03 = [ ...
     0.0    1.000000
     2.0    1.000000
     4.0    1.000000
     6.0    1.000000
     8.0    1.000000
    10.0    0.909297
    12.0    0.704745
    14.0    0.435396
    16.0    0.179982
    18.0    0.041985];

% 12 symbol received
per_case04 = [ ...
     0.0    1.000000
     2.0    1.000000
     4.0    1.000000
     6.0    1.000000
     8.0    0.806841
    10.0    0.477950
    12.0    0.170348
    14.0    0.024400
    16.0    0.001900
    18.0    0.000200];

% 13 symbol received
per_case05 = [ ...
     0.0    1.000000
     2.0    1.000000
     4.0    1.000000
     6.0    0.851380
     8.0    0.457763
    10.0    0.117149
    12.0    0.011400
    14.0    0.000500];

% 14 symbol received
per_case06 = [ ...
     0.0    1.000000
     2.0    1.000000
     4.0    0.966265
     6.0    0.612214
     8.0    0.168064
    10.0    0.017300
    12.0    0.000200];

figure
semilogy(per_case02(:,1),per_case02(:,2),'-ko'), hold on
semilogy(per_case03(:,1),per_case03(:,2),'-kx')
semilogy(per_case04(:,1),per_case04(:,2),'-ks')
semilogy(per_case05(:,1),per_case05(:,2),'-k^')
semilogy(per_case06(:,1),per_case06(:,2),'-kd'), hold off
grid minor
title({'PER performance of partially received subframe'; 'OTFS BW10 FC4000 DS0.1 TDLA MCS8 LEN6120 DDTONE TFEQMMSE'})
xlabel('SNR (dB)'), ylabel('PER')
legend('10/14 received', '11/14 received', '12/14 received', '13/14 received', '14/14 received')
axis([-1 19 1e-3 1])
