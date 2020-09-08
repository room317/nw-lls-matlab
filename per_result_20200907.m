% simulator: new_wave_test_r3
% commit: a48c3a6
% maintain bandwidth (10 MHz), single slot (14 ofdm symbols) simulation

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V3 DS0.1 TDLC MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case01 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        447   0.897092
     6.0        893   0.449048
     8.0       1707   0.234915
    10.0       3742   0.107162
    12.0      10000   0.038300
    14.0      10000   0.011500
    16.0      10000   0.002700
    18.0      10000   0.000300];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V120 DS0.1 TDLC MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case02 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        469   0.855011
     6.0        896   0.447545
     8.0       1640   0.244512
    10.0       3925   0.102166
    12.0      10000   0.034700
    14.0      10000   0.006000
    16.0      10000   0.001000
    18.0      10000   0.000100];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V200 DS0.1 TDLC MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case03 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        446   0.899103
     6.0        767   0.522816
     8.0       1497   0.267869
    10.0       4219   0.095046
    12.0      10000   0.025000
    14.0      10000   0.003500
    16.0      10000   0.000400];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V500 DS0.1 TDLC MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case04 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        414   0.968599
     6.0        667   0.601199
     8.0       1376   0.291424
    10.0       5994   0.066900
    12.0      10000   0.010600
    14.0      10000   0.000800];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V3 DS0.1 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case05 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        405   0.990123
     6.0        641   0.625585
     8.0       1537   0.260898
    10.0       5931   0.067611
    12.0      10000   0.014700
    14.0      10000   0.002900
    16.0      10000   0.000200];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V120 DS0.1 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case06 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        412   0.973301
     6.0        630   0.636508
     8.0       1646   0.243621
    10.0       6563   0.061100
    12.0      10000   0.010400
    14.0      10000   0.001600
    16.0      10000   0.000100];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V200 DS0.1 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case07 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        402   0.997512
     6.0        615   0.652033
     8.0       1662   0.241276
    10.0       8214   0.048819
    12.0      10000   0.007500
    14.0      10000   0.000900];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V500 DS0.1 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case08 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        401   1.000000
     6.0        526   0.762357
     8.0       1594   0.251568
    10.0      10000   0.032800
    12.0      10000   0.002300
    14.0      10000   0.000100];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V3 DS0.1 TDLC MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case09 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        467   0.858672
     6.0        903   0.444075
     8.0       1746   0.229668
    10.0       3936   0.101880
    12.0      10000   0.036500
    14.0      10000   0.011400
    16.0      10000   0.001800
    18.0      10000   0.000700];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V120 DS0.1 TDLC MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case10 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        457   0.877462
     6.0        889   0.451069
     8.0       1681   0.238548
    10.0       3960   0.101263
    12.0      10000   0.030900
    14.0      10000   0.005100
    16.0      10000   0.000700];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V200 DS0.1 TDLC MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case11 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        442   0.907240
     6.0        791   0.506953
     8.0       1537   0.260898
    10.0       4568   0.087785
    12.0      10000   0.021700
    14.0      10000   0.002100
    16.0      10000   0.000300];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V500 DS0.1 TDLC MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case12 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        405   0.990123
     6.0        613   0.654160
     8.0       1517   0.264338
    10.0       7038   0.056976
    12.0      10000   0.007000
    14.0      10000   0.000700];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V3 DS0.1 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case13 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        411   0.975669
     6.0        638   0.628527
     8.0       1675   0.239403
    10.0       5728   0.070007
    12.0      10000   0.011700
    14.0      10000   0.002900
    16.0      10000   0.000200];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V120 DS0.1 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case14 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        406   0.987685
     6.0        663   0.604827
     8.0       1498   0.267690
    10.0       6698   0.059869
    12.0      10000   0.009000
    14.0      10000   0.001100
    16.0      10000   0.000000];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V200 DS0.1 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case15 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        403   0.995037
     6.0        572   0.701049
     8.0       1811   0.221425
    10.0       7499   0.053474
    12.0      10000   0.006700
    14.0      10000   0.000300];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V500 DS0.1 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case16 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        402   0.997512
     6.0        520   0.771154
     8.0       1544   0.259715
    10.0      10000   0.029200
    12.0      10000   0.000800];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V3 DS0.01 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case17 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        605   0.662810
     6.0       1235   0.324696
     8.0       2180   0.183945
    10.0       4034   0.099405
    12.0       5852   0.068524
    14.0      10000   0.033400
    16.0      10000   0.016800
    18.0      10000   0.006700];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V120 DS0.01 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case18 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        540   0.742593
     6.0       1268   0.316246
     8.0       2165   0.185219
    10.0       3310   0.121148
    12.0       5741   0.069848
    14.0      10000   0.030300
    16.0      10000   0.012000
    18.0      10000   0.004100];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V200 DS0.01 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case19 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        521   0.769674
     6.0       1077   0.372331
     8.0       1866   0.214898
    10.0       3482   0.115164
    12.0       6530   0.061409
    14.0      10000   0.023700
    16.0      10000   0.006500
    18.0      10000   0.001500];

% RUNNING CASE: OFDM BW10 SCS60 FC4000 V500 DS0.01 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case20 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        436   0.919725
     6.0        801   0.500624
     8.0       1426   0.281206
    10.0       3525   0.113759
    12.0      10000   0.033400
    14.0      10000   0.007600
    16.0      10000   0.002100
    18.0      10000   0.000100];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V3 DS0.01 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case21 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        589   0.680815
     6.0       1385   0.289531
     8.0       2428   0.165157
    10.0       3464   0.115762
    12.0       6957   0.057640
    14.0      10000   0.032100
    16.0      10000   0.017600
    18.0      10000   0.008300];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V120 DS0.01 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case22 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        557   0.719928
     6.0       1217   0.329499
     8.0       2141   0.187296
    10.0       3692   0.108613
    12.0       6698   0.059869
    14.0      10000   0.024800
    16.0      10000   0.012200
    18.0      10000   0.003100];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V200 DS0.01 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case23 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        538   0.745353
     6.0       1097   0.365542
     8.0       1923   0.208528
    10.0       3453   0.116131
    12.0       8523   0.047049
    14.0      10000   0.018600
    16.0      10000   0.005900
    18.0      10000   0.001400];

% RUNNING CASE: OTFS BW10 SCS60 FC4000 V500 DS0.01 TDLE MCS8 LEN3537 SIM10000 SNR0 REAL TFEQMMSE
per_case24 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        443   0.905192
     6.0        819   0.489621
     8.0       1592   0.251884
    10.0       4193   0.095636
    12.0      10000   0.027700
    14.0      10000   0.006200
    16.0      10000   0.001400
    18.0      10000   0.000100];

% figure
% semilogy(per_case01(:,1),per_case01(:,3),'-ko'), hold on
% semilogy(per_case02(:,1),per_case02(:,3),'-kx')
% semilogy(per_case03(:,1),per_case03(:,3),'-ks')
% semilogy(per_case04(:,1),per_case04(:,3),'-k^'), hold off
% grid minor
% title('OFDM BW10 SCS60 FC4000 DS0.1 TDLC MCS8 LEN3537 REAL TFEQMMSE')
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case05(:,1),per_case05(:,3),'-ko'), hold on
% semilogy(per_case06(:,1),per_case06(:,3),'-kx')
% semilogy(per_case07(:,1),per_case07(:,3),'-ks')
% semilogy(per_case08(:,1),per_case08(:,3),'-k^'), hold off
% grid minor
% title('OFDM BW10 SCS60 FC4000 DS0.1 TDLE MCS8 LEN3537 REAL TFEQMMSE')
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case09(:,1),per_case09(:,3),'-ko'), hold on
% semilogy(per_case10(:,1),per_case10(:,3),'-kx')
% semilogy(per_case11(:,1),per_case11(:,3),'-ks')
% semilogy(per_case12(:,1),per_case12(:,3),'-k^'), hold off
% grid minor
% title('OTFS BW10 SCS60 FC4000 DS0.1 TDLC MCS8 LEN3537 REAL TFEQMMSE')
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case13(:,1),per_case13(:,3),'-ko'), hold on
% semilogy(per_case14(:,1),per_case14(:,3),'-kx')
% semilogy(per_case15(:,1),per_case15(:,3),'-ks')
% semilogy(per_case16(:,1),per_case16(:,3),'-k^'), hold off
% grid minor
% title('OTFS BW10 SCS60 FC4000 DS0.1 TDLE MCS8 LEN3537 REAL TFEQMMSE')
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case17(:,1),per_case17(:,3),'-ko'), hold on
% semilogy(per_case18(:,1),per_case18(:,3),'-kx')
% semilogy(per_case19(:,1),per_case19(:,3),'-ks')
% semilogy(per_case20(:,1),per_case20(:,3),'-k^'), hold off
% grid minor
% title('OFDM BW10 SCS60 FC4000 DS0.01 TDLE MCS8 LEN3537 REAL TFEQMMSE')
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case21(:,1),per_case21(:,3),'-ko'), hold on
% semilogy(per_case22(:,1),per_case22(:,3),'-kx')
% semilogy(per_case23(:,1),per_case23(:,3),'-ks')
% semilogy(per_case24(:,1),per_case24(:,3),'-k^'), hold off
% grid minor
% title('OTFS BW10 SCS60 FC4000 DS0.1 TDLE MCS8 LEN3537 REAL TFEQMMSE')
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])

figure
semilogy(per_case01(:,1),per_case01(:,3),'-bo'), hold on
semilogy(per_case02(:,1),per_case02(:,3),'-bx')
semilogy(per_case03(:,1),per_case03(:,3),'-bs')
semilogy(per_case04(:,1),per_case04(:,3),'-b^')
semilogy(per_case09(:,1),per_case09(:,3),'-ro')
semilogy(per_case10(:,1),per_case10(:,3),'-rx')
semilogy(per_case11(:,1),per_case11(:,3),'-rs')
semilogy(per_case12(:,1),per_case12(:,3),'-r^'), hold off
grid minor
title('BW10 SCS60 FC4000 DS0.1 TDLC MCS8 LEN3537 REAL TFEQMMSE')
xlabel('SNR (dB)'), ylabel('PER')
legend('OFDM V3', 'OFDM V120', 'OFDM V200', 'OFDM V500', 'OTFS V3', 'OTFS V120', 'OTFS V200', 'OTFS V500')
axis([-1 19 1e-3 1])

figure
semilogy(per_case05(:,1),per_case05(:,3),'-bo'), hold on
semilogy(per_case06(:,1),per_case06(:,3),'-bx')
semilogy(per_case07(:,1),per_case07(:,3),'-bs')
semilogy(per_case08(:,1),per_case08(:,3),'-b^')
semilogy(per_case13(:,1),per_case13(:,3),'-ro')
semilogy(per_case14(:,1),per_case14(:,3),'-rx')
semilogy(per_case15(:,1),per_case15(:,3),'-rs')
semilogy(per_case16(:,1),per_case16(:,3),'-r^'), hold off
grid minor
title('BW10 SCS60 FC4000 DS0.1 TDLE MCS8 LEN3537 REAL TFEQMMSE')
xlabel('SNR (dB)'), ylabel('PER')
legend('OFDM V3', 'OFDM V120', 'OFDM V200', 'OFDM V500', 'OTFS V3', 'OTFS V120', 'OTFS V200', 'OTFS V500')
axis([-1 19 1e-3 1])

figure
semilogy(per_case17(:,1),per_case17(:,3),'-bo'), hold on
semilogy(per_case18(:,1),per_case18(:,3),'-bx')
semilogy(per_case19(:,1),per_case19(:,3),'-bs')
semilogy(per_case20(:,1),per_case20(:,3),'-b^')
semilogy(per_case21(:,1),per_case21(:,3),'-ro')
semilogy(per_case22(:,1),per_case22(:,3),'-rx')
semilogy(per_case23(:,1),per_case23(:,3),'-rs')
semilogy(per_case24(:,1),per_case24(:,3),'-r^'), hold off
grid minor
title('BW10 SCS60 FC4000 DS0.01 TDLE MCS8 LEN3537 REAL TFEQMMSE')
xlabel('SNR (dB)'), ylabel('PER')
legend('OFDM V3', 'OFDM V120', 'OFDM V200', 'OFDM V500', 'OTFS V3', 'OTFS V120', 'OTFS V200', 'OTFS V500')
axis([-1 19 1e-3 1])
