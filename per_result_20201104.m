% simulator: new_wave_lls
% commit: 
% bandwidth 5 MHz, subcarrier spacing 60 khz, 4 slots, single user, 1 rbs x 4 slot per user, code rate 2/3
% comment: to compare inter-simulator per performance
%   1. awgn vs. tdl-c 100 ns vs. tdl-e 10 ns
%   2. ofdm

%% ofdm / single user / short transmit block / partial spreading

% RUNNING CASE: OFDM MU BW5 SCS60 SLOT4 FC4000 AWGN MCS8 LEN1748 SIM10000 SNR8 REAL TFEQMMSE
per_case01 = [ ...
     8.0        401   1.000000
     8.2        406   0.987685
     8.4        439   0.913440
     8.6        602   0.666113
     8.8       1206   0.332504
     9.0       3243   0.123651
     9.2      10000   0.031900
     9.4      10000   0.007900
     9.6      10000   0.002400];

per_case02 = [ ...
    10.0        100   1.000000
    12.0        100   0.730000
    14.0        100   0.400000
    16.0        100   0.160000
    18.0        100   0.070000
    20.0        100   0.000000];

per_case03 = [ ...
    10.0        100   1.000000
    12.0        100   0.820000
    14.0        100   0.430000
    16.0        100   0.210000
    18.0        100   0.050000
    20.0        100   0.000000];

% REF: BW 5 MHz TDL-C 100 ns 16 QAM CR 2/3
per_case11 = [ ...
    12.0      10000   0.5024
    14.0      10000   0.194
    16.0      10000   0.06463
    18.0      10000   0.01598
    20.0      10000   0.003117
    22.0      10000   0.0004745];

% REF: BW 5 MHz TDL-E 10 ns 16 QAM CR 2/3
per_case12 = [ ...
    10.0      10000   0.4463
    11.0      10000   0.1473
    12.0      10000   0.03149
    13.0      10000   0.004296
    14.0      10000   0.0003964];

% basic graphs

figure
semilogy(per_case01(:,1),per_case01(:,3),'-ko'), hold on
semilogy(per_case02(:,1),per_case02(:,3),'-kx')
semilogy(per_case03(:,1),per_case03(:,3),'-ks')
semilogy(per_case11(:,1),per_case11(:,3),'-rx')
semilogy(per_case12(:,1),per_case12(:,3),'-rs'), hold off
grid minor
title({'OFDM MU BW5 SCS60 SLOT4 FC4000 16QAM R2/3 LEN1748 REAL TFEQMMSE';'SINGLE USER (1x4)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('AWGN', 'TDL-C 100 ns', 'TDL-E 10 ns', 'REF: TDL-C 100 ns', 'REF: TDL-E 10 ns')
axis([8 22 1e-3 1])

% figure
% semilogy(per_case05(:,1),per_case05(:,3),'-ko'), hold on
% semilogy(per_case06(:,1),per_case06(:,3),'-kx')
% semilogy(per_case07(:,1),per_case07(:,3),'-ks')
% semilogy(per_case08(:,1),per_case08(:,3),'-k^'), hold off
% grid minor
% title({'OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN964 REAL TFEQMMSE';'SINGLE USER (3x1)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case09(:,1),per_case09(:,3),'-ko'), hold on
% semilogy(per_case10(:,1),per_case10(:,3),'-kx')
% semilogy(per_case11(:,1),per_case11(:,3),'-ks')
% semilogy(per_case12(:,1),per_case12(:,3),'-k^')
% semilogy(per_case09(:,1),per_case09(:,4),'-bo')
% semilogy(per_case10(:,1),per_case10(:,4),'-bx')
% semilogy(per_case11(:,1),per_case11(:,4),'-bs')
% semilogy(per_case12(:,1),per_case12(:,4),'-b^')
% semilogy(per_case09(:,1),per_case09(:,5),'-ro')
% semilogy(per_case10(:,1),per_case10(:,5),'-rx')
% semilogy(per_case11(:,1),per_case11(:,5),'-rs')
% semilogy(per_case12(:,1),per_case12(:,5),'-r^')
% semilogy(per_case09(:,1),per_case09(:,6),'-go')
% semilogy(per_case10(:,1),per_case10(:,6),'-gx')
% semilogy(per_case11(:,1),per_case11(:,6),'-gs')
% semilogy(per_case12(:,1),per_case12(:,6),'-g^'), hold off
% grid minor
% title({'OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN964 REAL TFEQMMSE'; '4 USERS (3x1)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case13(:,1),per_case13(:,3),'-ko'), hold on
% semilogy(per_case14(:,1),per_case14(:,3),'-kx')
% semilogy(per_case15(:,1),per_case15(:,3),'-ks')
% semilogy(per_case16(:,1),per_case16(:,3),'-k^')
% semilogy(per_case13(:,1),per_case13(:,4),'-bo')
% semilogy(per_case14(:,1),per_case14(:,4),'-bx')
% semilogy(per_case15(:,1),per_case15(:,4),'-bs')
% semilogy(per_case16(:,1),per_case16(:,4),'-b^')
% semilogy(per_case13(:,1),per_case13(:,5),'-ro')
% semilogy(per_case14(:,1),per_case14(:,5),'-rx')
% semilogy(per_case15(:,1),per_case15(:,5),'-rs')
% semilogy(per_case16(:,1),per_case16(:,5),'-r^')
% semilogy(per_case13(:,1),per_case13(:,6),'-go')
% semilogy(per_case14(:,1),per_case14(:,6),'-gx')
% semilogy(per_case15(:,1),per_case15(:,6),'-gs')
% semilogy(per_case16(:,1),per_case16(:,6),'-g^'), hold off
% grid minor
% title({'OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN964 REAL TFEQMMSE'; '4 USERS (3x1)'})
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
% title({'OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE'; 'SINGLE USER (3x1)'})
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
% title({'OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE'; 'SINGLE USER (3x1)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case25(:,1),per_case25(:,3),'-ko'), hold on
% semilogy(per_case26(:,1),per_case26(:,3),'-kx')
% semilogy(per_case27(:,1),per_case27(:,3),'-ks')
% semilogy(per_case28(:,1),per_case28(:,3),'-k^'), hold off
% grid minor
% title({'OFDM BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE'; 'SINGLE USER (3x1)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case29(:,1),per_case29(:,3),'-ko'), hold on
% semilogy(per_case30(:,1),per_case30(:,3),'-kx')
% semilogy(per_case31(:,1),per_case31(:,3),'-ks')
% semilogy(per_case32(:,1),per_case32(:,3),'-k^'), hold off
% grid minor
% title({'OFDM BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE'; 'SINGLE USER (3x1)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case33(:,1),per_case33(:,3),'-ko'), hold on
% semilogy(per_case34(:,1),per_case34(:,3),'-kx')
% semilogy(per_case35(:,1),per_case35(:,3),'-ks')
% semilogy(per_case36(:,1),per_case36(:,3),'-k^'), hold off
% grid minor
% title({'OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE'; 'SINGLE USER (FULL)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case37(:,1),per_case37(:,3),'-ko'), hold on
% semilogy(per_case38(:,1),per_case38(:,3),'-kx')
% semilogy(per_case39(:,1),per_case39(:,3),'-ks')
% semilogy(per_case40(:,1),per_case40(:,3),'-k^'), hold off
% grid minor
% title({'OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE'; 'SINGLE USER (FULL)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3', 'V120', 'V200', 'V500')
% axis([-1 19 1e-3 1])

% more graphs
% 
% figure
% semilogy(per_case17(:,1),per_case17(:,3),'-ko'), hold on
% semilogy(per_case09(:,1),per_case09(:,3),':ko')
% semilogy(per_case09(:,1),per_case09(:,4),':bo')
% semilogy(per_case09(:,1),per_case09(:,5),':ro')
% semilogy(per_case09(:,1),per_case09(:,6),':go'), hold off
% grid minor
% title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, LONG CB, 3 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case20(:,1),per_case20(:,3),'-ko'), hold on
% semilogy(per_case12(:,1),per_case12(:,3),':kx')
% semilogy(per_case12(:,1),per_case12(:,4),':bx')
% semilogy(per_case12(:,1),per_case12(:,5),':rx')
% semilogy(per_case12(:,1),per_case12(:,6),':gx'), hold off
% grid minor
% title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, LONG CB, 500 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case21(:,1),per_case21(:,3),'-ko'), hold on
% semilogy(per_case13(:,1),per_case13(:,3),':ko')
% semilogy(per_case13(:,1),per_case13(:,4),':bo')
% semilogy(per_case13(:,1),per_case13(:,5),':ro')
% semilogy(per_case13(:,1),per_case13(:,6),':go'), hold off
% grid minor
% title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, LONG CB, 3 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case24(:,1),per_case24(:,3),'-ko'), hold on
% semilogy(per_case16(:,1),per_case16(:,3),':kx')
% semilogy(per_case16(:,1),per_case16(:,4),':bx')
% semilogy(per_case16(:,1),per_case16(:,5),':rx')
% semilogy(per_case16(:,1),per_case16(:,6),':gx'), hold off
% grid minor
% title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, LONG CB, 500 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case01(:,1),per_case01(:,3),'-ko'), hold on
% semilogy(per_case02(:,1),per_case02(:,3),'-kx')
% semilogy(per_case03(:,1),per_case03(:,3),'-ks')
% semilogy(per_case04(:,1),per_case04(:,3),'-k^')
% semilogy(per_case17(:,1),per_case17(:,3),'-ro')
% semilogy(per_case18(:,1),per_case18(:,3),'-rx')
% semilogy(per_case19(:,1),per_case19(:,3),'-rs')
% semilogy(per_case20(:,1),per_case20(:,3),'-r^'), hold off
% grid minor
% title({'SHORT CODE BLOCK VS LONG CODE BLOCK, RESOURCE: 3x1, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3 (SHORT CB)', 'V120 (SHORT CB)', 'V200 (SHORT CB)', 'V500 (SHORT CB)', 'V3 (LONG CB)', 'V120 (LONG CB)', 'V200 (LONG CB)', 'V500 (LONG CB)')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case05(:,1),per_case05(:,3),'-ko'), hold on
% semilogy(per_case06(:,1),per_case06(:,3),'-kx')
% semilogy(per_case07(:,1),per_case07(:,3),'-ks')
% semilogy(per_case08(:,1),per_case08(:,3),'-k^')
% semilogy(per_case21(:,1),per_case21(:,3),'-ro')
% semilogy(per_case22(:,1),per_case22(:,3),'-rx')
% semilogy(per_case23(:,1),per_case23(:,3),'-rs')
% semilogy(per_case24(:,1),per_case24(:,3),'-r^'), hold off
% grid minor
% title({'SHORT CODE BLOCK VS LONG CODE BLOCK, RESOURCE: 3x1, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('V3 (SHORT CB)', 'V120 (SHORT CB)', 'V200 (SHORT CB)', 'V500 (SHORT CB)', 'V3 (LONG CB)', 'V120 (LONG CB)', 'V200 (LONG CB)', 'V500 (LONG CB)')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case17(:,1),per_case17(:,3),'-bo'), hold on
% semilogy(per_case18(:,1),per_case18(:,3),'-bx')
% semilogy(per_case19(:,1),per_case19(:,3),'-bs')
% semilogy(per_case20(:,1),per_case20(:,3),'-b^')
% semilogy(per_case25(:,1),per_case25(:,3),'-ro')
% semilogy(per_case26(:,1),per_case26(:,3),'-rx')
% semilogy(per_case27(:,1),per_case27(:,3),'-rs')
% semilogy(per_case28(:,1),per_case28(:,3),'-r^'), hold off
% grid minor
% title({'OFDM VS OTFS, RESOURCE: 3x1, LONG CB, SINGLE USER';'(BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('OTFS V3', 'OTFS V120', 'OTFS V200', 'OTFS V500', 'OFDM V3', 'OFDM V120', 'OFDM V200', 'OFDM V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case21(:,1),per_case21(:,3),'-bo'), hold on
% semilogy(per_case22(:,1),per_case22(:,3),'-bx')
% semilogy(per_case23(:,1),per_case23(:,3),'-bs')
% semilogy(per_case24(:,1),per_case24(:,3),'-b^')
% semilogy(per_case29(:,1),per_case29(:,3),'-ro')
% semilogy(per_case30(:,1),per_case30(:,3),'-rx')
% semilogy(per_case31(:,1),per_case31(:,3),'-rs')
% semilogy(per_case32(:,1),per_case32(:,3),'-r^'), hold off
% grid minor
% title({'OFDM VS OTFS, RESOURCE: 3x1, LONG CB, SINGLE USER';'(BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('OTFS V3', 'OTFS V120', 'OTFS V200', 'OTFS V500', 'OFDM V3', 'OFDM V120', 'OFDM V200', 'OFDM V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case17(:,1),per_case17(:,3),'-bo'), hold on
% semilogy(per_case18(:,1),per_case18(:,3),'-bx')
% semilogy(per_case19(:,1),per_case19(:,3),'-bs')
% semilogy(per_case20(:,1),per_case20(:,3),'-b^')
% semilogy(per_case33(:,1),per_case33(:,3),'-ro')
% semilogy(per_case34(:,1),per_case34(:,3),'-rx')
% semilogy(per_case35(:,1),per_case35(:,3),'-rs')
% semilogy(per_case36(:,1),per_case36(:,3),'-r^'), hold off
% grid minor
% title({'FULL SPREAD(11x4) VS PARTIAL SPREAD(3x1), LONG CB, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('PARTIAL V3', 'PARTIAL V120', 'PARTIAL V200', 'PARTIAL V500', 'FULL V3', 'FULL V120', 'FULL V200', 'FULL V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case21(:,1),per_case21(:,3),'-bo'), hold on
% semilogy(per_case22(:,1),per_case22(:,3),'-bx')
% semilogy(per_case23(:,1),per_case23(:,3),'-bs')
% semilogy(per_case24(:,1),per_case24(:,3),'-b^')
% semilogy(per_case37(:,1),per_case37(:,3),'-ro')
% semilogy(per_case38(:,1),per_case38(:,3),'-rx')
% semilogy(per_case39(:,1),per_case39(:,3),'-rs')
% semilogy(per_case40(:,1),per_case40(:,3),'-r^'), hold off
% grid minor
% title({'FULL SPREAD(11x4) VS PARTIAL SPREAD(3x1), LONG CB, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('PARTIAL V3', 'PARTIAL V120', 'PARTIAL V200', 'PARTIAL V500', 'FULL V3', 'FULL V120', 'FULL V200', 'FULL V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case17(:,1),per_case17(:,3),'-bo'), hold on
% semilogy(per_case18(:,1),per_case18(:,3),'-bx')
% semilogy(per_case19(:,1),per_case19(:,3),'-bs')
% semilogy(per_case20(:,1),per_case20(:,3),'-b^')
% semilogy(per_case41(:,1),per_case41(:,3),'-ro')
% semilogy(per_case42(:,1),per_case42(:,3),'-rx')
% semilogy(per_case43(:,1),per_case43(:,3),'-rs')
% semilogy(per_case44(:,1),per_case44(:,3),'-r^'), hold off
% grid minor
% title({'WRONG NOISE POWER SETTING VS CORRECT NOISE POWER SETTING, RESOURCE: 3x1, LONG CB, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('CORRECT V3', 'CORRECT V120', 'CORRECT V200', 'CORRECT V500', 'WRONG V3', 'WRONG V120', 'WRONG V200', 'WRONG V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case21(:,1),per_case21(:,3),'-bo'), hold on
% semilogy(per_case22(:,1),per_case22(:,3),'-bx')
% semilogy(per_case23(:,1),per_case23(:,3),'-bs')
% semilogy(per_case24(:,1),per_case24(:,3),'-b^')
% semilogy(per_case45(:,1),per_case45(:,3),'-ro')
% semilogy(per_case46(:,1),per_case46(:,3),'-rx')
% semilogy(per_case47(:,1),per_case47(:,3),'-rs')
% semilogy(per_case48(:,1),per_case48(:,3),'-r^'), hold off
% grid minor
% title({'WRONG NOISE POWER SETTING VS CORRECT NOISE POWER SETTING, RESOURCE: 3x1, LONG CB, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('CORRECT V3', 'CORRECT V120', 'CORRECT V200', 'CORRECT V500', 'WRONG V3', 'WRONG V120', 'WRONG V200', 'WRONG V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case25(:,1),per_case25(:,3),'-bo'), hold on
% semilogy(per_case26(:,1),per_case26(:,3),'-bx')
% semilogy(per_case27(:,1),per_case27(:,3),'-bs')
% semilogy(per_case28(:,1),per_case28(:,3),'-b^')
% semilogy(per_case49(:,1),per_case49(:,3),'-ro')
% semilogy(per_case50(:,1),per_case50(:,3),'-rx')
% semilogy(per_case51(:,1),per_case51(:,3),'-rs')
% semilogy(per_case52(:,1),per_case52(:,3),'-r^'), hold off
% grid minor
% title({'WRONG NOISE POWER SETTING VS CORRECT NOISE POWER SETTING, RESOURCE: 3x1, LONG CB, SINGLE USER';'(OFDM BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('CORRECT V3', 'CORRECT V120', 'CORRECT V200', 'CORRECT V500', 'WRONG V3', 'WRONG V120', 'WRONG V200', 'WRONG V500')
% axis([-1 19 1e-3 1])
% 
% figure
% semilogy(per_case29(:,1),per_case29(:,3),'-bo'), hold on
% semilogy(per_case30(:,1),per_case30(:,3),'-bx')
% semilogy(per_case31(:,1),per_case31(:,3),'-bs')
% semilogy(per_case32(:,1),per_case32(:,3),'-b^')
% semilogy(per_case53(:,1),per_case53(:,3),'-ro')
% semilogy(per_case54(:,1),per_case54(:,3),'-rx')
% semilogy(per_case55(:,1),per_case55(:,3),'-rs')
% semilogy(per_case56(:,1),per_case56(:,3),'-r^'), hold off
% grid minor
% title({'WRONG NOISE POWER SETTING VS CORRECT NOISE POWER SETTING, RESOURCE: 3x1, LONG CB, SINGLE USER';'(OFDM BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
% xlabel('SNR (dB)'), ylabel('PER')
% legend('CORRECT V3', 'CORRECT V120', 'CORRECT V200', 'CORRECT V500', 'WRONG V3', 'WRONG V120', 'WRONG V200', 'WRONG V500')
% axis([-1 19 1e-3 1])
