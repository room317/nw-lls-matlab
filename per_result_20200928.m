% simulator: new_wave_lls
% commit: 0143f17
% bandwidth 10 MHz, subcarrier spacing 60 khz, 4 slots, single user, 3 rbs x 1 slot per user
%   1. single user vs. multi users
%   2. single subframe vs. multi subframes
%   3. fully spreaded otfs vs. partially-spreaded otfs
%   4. otfs vs. ofdm

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case01 = [ ...
     0.0       1633   0.245560
     2.0       2545   0.157564
     4.0       4197   0.095544
     6.0       7186   0.055803
     8.0      10000   0.033900
    10.0      10000   0.015900
    12.0      10000   0.006600
    14.0      10000   0.003000
    16.0      10000   0.000900];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case02 = [ ...
     0.0       1502   0.266977
     2.0       2509   0.159825
     4.0       3934   0.101932
     6.0       7450   0.053826
     8.0      10000   0.026900
    10.0      10000   0.009500
    12.0      10000   0.003800
    14.0      10000   0.001300
    16.0      10000   0.000200];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case03 = [ ...
     0.0       1482   0.270580
     2.0       2525   0.158812
     4.0       3834   0.104591
     6.0       8688   0.046156
     8.0      10000   0.017800
    10.0      10000   0.005600
    12.0      10000   0.001300
    14.0      10000   0.000200];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case04 = [ ...
     0.0       1403   0.285816
     2.0       2961   0.135427
     4.0       7076   0.056670
     6.0      10000   0.021300
     8.0      10000   0.005800
    10.0      10000   0.002200
    12.0      10000   0.000400];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case05 = [ ...
     0.0       1880   0.213298
     2.0       3087   0.129900
     4.0       4936   0.081240
     6.0       7788   0.051489
     8.0      10000   0.037600
    10.0      10000   0.023000
    12.0      10000   0.017400
    14.0      10000   0.008700
    16.0      10000   0.005000
    18.0      10000   0.003000];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case06 = [ ...
     0.0       1748   0.229405
     2.0       3015   0.133002
     4.0       4041   0.099233
     6.0       6664   0.060174
     8.0      10000   0.036800
    10.0      10000   0.024500
    12.0      10000   0.010700
    14.0      10000   0.003800
    16.0      10000   0.001500
    18.0      10000   0.000700];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case07 = [ ...
     0.0       1713   0.234092
     2.0       2619   0.153112
     4.0       3818   0.105029
     6.0       5974   0.067124
     8.0      10000   0.036200
    10.0      10000   0.015200
    12.0      10000   0.007100
    14.0      10000   0.003800
    16.0      10000   0.001100
    18.0      10000   0.000600];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case08 = [ ...
     0.0       1414   0.283593
     2.0       2611   0.153581
     4.0       5903   0.067932
     6.0      10000   0.034700
     8.0      10000   0.013900
    10.0      10000   0.005400
    12.0      10000   0.001600
    14.0      10000   0.000700];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case09 = [ ...
     0.0       1033   0.398838  0.391094  0.389158  0.388190
     2.0       1848   0.219697  0.216991  0.246753  0.240801
     4.0       3311   0.142857  0.146481  0.123830  0.121111
     6.0       4696   0.085605  0.085392  0.088160  0.088586
     8.0       9101   0.045490  0.045160  0.044061  0.044720
    10.0      10000   0.022600  0.022500  0.020200  0.021400
    12.0      10000   0.011300  0.011500  0.010800  0.010900
    14.0      10000   0.004700  0.005100  0.004000  0.003700
    16.0      10000   0.001800  0.001200  0.002000  0.001400
    18.0      10000   0.000600  0.000800  0.000600  0.000400];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case10 = [ ...
     0.0        990   0.405051  0.448485  0.426263  0.441414
     2.0       1722   0.264808  0.232869  0.263066  0.255517
     4.0       2534   0.158248  0.161800  0.175612  0.172060
     6.0       4069   0.098550  0.101991  0.103957  0.100025
     8.0       8636   0.048518  0.053150  0.046434  0.050486
    10.0      10000   0.021000  0.020800  0.019000  0.021200
    12.0      10000   0.009000  0.006600  0.007100  0.006500
    14.0      10000   0.001200  0.002400  0.001900  0.002000
    16.0      10000   0.000200  0.000500  0.000100  0.000400];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case11 = [ ...
     0.0        893   0.462486  0.483763  0.449048  0.484882
     2.0       1425   0.303860  0.329123  0.281404  0.299649
     4.0       2170   0.184793  0.194009  0.193088  0.196774
     6.0       3861   0.103859  0.113960  0.109816  0.117327
     8.0       8140   0.054914  0.050246  0.049263  0.049386
    10.0      10000   0.019300  0.018700  0.018700  0.018700
    12.0      10000   0.005900  0.005100  0.005700  0.005700
    14.0      10000   0.000600  0.001400  0.000900  0.001000
    16.0      10000   0.000200  0.000200  0.000000  0.000500];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case12 = [ ...
     0.0        739   0.608931  0.573748  0.589986  0.542625
     2.0       1048   0.382634  0.399809  0.418893  0.395038
     4.0       1946   0.223535  0.242549  0.206064  0.226619
     6.0       4333   0.096238  0.094853  0.092546  0.097161
     8.0      10000   0.032200  0.034800  0.036900  0.034600
    10.0      10000   0.010400  0.011900  0.011200  0.010600
    12.0      10000   0.002800  0.004000  0.002900  0.004000
    14.0      10000   0.000600  0.000300  0.000300  0.000600];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case13 = [ ...
     0.0       1418   0.283498  0.282793  0.307475  0.309591
     2.0       2588   0.165379  0.166924  0.157651  0.154946
     4.0       4433   0.090458  0.090684  0.093390  0.093616
     6.0       7151   0.056775  0.056076  0.065585  0.065445
     8.0       9586   0.041832  0.041832  0.043814  0.044127
    10.0      10000   0.025400  0.025600  0.024700  0.023800
    12.0      10000   0.015300  0.016000  0.015000  0.014500
    14.0      10000   0.009400  0.009000  0.011600  0.012000
    16.0      10000   0.005700  0.006400  0.005100  0.005400
    18.0      10000   0.003000  0.003300  0.002800  0.002600];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case14 = [ ...
     0.0       1164   0.363402  0.344502  0.350515  0.344502
     2.0       2012   0.199304  0.201292  0.209245  0.218688
     4.0       3267   0.124579  0.137435  0.122743  0.131007
     6.0       5024   0.079817  0.084594  0.083997  0.087381
     8.0       7674   0.054079  0.053297  0.057727  0.052254
    10.0      10000   0.031900  0.033000  0.030600  0.033700
    12.0      10000   0.015200  0.016400  0.016400  0.019100
    14.0      10000   0.007700  0.006500  0.008100  0.006500
    16.0      10000   0.001900  0.001200  0.001400  0.002200
    18.0      10000   0.000500  0.000800  0.000700  0.000800];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case15 = [ ...
     0.0       1013   0.412636  0.401777  0.419546  0.395854
     2.0       1661   0.249849  0.251656  0.241421  0.261288
     4.0       2755   0.148094  0.152813  0.145554  0.152813
     6.0       4211   0.095227  0.101164  0.100689  0.107575
     8.0       6650   0.062707  0.060301  0.062256  0.062857
    10.0      10000   0.031200  0.031500  0.036800  0.032200
    12.0      10000   0.014700  0.012600  0.014000  0.013700
    14.0      10000   0.005500  0.004100  0.005100  0.005200
    16.0      10000   0.001600  0.001400  0.001800  0.001600
    18.0      10000   0.000200  0.000200  0.000600  0.000200];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case16 = [ ...
     0.0        738   0.550136  0.546070  0.543360  0.555556
     2.0       1107   0.362240  0.362240  0.364047  0.368564
     4.0       2032   0.197343  0.207677  0.208169  0.205217
     6.0       4120   0.105583  0.097330  0.105340  0.103641
     8.0       8719   0.045992  0.046106  0.048629  0.047712
    10.0      10000   0.017600  0.019000  0.016600  0.018300
    12.0      10000   0.006700  0.007200  0.005800  0.007100
    14.0      10000   0.001900  0.002800  0.002400  0.003200
    16.0      10000   0.001200  0.000900  0.000800  0.000700
    18.0      10000   0.000100  0.000200  0.000200  0.000100];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case17 = [ ...
     0.0       1753   0.228751
     2.0       3086   0.129942
     4.0       4717   0.085012
     6.0       8663   0.046289
     8.0      10000   0.019800
    10.0      10000   0.006800
    12.0      10000   0.001800
    14.0      10000   0.000300];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case18 = [ ...
     0.0       3856   0.103994
     2.0      10000   0.017100
     4.0      10000   0.002300
     6.0      10000   0.000100];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case19 = [ ...
     0.0       4252   0.094309
     2.0      10000   0.011100
     4.0      10000   0.000800];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case20 = [ ...
     0.0       7527   0.053275
     2.0      10000   0.002300
     4.0      10000   0.000100];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case21 = [ ...
     0.0       1966   0.203967
     2.0       3330   0.120420
     4.0       4597   0.087231
     6.0       7861   0.051011
     8.0      10000   0.034600
    10.0      10000   0.022200
    12.0      10000   0.010800
    14.0      10000   0.004400
    16.0      10000   0.001000
    18.0      10000   0.001000];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case22 = [ ...
     0.0       2864   0.140014
     2.0      10000   0.035100
     4.0      10000   0.008000
     6.0      10000   0.000900];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case23 = [ ...
     0.0       3141   0.127666
     2.0      10000   0.025400
     4.0      10000   0.003000
     6.0      10000   0.000200];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case24 = [ ...
     0.0       5989   0.066956
     2.0      10000   0.007000
     4.0      10000   0.000600];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case25 = [ ...
     0.0       1825   0.219726
     2.0       3127   0.128238
     4.0       5118   0.078351
     6.0       7896   0.050785
     8.0      10000   0.019800
    10.0      10000   0.009300
    12.0      10000   0.002800];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case26 = [ ...
     0.0       3739   0.107248
     2.0      10000   0.018100
     4.0      10000   0.002700
     6.0      10000   0.000100];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case27 = [ ...
     0.0       4045   0.099135
     2.0      10000   0.011800
     4.0      10000   0.000900];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case28 = [ ...
     0.0       6973   0.057508
     2.0      10000   0.002400
     4.0      10000   0.000000];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case29 = [ ...
     0.0       2138   0.187558
     2.0       3288   0.121959
     4.0       4635   0.086516
     6.0       7548   0.053127
     8.0      10000   0.037000
    10.0      10000   0.018700
    12.0      10000   0.010600
    14.0      10000   0.005500
    16.0      10000   0.002300
    18.0      10000   0.001000];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case30 = [ ...
     0.0       2973   0.134881
     2.0      10000   0.034300
     4.0      10000   0.006900
     6.0      10000   0.001100
     8.0      10000   0.000200];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case31 = [ ...
     0.0       2930   0.136860
     2.0      10000   0.026700
     4.0      10000   0.003700
     6.0      10000   0.000600];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case32 = [ ...
     0.0       5672   0.070698
     2.0      10000   0.006600
     4.0      10000   0.000600];

% RUNNING CASE: OTFS SU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case33 = [ ...
     0.0        401   1.000000
     2.0        707   0.567185
     4.0       1340   0.299254
     6.0       2937   0.136534
     8.0       7214   0.055586
    10.0      10000   0.016300
    12.0      10000   0.004200
    14.0      10000   0.000900];

% RUNNING CASE: OTFS SU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case34 = [ ...
     0.0        401   1.000000
     2.0        568   0.705986
     4.0       1447   0.277125
     6.0       5286   0.075861
     8.0      10000   0.011500
    10.0      10000   0.001500
    12.0      10000   0.000000];

% RUNNING CASE: OTFS SU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case35 = [ ...
     0.0        401   1.000000
     2.0        493   0.813387
     4.0       1614   0.248451
     6.0       8322   0.048186
     8.0      10000   0.005900
    10.0      10000   0.000300];

% RUNNING CASE: OTFS SU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case36 = [ ...
     0.0        401   1.000000
     2.0        466   0.860515
     4.0       2095   0.191408
     6.0      10000   0.013200
     8.0      10000   0.000100];

% RUNNING CASE: OTFS SU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case37 = [ ...
     0.0        411   0.975669
     2.0        868   0.461982
     4.0       1657   0.242004
     6.0       2546   0.157502
     8.0       4392   0.091302
    10.0       8021   0.049994
    12.0      10000   0.024300
    14.0      10000   0.010300
    16.0      10000   0.004800
    18.0      10000   0.001300];

% RUNNING CASE: OTFS SU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case38 = [ ...
     0.0        401   1.000000
     2.0        653   0.614089
     4.0       1387   0.289113
     6.0       3932   0.101984
     8.0       9145   0.043849
    10.0      10000   0.011600
    12.0      10000   0.002800
    14.0      10000   0.000400];

% RUNNING CASE: OTFS SU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case39 = [ ...
     0.0        401   1.000000
     2.0        579   0.692573
     4.0       1565   0.256230
     6.0       4739   0.084617
     8.0      10000   0.023000
    10.0      10000   0.004800
    12.0      10000   0.000900];

% RUNNING CASE: OTFS SU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case40 = [ ...
     0.0        401   1.000000
     2.0        502   0.798805
     4.0       1721   0.233004
     6.0      10000   0.039200
     8.0      10000   0.004700
    10.0      10000   0.000300];

% basic graphs

% figure
% semilogy(per_case01(:,1),per_case01(:,3),'-ko'), hold on
% semilogy(per_case02(:,1),per_case02(:,3),'-kx')
% semilogy(per_case03(:,1),per_case03(:,3),'-ks')
% semilogy(per_case04(:,1),per_case04(:,3),'-k^'), hold off
% grid minor
% title({'OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN964 REAL TFEQMMSE';'SINGLE USER (3x1)'})
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

figure
semilogy(per_case01(:,1),per_case01(:,3),'-ko'), hold on
semilogy(per_case09(:,1),per_case09(:,3),':ko')
semilogy(per_case09(:,1),per_case09(:,4),':bo')
semilogy(per_case09(:,1),per_case09(:,5),':ro')
semilogy(per_case09(:,1),per_case09(:,6),':go'), hold off
grid minor
title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, SHORT CB, 3 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN964 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
axis([-1 19 1e-3 1])

figure
semilogy(per_case04(:,1),per_case04(:,3),'-ko'), hold on
semilogy(per_case12(:,1),per_case12(:,3),':kx')
semilogy(per_case12(:,1),per_case12(:,4),':bx')
semilogy(per_case12(:,1),per_case12(:,5),':rx')
semilogy(per_case12(:,1),per_case12(:,6),':gx'), hold off
grid minor
title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, SHORT CB, 500 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN964 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
axis([-1 19 1e-3 1])

figure
semilogy(per_case05(:,1),per_case05(:,3),'-ko'), hold on
semilogy(per_case13(:,1),per_case13(:,3),':ko')
semilogy(per_case13(:,1),per_case13(:,4),':bo')
semilogy(per_case13(:,1),per_case13(:,5),':ro')
semilogy(per_case13(:,1),per_case13(:,6),':go'), hold off
grid minor
title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, SHORT CB, 3 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN964 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
axis([-1 19 1e-3 1])

figure
semilogy(per_case08(:,1),per_case08(:,3),'-ko'), hold on
semilogy(per_case16(:,1),per_case16(:,3),':kx')
semilogy(per_case16(:,1),per_case16(:,4),':bx')
semilogy(per_case16(:,1),per_case16(:,5),':rx')
semilogy(per_case16(:,1),per_case16(:,6),':gx'), hold off
grid minor
title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, SHORT CB, 500 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN964 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
axis([-1 19 1e-3 1])

figure
semilogy(per_case01(:,1),per_case01(:,3),'-ko'), hold on
semilogy(per_case02(:,1),per_case02(:,3),'-kx')
semilogy(per_case03(:,1),per_case03(:,3),'-ks')
semilogy(per_case04(:,1),per_case04(:,3),'-k^')
semilogy(per_case17(:,1),per_case17(:,3),'-ro')
semilogy(per_case18(:,1),per_case18(:,3),'-rx')
semilogy(per_case19(:,1),per_case19(:,3),'-rs')
semilogy(per_case20(:,1),per_case20(:,3),'-r^'), hold off
grid minor
title({'SHORT CODE BLOCK VS LONG CODE BLOCK, RESOURCE: 3x1, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('V3 (SHORT CB)', 'V120 (SHORT CB)', 'V200 (SHORT CB)', 'V500 (SHORT CB)', 'V3 (LONG CB)', 'V120 (LONG CB)', 'V200 (LONG CB)', 'V500 (LONG CB)')
axis([-1 19 1e-3 1])

figure
semilogy(per_case05(:,1),per_case05(:,3),'-ko'), hold on
semilogy(per_case06(:,1),per_case06(:,3),'-kx')
semilogy(per_case07(:,1),per_case07(:,3),'-ks')
semilogy(per_case08(:,1),per_case08(:,3),'-k^')
semilogy(per_case21(:,1),per_case21(:,3),'-ro')
semilogy(per_case22(:,1),per_case22(:,3),'-rx')
semilogy(per_case23(:,1),per_case23(:,3),'-rs')
semilogy(per_case24(:,1),per_case24(:,3),'-r^'), hold off
grid minor
title({'SHORT CODE BLOCK VS LONG CODE BLOCK, RESOURCE: 3x1, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('V3 (SHORT CB)', 'V120 (SHORT CB)', 'V200 (SHORT CB)', 'V500 (SHORT CB)', 'V3 (LONG CB)', 'V120 (LONG CB)', 'V200 (LONG CB)', 'V500 (LONG CB)')
axis([-1 19 1e-3 1])

figure
semilogy(per_case17(:,1),per_case17(:,3),'-bo'), hold on
semilogy(per_case18(:,1),per_case18(:,3),'-bx')
semilogy(per_case19(:,1),per_case19(:,3),'-bs')
semilogy(per_case20(:,1),per_case20(:,3),'-b^')
semilogy(per_case25(:,1),per_case25(:,3),'-ro')
semilogy(per_case26(:,1),per_case26(:,3),'-rx')
semilogy(per_case27(:,1),per_case27(:,3),'-rs')
semilogy(per_case28(:,1),per_case28(:,3),'-r^'), hold off
grid minor
title({'OFDM VS OTFS, RESOURCE: 3x1, LONG CB, SINGLE USER';'(BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('OTFS V3', 'OTFS V120', 'OTFS V200', 'OTFS V500', 'OFDM V3', 'OFDM V120', 'OFDM V200', 'OFDM V500')
axis([-1 19 1e-3 1])

figure
semilogy(per_case21(:,1),per_case21(:,3),'-bo'), hold on
semilogy(per_case22(:,1),per_case22(:,3),'-bx')
semilogy(per_case23(:,1),per_case23(:,3),'-bs')
semilogy(per_case24(:,1),per_case24(:,3),'-b^')
semilogy(per_case29(:,1),per_case29(:,3),'-ro')
semilogy(per_case30(:,1),per_case30(:,3),'-rx')
semilogy(per_case31(:,1),per_case31(:,3),'-rs')
semilogy(per_case32(:,1),per_case32(:,3),'-r^'), hold off
grid minor
title({'OFDM VS OTFS, RESOURCE: 3x1, LONG CB, SINGLE USER';'(BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('OTFS V3', 'OTFS V120', 'OTFS V200', 'OTFS V500', 'OFDM V3', 'OFDM V120', 'OFDM V200', 'OFDM V500')
axis([-1 19 1e-3 1])

figure
semilogy(per_case17(:,1),per_case17(:,3),'-bo'), hold on
semilogy(per_case18(:,1),per_case18(:,3),'-bx')
semilogy(per_case19(:,1),per_case19(:,3),'-bs')
semilogy(per_case20(:,1),per_case20(:,3),'-b^')
semilogy(per_case33(:,1),per_case33(:,3),'-ro')
semilogy(per_case34(:,1),per_case34(:,3),'-rx')
semilogy(per_case35(:,1),per_case35(:,3),'-rs')
semilogy(per_case36(:,1),per_case36(:,3),'-r^'), hold off
grid minor
title({'FULL SPREADING(50x4) VS PARTIAL SPREADING(3x1), LONG CB, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('PARTIAL V3', 'PARTIAL V120', 'PARTIAL V200', 'PARTIAL V500', 'FULL V3', 'FULL V120', 'FULL V200', 'FULL V500')
axis([-1 19 1e-3 1])

figure
semilogy(per_case21(:,1),per_case21(:,3),'-bo'), hold on
semilogy(per_case22(:,1),per_case22(:,3),'-bx')
semilogy(per_case23(:,1),per_case23(:,3),'-bs')
semilogy(per_case24(:,1),per_case24(:,3),'-b^')
semilogy(per_case37(:,1),per_case37(:,3),'-ro')
semilogy(per_case38(:,1),per_case38(:,3),'-rx')
semilogy(per_case39(:,1),per_case39(:,3),'-rs')
semilogy(per_case40(:,1),per_case40(:,3),'-r^'), hold off
grid minor
title({'FULL SPREADING(50x4) VS PARTIAL SPREADING(3x1), LONG CB, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('PARTIAL V3', 'PARTIAL V120', 'PARTIAL V200', 'PARTIAL V500', 'FULL V3', 'FULL V120', 'FULL V200', 'FULL V500')
axis([-1 19 1e-3 1])
