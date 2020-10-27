% simulator: new_wave_lls
% commit: 60cfc4a
% bandwidth 10 MHz, subcarrier spacing 60 khz, 4 slots, single user, 3 rbs x 1 slot per user
% comment: previous result(per_result_20200928) has wrong noise power setting.
%          this is noise power matched version.
%   1. single user vs. multi users
%   2. single subframe vs. multi subframes (short tb vs. long tb)
%   3. fully spreaded otfs vs. partially-spreaded otfs
%   4. otfs vs. ofdm

%% otfs / single user / short transmit block / partial spreading

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case01 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        518   0.774131
     6.0       1299   0.308699
     8.0       2286   0.175416
    10.0       3808   0.105305
    12.0       6388   0.062774
    14.0      10000   0.029900
    16.0      10000   0.015600
    18.0      10000   0.005900];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case02 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        561   0.714795
     6.0       1202   0.333611
     8.0       1987   0.201812
    10.0       3419   0.117286
    12.0       6467   0.062007
    14.0      10000   0.027700
    16.0      10000   0.008000
    18.0      10000   0.002100];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case03 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        511   0.784736
     6.0       1062   0.377589
     8.0       1997   0.200801
    10.0       3428   0.116978
    12.0       7587   0.052854
    14.0      10000   0.016500
    16.0      10000   0.006700
    18.0      10000   0.001500];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case04 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        425   0.943529
     6.0        766   0.523499
     8.0       1537   0.260898
    10.0       3506   0.114375
    12.0      10000   0.030300
    14.0      10000   0.008000
    16.0      10000   0.000600];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case05 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        665   0.603008
     6.0       1593   0.251726
     8.0       2681   0.149571
    10.0       4877   0.082223
    12.0       8052   0.049801
    14.0      10000   0.031800
    16.0      10000   0.021800
    18.0      10000   0.015600];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case06 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        615   0.652033
     6.0       1570   0.255414
     8.0       2595   0.154528
    10.0       4117   0.097401
    12.0       6522   0.061484
    14.0       9254   0.043333
    16.0      10000   0.022700
    18.0      10000   0.009800];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case07 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        544   0.737132
     6.0       1345   0.298141
     8.0       2284   0.175569
    10.0       3412   0.117526
    12.0       5976   0.067102
    14.0      10000   0.035800
    16.0      10000   0.019300
    18.0      10000   0.006600];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN964 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case08 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        454   0.883260
     6.0        815   0.492025
     8.0       1519   0.263989
    10.0       3403   0.117837
    12.0       8367   0.047926
    14.0      10000   0.017300
    16.0      10000   0.005800
    18.0      10000   0.001600];

%% otfs / multi user / long transmit block / partial spreading

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR4

per_case09 = [ ...
     0.0        401   1.000000  1.000000  1.000000  1.000000
     2.0        411   0.978102  0.975669  0.985401  0.985401
     4.0        778   0.515424  0.524422  0.537275  0.541131
     6.0       1483   0.271746  0.270398  0.277815  0.277815
     8.0       2342   0.171648  0.171221  0.182323  0.183177
    10.0       4042   0.100693  0.099208  0.101682  0.101682
    12.0       7032   0.057025  0.057309  0.059585  0.059158
    14.0      10000   0.028100  0.028900  0.025600  0.026500
    16.0      10000   0.012000  0.011600  0.011900  0.011700
    18.0      10000   0.004000  0.004200  0.003300  0.002600];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case10 = [ ...
     0.0        401   1.000000  1.000000  1.000000  1.000000
     2.0        401   1.000000  1.000000  1.000000  1.000000
     4.0        477   0.863732  0.840671  0.863732  0.855346
     6.0       1248   0.322115  0.321314  0.322115  0.341346
     8.0       7322   0.054766  0.058864  0.056132  0.060366
    10.0      10000   0.004500  0.005300  0.005300  0.007000
    12.0      10000   0.000300  0.000200  0.000200  0.000200];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case11 = [ ...
     0.0        401   1.000000  1.000000  1.000000  1.000000
     2.0        401   1.000000  1.000000  1.000000  1.000000
     4.0        471   0.851380  0.870488  0.885350  0.855626
     6.0       1094   0.366545  0.393053  0.378428  0.372943
     8.0       4813   0.083732  0.090796  0.083316  0.087056
    10.0      10000   0.009300  0.007100  0.005400  0.005700
    12.0      10000   0.000200  0.000100  0.000200  0.000600];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case12 = [ ...
     0.0        401   1.000000  1.000000  1.000000  1.000000
     2.0        401   1.000000  1.000000  1.000000  1.000000
     4.0        429   0.955711  0.939394  0.976690  0.934732
     6.0        855   0.490058  0.469006  0.481871  0.480702
     8.0       5130   0.083236  0.098441  0.078168  0.081871
    10.0      10000   0.003100  0.003500  0.004100  0.003500
    12.0      10000   0.000000  0.000000  0.000000  0.000000];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case13 = [ ...
     0.0        401   1.000000  1.000000  1.000000  1.000000
     2.0        402   1.000000  1.000000  1.000000  0.997512
     4.0        972   0.412551  0.422840  0.434156  0.433128
     6.0       1924   0.225052  0.225052  0.208420  0.209979
     8.0       3067   0.133355  0.133355  0.131073  0.130747
    10.0       5029   0.082323  0.084908  0.080334  0.079738
    12.0       7074   0.057676  0.058666  0.056686  0.057959
    14.0      10000   0.038900  0.038900  0.036600  0.035500
    16.0      10000   0.023000  0.022500  0.022400  0.022000
    18.0      10000   0.010300  0.010300  0.009400  0.009300];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case14 = [ ...
     0.0        401   1.000000  1.000000  1.000000  1.000000
     2.0        401   1.000000  1.000000  1.000000  1.000000
     4.0        515   0.778641  0.798058  0.801942  0.780583
     6.0       1431   0.292802  0.324249  0.280224  0.306778
     8.0       5831   0.074258  0.075802  0.072372  0.068770
    10.0      10000   0.012800  0.013700  0.012000  0.011500
    12.0      10000   0.002000  0.002000  0.001600  0.001300
    14.0      10000   0.000000  0.000000  0.000100  0.000000];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case15 = [ ...
     0.0        401   1.000000  1.000000  1.000000  1.000000
     2.0        401   1.000000  1.000000  1.000000  1.000000
     4.0        490   0.828571  0.848980  0.818367  0.846939
     6.0       1213   0.330585  0.383347  0.347898  0.353669
     8.0       4970   0.081690  0.084507  0.080684  0.083903
    10.0      10000   0.013600  0.012700  0.013000  0.013200
    12.0      10000   0.001000  0.001200  0.001200  0.001000
    14.0      10000   0.000100  0.000100  0.000000  0.000200];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR4
per_case16 = [ ...
     0.0        401   1.000000  1.000000  1.000000  1.000000
     2.0        401   1.000000  1.000000  1.000000  1.000000
     4.0        436   0.931193  0.926606  0.928899  0.919725
     6.0        910   0.448352  0.440659  0.457143  0.446154
     8.0       4528   0.092314  0.090327  0.093198  0.088560
    10.0      10000   0.007500  0.007700  0.006000  0.007100
    12.0      10000   0.000300  0.000500  0.000400  0.000600];

%% otfs / single user / long transmit block / partial spreading

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case17 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        853   0.470106
     6.0       1554   0.258044
     8.0       2487   0.161238
    10.0       4507   0.088973
    12.0       8841   0.045357
    14.0      10000   0.019400
    16.0      10000   0.006900
    18.0      10000   0.002800];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case18 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        460   0.871739
     6.0       1632   0.245711
     8.0      10000   0.037100
    10.0      10000   0.003500
    12.0      10000   0.000000];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case19 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        443   0.905192
     6.0       1560   0.257051
     8.0      10000   0.033300
    10.0      10000   0.001300
    12.0      10000   0.000000];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case20 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        412   0.973301
     6.0       1134   0.353616
     8.0      10000   0.018200
    10.0      10000   0.000200];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case21 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        958   0.418580
     6.0       2079   0.192881
     8.0       3188   0.125784
    10.0       4792   0.083681
    12.0       7310   0.054856
    14.0      10000   0.035200
    16.0      10000   0.019800
    18.0      10000   0.008700];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case22 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        545   0.735780
     6.0       1554   0.258044
     8.0       7203   0.055671
    10.0      10000   0.009500
    12.0      10000   0.002100
    14.0      10000   0.000100];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case23 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        488   0.821721
     6.0       1470   0.272789
     8.0       8030   0.049938
    10.0      10000   0.006200
    12.0      10000   0.000900];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case24 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        413   0.970944
     6.0       1254   0.319777
     8.0      10000   0.032000
    10.0      10000   0.001600
    12.0      10000   0.000000];

%% ofdm / single user / long transmit block / partial spreading

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case25 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        852   0.470657
     6.0       1560   0.257051
     8.0       2457   0.163207
    10.0       4429   0.090540
    12.0       8448   0.047467
    14.0      10000   0.021600
    16.0      10000   0.008000
    18.0      10000   0.003800];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case26 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        468   0.856838
     6.0       1612   0.248759
     8.0      10000   0.037900
    10.0      10000   0.003300
    12.0      10000   0.000100];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case27 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        459   0.873638
     6.0       1514   0.264861
     8.0      10000   0.031200
    10.0      10000   0.001900
    12.0      10000   0.000100];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case28 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        410   0.978049
     6.0       1292   0.310372
     8.0      10000   0.019900
    10.0      10000   0.000200];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case29 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        934   0.429336
     6.0       1969   0.203657
     8.0       3199   0.125352
    10.0       4121   0.097306
    12.0       7576   0.052930
    14.0      10000   0.039300
    16.0      10000   0.019600
    18.0      10000   0.011100];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case30 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        525   0.763810
     6.0       1619   0.247684
     8.0       6307   0.063580
    10.0      10000   0.010000
    12.0      10000   0.001100
    14.0      10000   0.000000];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case31 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        499   0.803607
     6.0       1442   0.278086
     8.0       7918   0.050644
    10.0      10000   0.007000
    12.0      10000   0.000500];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case32 = [ ...
     0.0        401   1.000000
     2.0        401   1.000000
     4.0        424   0.945755
     6.0       1353   0.296378
     8.0      10000   0.029900
    10.0      10000   0.001000
    12.0      10000   0.000000];

%% otfs / single user / long transmit block / full spreading

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

%% otfs / single user / long transmit block / partial spreading / wrong noise power

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case41 = [ ...
     0.0       1753   0.228751
     2.0       3086   0.129942
     4.0       4717   0.085012
     6.0       8663   0.046289
     8.0      10000   0.019800
    10.0      10000   0.006800
    12.0      10000   0.001800
    14.0      10000   0.000300];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE USR1
per_case42 = [ ...
     0.0       3856   0.103994
     2.0      10000   0.017100
     4.0      10000   0.002300
     6.0      10000   0.000100];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE USR1
per_case43 = [ ...
     0.0       4252   0.094309
     2.0      10000   0.011100
     4.0      10000   0.000800];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE USR1
per_case44 = [ ...
     0.0       7527   0.053275
     2.0      10000   0.002300
     4.0      10000   0.000100];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE USR1
per_case45 = [ ...
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

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE USR1
per_case46 = [ ...
     0.0       2864   0.140014
     2.0      10000   0.035100
     4.0      10000   0.008000
     6.0      10000   0.000900];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case47 = [ ...
     0.0       3141   0.127666
     2.0      10000   0.025400
     4.0      10000   0.003000
     6.0      10000   0.000200];

% RUNNING CASE: OTFS MU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE USR1
per_case48 = [ ...
     0.0       5989   0.066956
     2.0      10000   0.007000
     4.0      10000   0.000600];

%% ofdm / single user / long transmit block / partial spreading / wrong noise power

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V3 DS0.1 TDLC MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE
per_case49 = [ ...
     0.0       1825   0.219726
     2.0       3127   0.128238
     4.0       5118   0.078351
     6.0       7896   0.050785
     8.0      10000   0.019800
    10.0      10000   0.009300
    12.0      10000   0.002800];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V120 DS0.1 TDLC MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE
per_case50 = [ ...
     0.0       3739   0.107248
     2.0      10000   0.018100
     4.0      10000   0.002700
     6.0      10000   0.000100];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V200 DS0.1 TDLC MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE
per_case51 = [ ...
     0.0       4045   0.099135
     2.0      10000   0.011800
     4.0      10000   0.000900];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V500 DS0.1 TDLC MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE
per_case52 = [ ...
     0.0       6973   0.057508
     2.0      10000   0.002400
     4.0      10000   0.000000];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V3 DS0.01 TDLE MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE
per_case53 = [ ...
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

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V120 DS0.01 TDLE MCS8 LEN6120 SIM10 SNR0 REAL TFEQMMSE
per_case54 = [ ...
     0.0       2973   0.134881
     2.0      10000   0.034300
     4.0      10000   0.006900
     6.0      10000   0.001100
     8.0      10000   0.000200];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V200 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case55 = [ ...
     0.0       2930   0.136860
     2.0      10000   0.026700
     4.0      10000   0.003700
     6.0      10000   0.000600];

% RUNNING CASE: OFDM MU BW10 SCS60 SLOT4 FC4000 V500 DS0.01 TDLE MCS8 LEN6120 SIM10000 SNR0 REAL TFEQMMSE
per_case56 = [ ...
     0.0       5672   0.070698
     2.0      10000   0.006600
     4.0      10000   0.000600];

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
semilogy(per_case17(:,1),per_case17(:,3),'-ko'), hold on
semilogy(per_case09(:,1),per_case09(:,3),':ko')
semilogy(per_case09(:,1),per_case09(:,4),':bo')
semilogy(per_case09(:,1),per_case09(:,5),':ro')
semilogy(per_case09(:,1),per_case09(:,6),':go'), hold off
grid minor
title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, LONG CB, 3 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
axis([-1 19 1e-3 1])

figure
semilogy(per_case20(:,1),per_case20(:,3),'-ko'), hold on
semilogy(per_case12(:,1),per_case12(:,3),':kx')
semilogy(per_case12(:,1),per_case12(:,4),':bx')
semilogy(per_case12(:,1),per_case12(:,5),':rx')
semilogy(per_case12(:,1),per_case12(:,6),':gx'), hold off
grid minor
title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, LONG CB, 500 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
axis([-1 19 1e-3 1])

figure
semilogy(per_case21(:,1),per_case21(:,3),'-ko'), hold on
semilogy(per_case13(:,1),per_case13(:,3),':ko')
semilogy(per_case13(:,1),per_case13(:,4),':bo')
semilogy(per_case13(:,1),per_case13(:,5),':ro')
semilogy(per_case13(:,1),per_case13(:,6),':go'), hold off
grid minor
title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, LONG CB, 3 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('SINGLE-USER', 'USER 1 OF 4', 'USER 2 OF 4', 'USER 3 OF 4', 'USER 4 OF 4')
axis([-1 19 1e-3 1])

figure
semilogy(per_case24(:,1),per_case24(:,3),'-ko'), hold on
semilogy(per_case16(:,1),per_case16(:,3),':kx')
semilogy(per_case16(:,1),per_case16(:,4),':bx')
semilogy(per_case16(:,1),per_case16(:,5),':rx')
semilogy(per_case16(:,1),per_case16(:,6),':gx'), hold off
grid minor
title({'SINGLE-USER VS MULTI-USER, RESOURCE: 3x1, LONG CB, 500 km/h';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
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

figure
semilogy(per_case17(:,1),per_case17(:,3),'-bo'), hold on
semilogy(per_case18(:,1),per_case18(:,3),'-bx')
semilogy(per_case19(:,1),per_case19(:,3),'-bs')
semilogy(per_case20(:,1),per_case20(:,3),'-b^')
semilogy(per_case41(:,1),per_case41(:,3),'-ro')
semilogy(per_case42(:,1),per_case42(:,3),'-rx')
semilogy(per_case43(:,1),per_case43(:,3),'-rs')
semilogy(per_case44(:,1),per_case44(:,3),'-r^'), hold off
grid minor
title({'WRONG NOISE POWER SETTING VS CORRECT NOISE POWER SETTING, RESOURCE: 3x1, LONG CB, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('CORRECT V3', 'CORRECT V120', 'CORRECT V200', 'CORRECT V500', 'WRONG V3', 'WRONG V120', 'WRONG V200', 'WRONG V500')
axis([-1 19 1e-3 1])

figure
semilogy(per_case21(:,1),per_case21(:,3),'-bo'), hold on
semilogy(per_case22(:,1),per_case22(:,3),'-bx')
semilogy(per_case23(:,1),per_case23(:,3),'-bs')
semilogy(per_case24(:,1),per_case24(:,3),'-b^')
semilogy(per_case45(:,1),per_case45(:,3),'-ro')
semilogy(per_case46(:,1),per_case46(:,3),'-rx')
semilogy(per_case47(:,1),per_case47(:,3),'-rs')
semilogy(per_case48(:,1),per_case48(:,3),'-r^'), hold off
grid minor
title({'WRONG NOISE POWER SETTING VS CORRECT NOISE POWER SETTING, RESOURCE: 3x1, LONG CB, SINGLE USER';'(OTFS BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('CORRECT V3', 'CORRECT V120', 'CORRECT V200', 'CORRECT V500', 'WRONG V3', 'WRONG V120', 'WRONG V200', 'WRONG V500')
axis([-1 19 1e-3 1])

figure
semilogy(per_case25(:,1),per_case25(:,3),'-bo'), hold on
semilogy(per_case26(:,1),per_case26(:,3),'-bx')
semilogy(per_case27(:,1),per_case27(:,3),'-bs')
semilogy(per_case28(:,1),per_case28(:,3),'-b^')
semilogy(per_case49(:,1),per_case49(:,3),'-ro')
semilogy(per_case50(:,1),per_case50(:,3),'-rx')
semilogy(per_case51(:,1),per_case51(:,3),'-rs')
semilogy(per_case52(:,1),per_case52(:,3),'-r^'), hold off
grid minor
title({'WRONG NOISE POWER SETTING VS CORRECT NOISE POWER SETTING, RESOURCE: 3x1, LONG CB, SINGLE USER';'(OFDM BW10 SCS60 SLOT4 FC4000 DS0.1 TDLC MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('CORRECT V3', 'CORRECT V120', 'CORRECT V200', 'CORRECT V500', 'WRONG V3', 'WRONG V120', 'WRONG V200', 'WRONG V500')
axis([-1 19 1e-3 1])

figure
semilogy(per_case29(:,1),per_case29(:,3),'-bo'), hold on
semilogy(per_case30(:,1),per_case30(:,3),'-bx')
semilogy(per_case31(:,1),per_case31(:,3),'-bs')
semilogy(per_case32(:,1),per_case32(:,3),'-b^')
semilogy(per_case53(:,1),per_case53(:,3),'-ro')
semilogy(per_case54(:,1),per_case54(:,3),'-rx')
semilogy(per_case55(:,1),per_case55(:,3),'-rs')
semilogy(per_case56(:,1),per_case56(:,3),'-r^'), hold off
grid minor
title({'WRONG NOISE POWER SETTING VS CORRECT NOISE POWER SETTING, RESOURCE: 3x1, LONG CB, SINGLE USER';'(OFDM BW10 SCS60 SLOT4 FC4000 DS0.01 TDLE MCS8 LEN6120 REAL TFEQMMSE)'})
xlabel('SNR (dB)'), ylabel('PER')
legend('CORRECT V3', 'CORRECT V120', 'CORRECT V200', 'CORRECT V500', 'WRONG V3', 'WRONG V120', 'WRONG V200', 'WRONG V500')
axis([-1 19 1e-3 1])
