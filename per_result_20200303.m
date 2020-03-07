% Simulator: new_wave_test_r3

% 1) OFDM with perfect tf-domain channel estimation
% 2) OTFS with perfect tf-domain channel estimation
% 3) OFDM with tf-domain channel estimation using 2 pilots
% 4) OTFS with dd-domain channel estimation using impulse train pilots
% 5) OFDM with real channel equalization
% 6) OTFS with real channel equalization

per_1 = [...
      0.0        401   1.000000
      2.0        407   0.985258
      4.0        702   0.571225
      6.0       1151   0.348393
      8.0       2445   0.164008
     10.0       6122   0.065501
     12.0      10000   0.019300
     14.0      10000   0.005100
     16.0      10000   0.001200
     18.0      10000   0.000100
      0.0        401   1.000000
      2.0        414   0.968599
      4.0        589   0.680815
      6.0       1058   0.379017
      8.0       2277   0.176109
     10.0       9005   0.044531
     12.0      10000   0.005900
     14.0      10000   0.000700
      0.0        401   1.000000
      2.0        404   0.992574
      4.0        554   0.723827
      6.0       1033   0.388190
      8.0       2785   0.143986
     10.0      10000   0.027900
     12.0      10000   0.004200
     14.0      10000   0.000100
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        494   0.811741
      6.0       1182   0.339255
      8.0       4894   0.081937
     10.0      10000   0.007400
     12.0      10000   0.000100
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        525   0.763810
      6.0       1602   0.250312
      8.0       6543   0.061287
     10.0      10000   0.009600
     12.0      10000   0.001100
     14.0      10000   0.000000
      0.0        401   1.000000
      2.0        402   0.997512
      4.0        514   0.780156
      6.0       1484   0.270216
      8.0       8292   0.048360
     10.0      10000   0.004200
     12.0      10000   0.000000
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        491   0.816701
      6.0       1590   0.252201
      8.0      10000   0.027900
     10.0      10000   0.001100
     12.0      10000   0.000000
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        445   0.901124
      6.0       1867   0.214783
      8.0      10000   0.008700
     10.0      10000   0.000000];

per_ofdm_tdla_v3_tf_perfect = per_1(1:10, [1 3]);
per_ofdm_tdla_v120_tf_perfect = per_1(11:18, [1 3]);
per_ofdm_tdla_v200_tf_perfect = per_1(19:26, [1 3]);
per_ofdm_tdla_v500_tf_perfect = per_1(27:33, [1 3]);
per_ofdm_tdld_v3_tf_perfect = per_1(34:41, [1 3]);
per_ofdm_tdld_v120_tf_perfect = per_1(42:48, [1 3]);
per_ofdm_tdld_v200_tf_perfect = per_1(49:55, [1 3]);
per_ofdm_tdld_v500_tf_perfect = per_1(56:61, [1 3]);

figure
semilogy(per_ofdm_tdla_v3_tf_perfect(:,1),per_ofdm_tdla_v3_tf_perfect(:,2),'-ko'), hold on
semilogy(per_ofdm_tdla_v120_tf_perfect(:,1),per_ofdm_tdla_v120_tf_perfect(:,2),'-kx')
semilogy(per_ofdm_tdla_v200_tf_perfect(:,1),per_ofdm_tdla_v200_tf_perfect(:,2),'-ks')
semilogy(per_ofdm_tdla_v500_tf_perfect(:,1),per_ofdm_tdla_v500_tf_perfect(:,2),'-k^'), hold off
grid minor
title('OFDM, TDL-A, 16QAM, Rate: 0.48, Perfect Ch. Est.')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

figure
semilogy(per_ofdm_tdld_v3_tf_perfect(:,1),per_ofdm_tdld_v3_tf_perfect(:,2),'-ko'), hold on
semilogy(per_ofdm_tdld_v120_tf_perfect(:,1),per_ofdm_tdld_v120_tf_perfect(:,2),'-kx')
semilogy(per_ofdm_tdld_v200_tf_perfect(:,1),per_ofdm_tdld_v200_tf_perfect(:,2),'-ks')
semilogy(per_ofdm_tdld_v500_tf_perfect(:,1),per_ofdm_tdld_v500_tf_perfect(:,2),'-k^'), hold off
grid minor
title('OFDM, TDL-D, 16QAM, Rate: 0.48, Perfect Ch. Est.')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

per_2 = [...
      0.0        401   1.000000
      2.0        402   0.997512
      4.0        639   0.627543
      6.0       1186   0.338111
      8.0       2457   0.163207
     10.0       5709   0.070240
     12.0      10000   0.021100
     14.0      10000   0.005200
     16.0      10000   0.001800
     18.0      10000   0.000300
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        555   0.722523
      6.0       1069   0.375117
      8.0       3341   0.120024
     10.0      10000   0.022900
     12.0      10000   0.002100
     14.0      10000   0.000300
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        474   0.845992
      6.0       1130   0.354867
      8.0       4718   0.084994
     10.0      10000   0.013800
     12.0      10000   0.000600
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        419   0.957041
      6.0       1195   0.335565
      8.0      10000   0.032100
     10.0      10000   0.001100
     12.0      10000   0.000000
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        521   0.769674
      6.0       1330   0.301504
      8.0       5283   0.075904
     10.0      10000   0.011500
     12.0      10000   0.000900
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        438   0.915525
      6.0       1489   0.269308
      8.0      10000   0.022000
     10.0      10000   0.000900
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        427   0.939110
      6.0       1697   0.236299
      8.0      10000   0.007500
     10.0      10000   0.000100
      0.0        401   1.000000
      2.0        401   1.000000
      4.0        405   0.990123
      6.0       2198   0.182439
      8.0      10000   0.000600];

per_otfs_tdla_v3_tf_perfect = per_2(1:10, [1 3]);
per_otfs_tdla_v120_tf_perfect = per_2(11:18, [1 3]);
per_otfs_tdla_v200_tf_perfect = per_2(19:25, [1 3]);
per_otfs_tdla_v500_tf_perfect = per_2(26:32, [1 3]);
per_otfs_tdld_v3_tf_perfect = per_2(33:39, [1 3]);
per_otfs_tdld_v120_tf_perfect = per_2(40:45, [1 3]);
per_otfs_tdld_v200_tf_perfect = per_2(46:51, [1 3]);
per_otfs_tdld_v500_tf_perfect = per_2(52:56, [1 3]);

figure
semilogy(per_otfs_tdla_v3_tf_perfect(:,1),per_otfs_tdla_v3_tf_perfect(:,2),'-ko'), hold on
semilogy(per_otfs_tdla_v120_tf_perfect(:,1),per_otfs_tdla_v120_tf_perfect(:,2),'-kx')
semilogy(per_otfs_tdla_v200_tf_perfect(:,1),per_otfs_tdla_v200_tf_perfect(:,2),'-ks')
semilogy(per_otfs_tdla_v500_tf_perfect(:,1),per_otfs_tdla_v500_tf_perfect(:,2),'-k^'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Rate: 0.48, Perfect Ch. Est.')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

figure
semilogy(per_otfs_tdld_v3_tf_perfect(:,1),per_otfs_tdld_v3_tf_perfect(:,2),'-ko'), hold on
semilogy(per_otfs_tdld_v120_tf_perfect(:,1),per_otfs_tdld_v120_tf_perfect(:,2),'-kx')
semilogy(per_otfs_tdld_v200_tf_perfect(:,1),per_otfs_tdld_v200_tf_perfect(:,2),'-ks')
semilogy(per_otfs_tdld_v500_tf_perfect(:,1),per_otfs_tdld_v500_tf_perfect(:,2),'-k^'), hold off
grid minor
title('OTFS, TDL-D, 16QAM, Rate: 0.48, Perfect Ch. Est.')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

per_3 = [...
];

per_ofdm_tdla_v3_tf_pilot = per_3(1:10, [1 3]);
per_ofdm_tdla_v120_tf_pilot = per_3(11:20, [1 3]);
per_ofdm_tdla_v200_tf_pilot = per_3(11:30, [1 3]);
per_ofdm_tdla_v500_tf_pilot = per_3(31:40, [1 3]);
per_ofdm_tdld_v3_tf_pilot = per_3(41:50, [1 3]);
per_ofdm_tdld_v120_tf_pilot = per_3(51:60, [1 3]);
per_ofdm_tdld_v200_tf_pilot = per_3(61:70, [1 3]);
per_ofdm_tdld_v500_tf_pilot = per_3(71:80, [1 3]);

figure
semilogy(per_ofdm_tdla_v3_tf_pilot(:,1),per_ofdm_tdla_v3_tf_pilot(:,2),'-ko'), hold on
semilogy(per_ofdm_tdla_v120_tf_pilot(:,1),per_ofdm_tdla_v120_tf_pilot(:,2),'-kx')
semilogy(per_ofdm_tdla_v200_tf_pilot(:,1),per_ofdm_tdla_v200_tf_pilot(:,2),'-ks')
semilogy(per_ofdm_tdla_v500_tf_pilot(:,1),per_ofdm_tdla_v500_tf_pilot(:,2),'-k^'), hold off
grid minor
title('OFDM, TDL-A, 16QAM, Rate: 0.48, Ch. Est. with Pilots')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

figure
semilogy(per_ofdm_tdld_v3_tf_pilot(:,1),per_ofdm_tdld_v3_tf_pilot(:,2),'-ko'), hold on
semilogy(per_ofdm_tdld_v120_tf_pilot(:,1),per_ofdm_tdld_v120_tf_pilot(:,2),'-kx')
semilogy(per_ofdm_tdld_v200_tf_pilot(:,1),per_ofdm_tdld_v200_tf_pilot(:,2),'-ks')
semilogy(per_ofdm_tdld_v500_tf_pilot(:,1),per_ofdm_tdld_v500_tf_pilot(:,2),'-k^'), hold off
grid minor
title('OFDM, TDL-D, 16QAM, Rate: 0.48, Ch. Est. with Pilots')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

per_4 = [...
];

per_otfs_tdla_v3_dd_pilot = per_4(1:10, [1 3]);
per_otfs_tdla_v120_dd_pilot = per_4(11:20, [1 3]);
per_otfs_tdla_v200_dd_pilot = per_4(21:30, [1 3]);
per_otfs_tdla_v500_dd_pilot = per_4(31:40, [1 3]);
per_otfs_tdld_v3_dd_pilot = per_4(41:50, [1 3]);
per_otfs_tdld_v120_dd_pilot = per_4(51:60, [1 3]);
per_otfs_tdld_v200_dd_pilot = per_4(61:70, [1 3]);
per_otfs_tdld_v500_dd_pilot = per_4(71:80, [1 3]);

figure
semilogy(per_otfs_tdla_v3_dd_pilot(:,1),per_otfs_tdla_v3_dd_pilot(:,2),'-ko'), hold on
semilogy(per_otfs_tdla_v120_dd_pilot(:,1),per_otfs_tdla_v120_dd_pilot(:,2),'-kx')
semilogy(per_otfs_tdla_v200_dd_pilot(:,1),per_otfs_tdla_v200_dd_pilot(:,2),'-ks')
semilogy(per_otfs_tdla_v500_dd_pilot(:,1),per_otfs_tdla_v500_dd_pilot(:,2),'-k^'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Rate: 0.48, Ch. Est. with Pilots')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

figure
semilogy(per_otfs_tdld_v3_dd_pilot(:,1),per_otfs_tdld_v3_dd_pilot(:,2),'-ko'), hold on
semilogy(per_otfs_tdld_v120_dd_pilot(:,1),per_otfs_tdld_v120_dd_pilot(:,2),'-kx')
semilogy(per_otfs_tdld_v200_dd_pilot(:,1),per_otfs_tdld_v200_dd_pilot(:,2),'-ks')
semilogy(per_otfs_tdld_v500_dd_pilot(:,1),per_otfs_tdld_v500_dd_pilot(:,2),'-k^'), hold off
grid minor
title('OTFS, TDL-D, 16QAM, Rate: 0.48, Ch. Est. with Pilots')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

per_5 = [...
];

per_ofdm_tdla_v3_real_ch = per_5(1:10, [1 3]);
per_ofdm_tdla_v120_real_ch = per_5(11:20, [1 3]);
per_ofdm_tdla_v200_real_ch = per_5(21:30, [1 3]);
per_ofdm_tdla_v500_real_ch = per_5(31:40, [1 3]);
per_ofdm_tdld_v3_real_ch = per_5(41:50, [1 3]);
per_ofdm_tdld_v120_real_ch = per_5(51:60, [1 3]);
per_ofdm_tdld_v200_real_ch = per_5(61:70, [1 3]);
per_ofdm_tdld_v500_real_ch = per_5(71:80, [1 3]);

figure
semilogy(per_ofdm_tdla_v3_real_ch(:,1),per_ofdm_tdla_v3_real_ch(:,2),'-ko'), hold on
semilogy(per_ofdm_tdla_v120_real_ch(:,1),per_ofdm_tdla_v120_real_ch(:,2),'-kx')
semilogy(per_ofdm_tdla_v200_real_ch(:,1),per_ofdm_tdla_v200_real_ch(:,2),'-ks')
semilogy(per_ofdm_tdla_v500_real_ch(:,1),per_ofdm_tdla_v500_real_ch(:,2),'-k^'), hold off
grid minor
title('OFDM, TDL-A, 16QAM, Rate: 0.48, Real Ch.')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

figure
semilogy(per_ofdm_tdld_v3_real_ch(:,1),per_ofdm_tdld_v3_real_ch(:,2),'-ko'), hold on
semilogy(per_ofdm_tdld_v120_real_ch(:,1),per_ofdm_tdld_v120_real_ch(:,2),'-kx')
semilogy(per_ofdm_tdld_v200_real_ch(:,1),per_ofdm_tdld_v200_real_ch(:,2),'-ks')
semilogy(per_ofdm_tdld_v500_real_ch(:,1),per_ofdm_tdld_v500_real_ch(:,2),'-k^'), hold off
grid minor
title('OFDM, TDL-D, 16QAM, Rate: 0.48, Real Ch.')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

per_6 = [...
];

per_otfs_tdla_v3_real_ch = per_6(1:10, [1 3]);
per_otfs_tdla_v120_real_ch = per_6(11:20, [1 3]);
per_otfs_tdla_v200_real_ch = per_6(21:30, [1 3]);
per_otfs_tdla_v500_real_ch = per_6(31:39, [1 3]);
per_otfs_tdld_v3_real_ch = per_6(40:49, [1 3]);
per_otfs_tdld_v120_real_ch = per_6(50:59, [1 3]);
per_otfs_tdld_v200_real_ch = per_6(60:69, [1 3]);
per_otfs_tdld_v500_real_ch = per_6(70:77, [1 3]);

figure
semilogy(per_otfs_tdla_v3_real_ch(:,1),per_otfs_tdla_v3_real_ch(:,2),'-ko'), hold on
semilogy(per_otfs_tdla_v120_real_ch(:,1),per_otfs_tdla_v120_real_ch(:,2),'-kx')
semilogy(per_otfs_tdla_v200_real_ch(:,1),per_otfs_tdla_v200_real_ch(:,2),'-ks')
semilogy(per_otfs_tdla_v500_real_ch(:,1),per_otfs_tdla_v500_real_ch(:,2),'-k^'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Rate: 0.48, Real Ch.')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')

figure
semilogy(per_otfs_tdld_v3_real_ch(:,1),per_otfs_tdld_v3_real_ch(:,2),'-ko'), hold on
semilogy(per_otfs_tdld_v120_real_ch(:,1),per_otfs_tdld_v120_real_ch(:,2),'-kx')
semilogy(per_otfs_tdld_v200_real_ch(:,1),per_otfs_tdld_v200_real_ch(:,2),'-ks')
semilogy(per_otfs_tdld_v500_real_ch(:,1),per_otfs_tdld_v500_real_ch(:,2),'-k^'), hold off
grid minor
title('OTFS, TDL-D, 16QAM, Rate: 0.48, Real Ch.')
xlabel('SNR (dB)'), ylabel('PER'), legend('3 km/h', '120 km/h', '200 km/h', '500 km/h')
