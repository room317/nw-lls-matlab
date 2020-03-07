% Simulator: test_ch_conv_script
% TDL-A
% 0.1 us rms delay spread
% 3 km/h mobility

% [SNR, NUM PKTS, PER TFEQ, PER DDEQ]
% 1) OTFS with cyclic postfix / synch delayed (perfect tf-domain channel estimation)
% 2) OTFS with cyclic prefix / perfect synch (perfect tf-domain channel estimation)

% [SNR, NUM PKTS, SER TFEQ, SER DDEQ, PER TFEQ, PER DDEQ]
% 3) OTFS with cyclic postfix / synch delayed (channel estimation with pilot)
% 4) OTFS with cyclic prefix / perfect synch (channel estimation with pilot)

% [SNR, NUM PKTS, PER TFEQ, PER DDEQ]

% with cyclic postfix / synch delayed to prevent inter-symbol interference
% whole resource as pilot
per_1 = [...
    10.0          401     1.000000     1.000000
    11.0          401     1.000000     1.000000
    12.0          401     1.000000     1.000000
    13.0          402     1.000000     0.997512
    14.0          435     0.960920     0.921839
    15.0          689     0.728592     0.582003
    16.0         1000     0.384000     0.181000
    17.0         1000     0.199000     0.065000
    18.0         1000     0.120000     0.031000
    19.0         1000     0.099000     0.041000];

% with cyclic prefix / perfect synch
% whole resource as pilot
per_2 = [...
    10.0          401     1.000000     1.000000
    11.0          401     1.000000     1.000000
    12.0          401     1.000000     1.000000
    13.0          401     1.000000     1.000000
    14.0          432     0.956019     0.928241
    15.0          744     0.697581     0.538978
    16.0         1000     0.353000     0.162000
    17.0         1000     0.204000     0.049000
    18.0         1000     0.125000     0.025000
    19.0         1000     0.078000     0.016000];

per_otfs_tdla_v3_tfeq_postfix = per_1(:, [1 3]);
per_otfs_tdla_v3_ddeq_postfix = per_1(:, [1 4]);
per_otfs_tdla_v3_tfeq_prefix = per_2(:, [1 3]);
per_otfs_tdla_v3_ddeq_prefix = per_2(:, [1 4]);

figure
semilogy(per_otfs_tdla_v3_tfeq_postfix(:,1),per_otfs_tdla_v3_tfeq_postfix(:,2),'-ko'), hold on
semilogy(per_otfs_tdla_v3_tfeq_prefix(:,1),per_otfs_tdla_v3_tfeq_prefix(:,2),'-kx'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Uncoded, Perfect CH EST, F-T EQ')
xlabel('SNR (dB)'), ylabel('PER'), legend('cyclic postfix', 'cyclic prefix')

figure
semilogy(per_otfs_tdla_v3_ddeq_postfix(:,1),per_otfs_tdla_v3_ddeq_postfix(:,2),'-ko'), hold on
semilogy(per_otfs_tdla_v3_ddeq_prefix(:,1),per_otfs_tdla_v3_ddeq_prefix(:,2),'-kx'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Uncoded, Perfect CH EST, D-D EQ')
xlabel('SNR (dB)'), ylabel('PER'), legend('cyclic postfix', 'cyclic prefix')

per_3 = [...
    10.0   1167  0.1556  0.2976  0.3436  1.0000
    12.0   1744  0.1018  0.2187  0.2299  1.0000
    14.0   2660  0.0664  0.1545  0.1508  1.0000
    16.0   4001  0.0432  0.1010  0.1002  1.0000
    18.0   6127  0.0267  0.0588  0.0654  1.0000
    20.0   7769  0.0197  0.0331  0.0516  1.0000
    22.0  10000  0.0134  0.0182  0.0392  0.9238
    24.0  10000  0.0085  0.0132  0.0264  0.5460
    26.0  10000  0.0050  0.0108  0.0195  0.3610
    28.0  10000  0.0032  0.0097  0.0138  0.3209
    30.0  10000  0.0015  0.0090  0.0107  0.2898
    32.0  10000  0.0011  0.0094  0.0092  0.2806
    34.0  10000  0.0005  0.0098  0.0062  0.2675
    36.0  10000  0.0002  0.0085  0.0033  0.2659
    38.0  10000  0.0000  0.0083  0.0013  0.2650];

per_4 = [...
    10.0   1183  0.1453  0.7339  0.3390  1.0000
    12.0   1858  0.0961  0.7190  0.2158  1.0000
    14.0   2664  0.0675  0.7122  0.1505  1.0000
    16.0   3992  0.0419  0.7078  0.1005  1.0000
    18.0   5952  0.0287  0.7035  0.0674  1.0000
    20.0   7972  0.0202  0.7024  0.0503  1.0000
    22.0  10000  0.0124  0.7013  0.0378  1.0000
    24.0  10000  0.0075  0.7003  0.0239  1.0000
    26.0  10000  0.0052  0.7005  0.0198  1.0000
    28.0  10000  0.0031  0.6996  0.0151  1.0000
    30.0  10000  0.0019  0.6993  0.0107  1.0000
    32.0  10000  0.0010  0.6992  0.0078  1.0000
    34.0  10000  0.0007  0.6986  0.0062  1.0000
    36.0  10000  0.0002  0.6987  0.0039  1.0000
    38.0  10000  0.0001  0.6992  0.0023  1.0000];

ser_otfs_tdla_v3_tfeq_postfix_pilot = per_3(:, [1 3]);
ser_otfs_tdla_v3_ddeq_postfix_pilot = per_3(:, [1 4]);
ser_otfs_tdla_v3_tfeq_prefix_pilot = per_4(:, [1 3]);
ser_otfs_tdla_v3_ddeq_prefix_pilot = per_4(:, [1 4]);
per_otfs_tdla_v3_tfeq_postfix_pilot = per_3(:, [1 5]);
per_otfs_tdla_v3_ddeq_postfix_pilot = per_3(:, [1 6]);
per_otfs_tdla_v3_tfeq_prefix_pilot = per_4(:, [1 5]);
per_otfs_tdla_v3_ddeq_prefix_pilot = per_4(:, [1 6]);

figure
semilogy(ser_otfs_tdla_v3_tfeq_postfix_pilot(:,1),ser_otfs_tdla_v3_tfeq_postfix_pilot(:,2),'-ko'), hold on
semilogy(ser_otfs_tdla_v3_tfeq_prefix_pilot(:,1),ser_otfs_tdla_v3_tfeq_prefix_pilot(:,2),'-kx'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Uncoded, CH EST with Pilot, F-T EQ')
xlabel('SNR (dB)'), ylabel('SER'), legend('cyclic postfix', 'cyclic prefix')

figure
semilogy(ser_otfs_tdla_v3_ddeq_postfix_pilot(:,1),ser_otfs_tdla_v3_ddeq_postfix_pilot(:,2),'-ko'), hold on
semilogy(ser_otfs_tdla_v3_ddeq_prefix_pilot(:,1),ser_otfs_tdla_v3_ddeq_prefix_pilot(:,2),'-kx'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Uncoded, CH EST with Pilot, D-D EQ')
xlabel('SNR (dB)'), ylabel('SER'), legend('cyclic postfix', 'cyclic prefix')

figure
semilogy(per_otfs_tdla_v3_tfeq_postfix_pilot(:,1),per_otfs_tdla_v3_tfeq_postfix_pilot(:,2),'-ko'), hold on
semilogy(per_otfs_tdla_v3_tfeq_prefix_pilot(:,1),per_otfs_tdla_v3_tfeq_prefix_pilot(:,2),'-kx'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Uncoded, CH EST with Pilot, F-T EQ')
xlabel('SNR (dB)'), ylabel('PER'), legend('cyclic postfix', 'cyclic prefix')

figure
semilogy(per_otfs_tdla_v3_ddeq_postfix_pilot(:,1),per_otfs_tdla_v3_ddeq_postfix_pilot(:,2),'-ko'), hold on
semilogy(per_otfs_tdla_v3_ddeq_prefix_pilot(:,1),per_otfs_tdla_v3_ddeq_prefix_pilot(:,2),'-kx'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Uncoded, CH EST with Pilot, D-D EQ')
xlabel('SNR (dB)'), ylabel('PER'), legend('cyclic postfix', 'cyclic prefix')

