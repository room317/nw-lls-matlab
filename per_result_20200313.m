% Simulator: test_ch_conv_script
% TDL-A
% 0.1 us rms delay spread
% 3 km/h mobility

% [SNR, NUM PKTS, SER TFEQ, SER DDEQ, PER TFEQ, PER DDEQ]
% 1) OTFS with cyclic postfix / pilot size 16 (channel estimation with pilot)
% 2) OTFS with cyclic postfix / pilot size 20 (channel estimation with pilot)
% 3) OTFS with cyclic postfix / pilot size 24 (channel estimation with pilot)
% 4) OTFS with cyclic postfix / pilot size 28 (channel estimation with pilot)

% with cyclic postfix / synch delayed to prevent inter-symbol interference
% pilot size: 16x14 / guard: 8x14
per_pilot16 = [...
    10.0    100  0.1303  0.2895  0.3800  1.0000
    12.0    100  0.1187  0.2234  0.2100  1.0000
    14.0    100  0.0762  0.1482  0.1500  1.0000
    16.0    100  0.0608  0.0973  0.1300  1.0000
    18.0    100  0.0269  0.0641  0.0500  1.0000
    20.0    100  0.0365  0.0445  0.0900  1.0000
    22.0    100  0.0202  0.0230  0.0700  0.9200
    24.0    100  0.0181  0.0260  0.0400  0.5800
    26.0    100  0.0053  0.0143  0.0200  0.3500
    28.0    100  0.0000  0.0030  0.0000  0.3500];

% with cyclic postfix / synch delayed to prevent inter-symbol interference
% pilot size: 20x14 / guard: 10x14
per_pilot20 = [...
    10.0   1000  0.1629  0.2190  0.3470  1.0000
    12.0   1000  0.1141  0.1355  0.2350  1.0000
    14.0   1000  0.0676  0.0740  0.1470  1.0000
    16.0   1000  0.0535  0.0393  0.1060  1.0000
    18.0   1000  0.0215  0.0146  0.0610  0.9450
    20.0   1000  0.0132  0.0071  0.0450  0.4790
    22.0   1000  0.0090  0.0061  0.0280  0.2230
    24.0   1000  0.0076  0.0063  0.0280  0.1500
    26.0   1000  0.0047  0.0057  0.0150  0.1470
    28.0   1000  0.0034  0.0080  0.0170  0.1430];
  
% with cyclic postfix / synch delayed to prevent inter-symbol interference
% pilot size: 24x14 / guard: 12x14
per_pilot24 = [...
    10.0   1000  0.1580  0.1531  0.3220  1.0000
    12.0   1000  0.1008  0.0735  0.2240  1.0000
    14.0   1000  0.0710  0.0303  0.1550  0.9980
    16.0   1000  0.0413  0.0126  0.0890  0.7350
    18.0   1000  0.0253  0.0062  0.0520  0.2570
    20.0   1000  0.0200  0.0080  0.0420  0.1180
    22.0   1000  0.0147  0.0069  0.0350  0.1080
    24.0   1000  0.0065  0.0023  0.0200  0.0940
    26.0   1000  0.0027  0.0025  0.0090  0.0980
    28.0   1000  0.0016  0.0028  0.0070  0.0900];

% with cyclic postfix / synch delayed to prevent inter-symbol interference
% pilot size: 28x14 / guard: 14x14
per_pilot28 = [...
    10.0   1000  0.1655  0.0943  0.3230  1.0000
    12.0   1000  0.1100  0.0357  0.2200  0.9610
    14.0   1000  0.0794  0.0148  0.1460  0.4640
    16.0   1000  0.0462  0.0085  0.1010  0.1480
    18.0   1000  0.0273  0.0049  0.0580  0.0830
    20.0   1000  0.0182  0.0032  0.0410  0.0740
    22.0   1000  0.0139  0.0048  0.0330  0.0670
    24.0   1000  0.0092  0.0037  0.0200  0.0670
    26.0   1000  0.0085  0.0041  0.0330  0.0670
    28.0   1000  0.0023  0.0026  0.0160  0.0610];

ser_otfs_tdla_v3_pilot16 = per_pilot16(:, [1 4]);
per_otfs_tdla_v3_pilot16 = per_pilot16(:, [1 6]);
ser_otfs_tdla_v3_pilot20 = per_pilot20(:, [1 4]);
per_otfs_tdla_v3_pilot20 = per_pilot20(:, [1 6]);
ser_otfs_tdla_v3_pilot24 = per_pilot24(:, [1 4]);
per_otfs_tdla_v3_pilot24 = per_pilot24(:, [1 6]);
ser_otfs_tdla_v3_pilot28 = per_pilot28(:, [1 4]);
per_otfs_tdla_v3_pilot28 = per_pilot28(:, [1 6]);

figure
semilogy(ser_otfs_tdla_v3_pilot16(:,1),ser_otfs_tdla_v3_pilot16(:,2),'-ko'), hold on
semilogy(ser_otfs_tdla_v3_pilot20(:,1),ser_otfs_tdla_v3_pilot20(:,2),'-kx')
semilogy(ser_otfs_tdla_v3_pilot24(:,1),ser_otfs_tdla_v3_pilot24(:,2),'-ks')
semilogy(ser_otfs_tdla_v3_pilot28(:,1),ser_otfs_tdla_v3_pilot28(:,2),'-k^'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Uncoded, CH EST with Pilot, D-D EQ')
xlabel('SNR (dB)'), ylabel('SER'), legend('pilot size: 16x14', 'pilot size: 20x14', 'pilot size: 24x14', 'pilot size: 28x14')

figure
semilogy(per_otfs_tdla_v3_pilot16(:,1),per_otfs_tdla_v3_pilot16(:,2),'-ko'), hold on
semilogy(per_otfs_tdla_v3_pilot20(:,1),per_otfs_tdla_v3_pilot20(:,2),'-kx')
semilogy(per_otfs_tdla_v3_pilot24(:,1),per_otfs_tdla_v3_pilot24(:,2),'-ks')
semilogy(per_otfs_tdla_v3_pilot28(:,1),per_otfs_tdla_v3_pilot28(:,2),'-k^'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Uncoded, CH EST with Pilot, D-D EQ')
xlabel('SNR (dB)'), ylabel('PER'), legend('pilot size: 16x14', 'pilot size: 20x14', 'pilot size: 24x14', 'pilot size: 28x14')
