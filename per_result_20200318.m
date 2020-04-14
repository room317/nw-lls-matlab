% Simulator: test_ch_conv_script
% TDL-A
% 0.01 us rms delay spread
% 3 km/h mobility

% [SNR, NUM PKTS, SER TFEQ, SER DDEQ, PER TFEQ, PER DDEQ]
% 1) OTFS with cyclic postfix / pilot size 16 (channel estimation with pilot)
% 2) OTFS with cyclic postfix / pilot size 20 (channel estimation with pilot)
% 3) OTFS with cyclic postfix / pilot size 24 (channel estimation with pilot)
% 4) OTFS with cyclic postfix / pilot size 28 (channel estimation with pilot)

% with cyclic postfix / synch delayed to prevent inter-symbol interference
% pilot size: 16x14 / guard: 8x14
per_pilot16 = [...
    10.0   1000  0.1554  0.2894  0.3440  1.0000
    12.0   1000  0.0901  0.2101  0.1950  1.0000
    14.0   1000  0.0692  0.1442  0.1430  1.0000
    16.0   1000  0.0368  0.0887  0.0770  1.0000
    18.0   1000  0.0254  0.0446  0.0560  1.0000
    20.0   1000  0.0169  0.0167  0.0420  1.0000
    22.0   1000  0.0133  0.0063  0.0320  0.7790
    24.0   1000  0.0060  0.0010  0.0160  0.1670
    26.0   1000  0.0069  0.0013  0.0190  0.0250
    28.0   1000  0.0049  0.0013  0.0090  0.0100];

% with cyclic postfix / synch delayed to prevent inter-symbol interference
% pilot size: 20x14 / guard: 10x14
per_pilot20 = [...
    10.0   1000  0.1503  0.2106  0.3310  1.0000
    12.0   1000  0.1000  0.1247  0.2140  1.0000
    14.0   1000  0.0610  0.0637  0.1290  1.0000
    16.0   1000  0.0443  0.0266  0.0970  1.0000
    18.0   1000  0.0251  0.0078  0.0530  0.9040
    20.0   1000  0.0194  0.0017  0.0400  0.2850
    22.0   1000  0.0098  0.0020  0.0250  0.0380
    24.0   1000  0.0088  0.0017  0.0270  0.0140
    26.0   1000  0.0084  0.0008  0.0250  0.0090
    28.0   1000  0.0028  0.0004  0.0120  0.0060];
  
% with cyclic postfix / synch delayed to prevent inter-symbol interference
% pilot size: 24x14 / guard: 12x14
per_pilot24 = [...
    10.0   1000  0.1503  0.1424  0.3270  1.0000
    12.0   1000  0.1045  0.0658  0.2110  1.0000
    14.0   1000  0.0747  0.0246  0.1490  0.9910
    16.0   1000  0.0412  0.0064  0.0900  0.6850
    18.0   1000  0.0283  0.0018  0.0680  0.1270
    20.0   1000  0.0157  0.0006  0.0370  0.0160
    22.0   1000  0.0107  0.0001  0.0290  0.0040
    24.0   1000  0.0125  0.0001  0.0220  0.0050
    26.0   1000  0.0044  0.0007  0.0130  0.0030
    28.0   1000  0.0070  0.0003  0.0160  0.0050];

% with cyclic postfix / synch delayed to prevent inter-symbol interference
% pilot size: 28x14 / guard: 14x14
per_pilot28 = [...
    10.0   1000  0.1381  0.0806  0.2880  1.0000
    12.0   1000  0.1010  0.0305  0.1940  0.9510
    14.0   1000  0.0630  0.0082  0.1280  0.4130
    16.0   1000  0.0425  0.0036  0.0790  0.0480
    18.0   1000  0.0296  0.0008  0.0600  0.0080
    20.0   1000  0.0203  0.0018  0.0390  0.0060
    22.0   1000  0.0131  0.0013  0.0310  0.0030
    24.0   1000  0.0051  0.0000  0.0190  0.0000
    26.0   1000  0.0057  0.0004  0.0140  0.0010
    28.0   1000  0.0042  0.0000  0.0080  0.0000];

% with cyclic postfix / synch delayed to prevent inter-symbol interference
% pilot size: 28x14 / guard: 18x14
per_pilot28a = [...
    10.0   1000  0.1470  0.0807  0.2960  1.0000
    12.0   1000  0.1059  0.0265  0.2120  0.9510
    14.0   1000  0.0688  0.0063  0.1280  0.3980
    16.0   1000  0.0451  0.0030  0.0890  0.0530
    18.0   1000  0.0278  0.0016  0.0580  0.0110
    20.0   1000  0.0199  0.0009  0.0470  0.0070
    22.0   1000  0.0136  0.0006  0.0280  0.0040
    24.0   1000  0.0059  0.0004  0.0150  0.0020
    26.0   1000  0.0054  0.0005  0.0140  0.0040
    28.0   1000  0.0060  0.0002  0.0160  0.0030];

ser_otfs_tdla_v3_pilot16 = per_pilot16(:, [1 4]);
per_otfs_tdla_v3_pilot16 = per_pilot16(:, [1 6]);
ser_otfs_tdla_v3_pilot20 = per_pilot20(:, [1 4]);
per_otfs_tdla_v3_pilot20 = per_pilot20(:, [1 6]);
ser_otfs_tdla_v3_pilot24 = per_pilot24(:, [1 4]);
per_otfs_tdla_v3_pilot24 = per_pilot24(:, [1 6]);
ser_otfs_tdla_v3_pilot28 = per_pilot28(:, [1 4]);
per_otfs_tdla_v3_pilot28 = per_pilot28(:, [1 6]);
ser_otfs_tdla_v3_pilot28a = per_pilot28a(:, [1 4]);
per_otfs_tdla_v3_pilot28a = per_pilot28a(:, [1 6]);

figure
semilogy(ser_otfs_tdla_v3_pilot16(:,1),ser_otfs_tdla_v3_pilot16(:,2),'-ko'), hold on
semilogy(ser_otfs_tdla_v3_pilot20(:,1),ser_otfs_tdla_v3_pilot20(:,2),'-kx')
semilogy(ser_otfs_tdla_v3_pilot24(:,1),ser_otfs_tdla_v3_pilot24(:,2),'-ks')
semilogy(ser_otfs_tdla_v3_pilot28(:,1),ser_otfs_tdla_v3_pilot28(:,2),'-k^')
semilogy(ser_otfs_tdla_v3_pilot28a(:,1),ser_otfs_tdla_v3_pilot28a(:,2),'-kd'), hold off
grid minor
title('2D-Spread Mod/Demod, TDL-A, 16QAM, Uncoded')
xlabel('SNR (dB)'), ylabel('SER')
legend('pilot resource: 25%', 'pilot resource: 31%', 'pilot resource: 38%', 'pilot resource: 44%', 'pilot resource: 56%')
axis([10 24 5e-3 1])

figure
semilogy(per_otfs_tdla_v3_pilot16(:,1),per_otfs_tdla_v3_pilot16(:,2),'-ko'), hold on
semilogy(per_otfs_tdla_v3_pilot20(:,1),per_otfs_tdla_v3_pilot20(:,2),'-kx')
semilogy(per_otfs_tdla_v3_pilot24(:,1),per_otfs_tdla_v3_pilot24(:,2),'-ks')
semilogy(per_otfs_tdla_v3_pilot28(:,1),per_otfs_tdla_v3_pilot28(:,2),'-k^')
semilogy(per_otfs_tdla_v3_pilot28a(:,1),per_otfs_tdla_v3_pilot28a(:,2),'-kd'), hold off
grid minor
title('OTFS, TDL-A, 16QAM, Uncoded, CH EST with Pilot, D-D EQ')
xlabel('SNR (dB)'), ylabel('PER')
legend('pilot resource: 25%', 'pilot resource: 31%', 'pilot resource: 38%', 'pilot resource: 44%', 'pilot resource: 56%')
axis([10 28 1e-2 1])
