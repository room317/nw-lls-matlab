% simulator: new_wave_lls
% commit: 

% RUNNING CASE: CASE: OTFS BW10 FC4000 V120 DS0.1 TDLA MCS8 LEN6120 SIM100 SNR0 DDTONE TFEQMMSE
% normal impulse
per_case01 = [ ...
  0.0       1.000000     1.297
  2.0       1.000000     1.030
  4.0       0.990000     0.796
  6.0       0.610000     0.665
  8.0       0.170000     0.548
 10.0       0.010000     0.443
 12.0       0.000000     0.382
 14.0       0.000000     0.316
 16.0       0.000000     0.272
 18.0       0.000000     0.228];

% reduced power impulse
per_case02 = [ ...
  0.0       1.000000    11.827
  2.0       1.000000     8.888
  4.0       1.000000     7.503
  6.0       1.000000     6.266
  8.0       1.000000     4.869
 10.0       1.000000     4.341
 12.0       1.000000     3.569
 14.0       1.000000     3.145
 16.0       1.000000     2.929
 18.0       1.000000     2.731];

% zadoff-chu spreading
per_case03 = [ ...
  0.0       1.000000     0.723
  2.0       1.000000     0.633
  4.0       0.990000     0.501
  6.0       0.740000     0.458
  8.0       0.360000     0.400
 10.0       0.090000     0.376
 12.0       0.060000     0.364
 14.0       0.070000     0.344
 16.0       0.210000     0.329
 18.0       0.440000     0.333];

figure
semilogy(per_case01(:,1),per_case01(:,2),'-ko'), hold on
semilogy(per_case03(:,1),per_case03(:,2),'-kx'), hold off
grid minor
title({'PER performance wrt pilot mapping'; 'OTFS BW10 FC4000 V120 DS0.1 TDLA MCS8 LEN6120 TFEQMMSE'})
xlabel('SNR (dB)'), ylabel('PER')
legend('Normalized pilot tone', 'Zadoff-Chu tone spreading')
axis([-1 19 1e-3 1])

figure
semilogy(per_case01(:,1),per_case01(:,3),'-ko'), hold on
% semilogy(per_case02(:,1),per_case02(:,3),'-kx')
semilogy(per_case03(:,1),per_case03(:,3),'-kx'), hold off
grid minor
title({'Channel estimation RMSE performance wrt pilot mapping'; 'OTFS BW10 FC4000 V120 DS0.1 TDLA MCS8 LEN6120 TFEQMMSE'})
xlabel('SNR (dB)'), ylabel('Channel estimation RMSE')
legend('Normalized pilot tone', 'Zadoff-Chu tone spreading')
axis([-1 19 1e-1 2])
