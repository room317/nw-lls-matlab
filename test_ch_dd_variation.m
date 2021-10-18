% test channel variation in delay-doppler domain 
% dump 'ch_eff_usr_dd' first

% check variables
if ~exist('ch_eff_usr_dd', 'var')
    error('Dump variables first!')
end

% set parameters
nsym = 14;

% ch_autocorr = ch_eff_usr_dd;
ch_autocorr = b;
ch_conv = zeros(size(ch_autocorr));
ch_est = zeros(size(ch_autocorr));
nsubc = size(ch_autocorr, 1)/nsym;

for i = 1:nsym
    for j = 1:nsubc
        ch_conv(:, (j-1)*nsym+i) = reshape(circshift(reshape(ch_autocorr(:, (i-1)*nsubc+j), nsubc, nsym), [-(j-1) -(i-1)]), [], 1);
    end
end

ch_est_idx = 1;
for i = 1:nsym
    for j = 1:nsubc
        ch_est(:, (i-1)*nsubc+j) = reshape(circshift(reshape(ch_conv(:, ch_est_idx), nsubc, nsym), [(j-1) (i-1)]), [], 1);
    end
end

figure, mesh(1:nsubc*nsym, 1:nsubc*nsym, abs(ch_autocorr))
figure, mesh(1:nsubc*nsym, 1:nsubc*nsym, abs(ch_est))


% fig1 = figure;
% fig2 = figure;
% for i = 1:nsym:size(ch_conv, 2)
%     figure(fig1)
%     subplot(2, 1, 1), plot(real(ch_conv(:, i)), '-b'), axis([0 size(ch_conv, 2) -1 1]), grid minor
%     subplot(2, 1, 2), plot(imag(ch_conv(:, i)), '-r'), axis([0 size(ch_conv, 2) -1 1]), grid minor
%     
%     figure(fig2)
%     plot(abs(ch_conv(:, i)), '-k'), axis([0 size(ch_conv, 2) 0 1]), grid minor
% end

figure
plot(real(ch_conv(:, 1)), '-b.'), hold on, plot(real(ch_conv(:, end)), '-r.'), hold off, grid minor
legend('channel gain @ position 1', 'channel gain @ position 2')
xlabel('correlation channel index'), ylabel('Channel gain (real)')
