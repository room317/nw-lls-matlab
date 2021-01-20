% generate sequence
[Ga, Gb] = wlanGolaySequence(32);

% system
x1 = Ga/sqrt(32);
x2 = Gb/sqrt(32);
h = complex(randn(16, 1), randn(16, 1));
% h = 1;
y1 = conv(h, x1);
y2 = conv(h, x2);

% demodulation
p1 = Ga/sqrt(32);
p2 = Gb/sqrt(32);
z = (xcorr(y1, p1)+xcorr(y2, p2));
hest = z(32+16-1:end)/2;

figure,
subplot(2, 2, 1), plot(xcorr(Ga, Ga), '-k'), grid minor, title('Autocorrelation of Ga(32)'), axis([0 64 -5 35]), text(40, 20, '(\itG_a \otimes \itG_a)')
subplot(2, 2, 2), plot(xcorr(Gb, Gb), '-k'), grid minor, title('Autocorrelation of Gb(32)'), axis([0 64 -5 35]), text(40, 20, '(\itG_b \otimes \itG_b)')
% subplot(2, 2, [3, 4]), plot(xcorr(Ga, Ga)+xcorr(Gb, Gb), '-b.'), grid minor, title('Sum of autocorrelations of two Complementary Golay Sequences(32)'), axis([0 64 -10 70])
pos = [0.35 0.1 0.35 0.35];
subplot('Position', pos), plot(xcorr(Ga, Ga)+xcorr(Gb, Gb), '-k'), grid minor, title('Sum of autocorrelations of two Complementary Golay Sequences(32)'), axis([0 64 -10 70]), text(40, 40, '(\itG_a \otimes \itG_a)+(\itG_b \otimes \itG_b)')

figure
subplot(2, 1, 1), plot(real(h), '-b.'), hold on, plot(real(hest), ':r.'), hold off, grid minor, legend('channel', 'test'), title('real')
subplot(2, 1, 2), plot(imag(h), '-b.'), hold on, plot(imag(hest), ':r.'), hold off, grid minor, legend('channel', 'test'), title('imag')

% % check xcorr
% a1 = xcorr(Ga, 1i*Ga);
% a2 = xcorr(1i*Ga, Ga);
% a3 = 1i*xcorr(Ga, Ga); % % -a1 = a2 = a3
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check xcorr
% a1 = xcorr(Ga, 1i*Gb);
% a2 = xcorr(1i*Ga, Gb);
% a3 = 1i*xcorr(Ga, Gb); % -a1 = a2 = a3
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check xcorr
% a1 = xcorr(Ga+1i*Gb, Ga+Gb);
% a2 = xcorr(Ga, Ga)+xcorr(Ga, Gb)+xcorr(1i*Gb, Ga)+xcorr(1i*Gb, Gb);
% a3 = xcorr(Ga, Ga)+xcorr(Ga, Gb)+1i*xcorr(Gb, Ga)+1i*xcorr(Gb, Gb); % a1 = a2 = a3
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check xcorr
% a1 = xcorr(Ga+Gb, Ga+1i*Gb);
% a2 = xcorr(Ga, Ga)+xcorr(Ga, 1i*Gb)+xcorr(Gb, Ga)+xcorr(Gb, 1i*Gb);
% a3 = xcorr(Ga, Ga)-1i*xcorr(Ga, Gb)+xcorr(Gb, Ga)-1i*xcorr(Gb, Gb); % a1 = a2 = a3
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check xcorr
% a1 = xcorr(Ga+1i*Gb, Ga+1i*Gb);
% a2 = xcorr(Ga, Ga)+xcorr(Ga, 1i*Gb)+xcorr(1i*Gb, Ga)+xcorr(1i*Gb, 1i*Gb);
% a3 = xcorr(Ga, Ga)-1i*xcorr(Ga, Gb)+1i*xcorr(Gb, Ga)+xcorr(Gb, Gb); % a1 = a2 = a3
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check xcorr
% a1 = xcorr(Ga+1i*Gb, 1i*Ga+Gb);
% a2 = xcorr(Ga, 1i*Ga)+xcorr(Ga, Gb)+xcorr(1i*Gb, 1i*Ga)+xcorr(1i*Gb, Gb);
% a3 = -1i*xcorr(Ga, Ga)+xcorr(Ga, Gb)+xcorr(Gb, Ga)+1i*xcorr(Gb, Gb); % a1 = a2 = a3
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check xcorr
% a1 = xcorr(1i*Ga+1*Gb, 1*Ga-1i*Gb);
% a2 = xcorr(1i*Ga, 1*Ga)+xcorr(1i*Ga, -1i*Gb)+xcorr(1*Gb, 1*Ga)+xcorr(1*Gb, -1i*Gb);
% a3 = 1i*xcorr(Ga, Ga)-xcorr(Ga, Gb)+xcorr(Gb, Ga)+1i*xcorr(Gb, Gb); % a1 = a2 = a3
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check xcorr
% a1 = xcorr(Ga, Gb);
% a2 = flipud(xcorr(Gb, Ga));
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), hold off, grid minor, title('imag')

% % check difference between conv and xcorr
% a1 = xcorr(Ga, Gb);
% a2 = flipud(conv(flipud(Ga), Gb));
% a3 = flipud(conv(Gb, flipud(Ga)));
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check difference between conv and xcorr
% a1 = xcorr(Ga, Gb);
% a2 = conv(Ga, flipud(Gb));
% a3 = conv(flipud(Gb), Ga);
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check difference between conv and xcorr
% a1 = xcorr(Ga, Gb);
% a2 = xcorr(flipud(Gb), flipud(Ga));
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), hold off, grid minor, title('imag')

% % check difference between conv and xcorr
% a1 = xcorr(Ga, flipud(Gb));
% a2 = xcorr(Gb, flipud(Ga));
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), hold off, grid minor, title('imag')

% % check difference between conv and xcorr
% a1 = conv(Ga, Gb);
% a2 = conv(flipud(Ga), flipud(Gb));
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), hold off, grid minor, title('imag')

% % check difference between conv and xcorr
% a1 = conv(flipud(Ga), Gb);
% a2 = conv(Ga, flipud(Gb));
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), hold off, grid minor, title('imag')

% % check difference between conv and xcorr
% a1 = xcorr(Ga, Ga)+xcorr(Gb, Gb);
% a2 = conv(Ga, flipud(Ga))+conv(Gb, flipud(Gb));
% a3 = conv(flipud(Ga), Ga)+conv(flipud(Gb), Gb);
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check conv
% a1 = conv(Ga, flipud(1i*Gb));
% a2 = 1i*conv(Ga, flipud(Gb));
% a3 = -1i*conv(Ga, flipud(Gb)); % a1 = a2 = -a3
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check conv
% a1 = conv(1i*Ga+1*Gb, flipud(conj(1*Ga-1i*Gb)));
% a2 = conv(1i*Ga, flipud(1*Ga))+conv(1i*Ga, flipud(conj(-1i*Gb)))+conv(1*Gb, flipud(1*Ga))+conv(1*Gb, flipud(conj(-1i*Gb)));
% a3 = 1i*conv(Ga, flipud(Ga))-conv(Ga, flipud(Gb))+conv(Gb, flipud(Ga))+1i*conv(Gb, flipud(Gb));
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')

% % check xcorr
% a1 = xcorr(1i*Ga+1*Gb, 1*Ga-1i*Gb);
% a2 = conv(1i*Ga+1*Gb, flipud(conj(1*Ga-1i*Gb)));
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), hold off, grid minor, title('imag')

% % check xcorr
% a1 = conv(flipud(1i*Ga)+1*Gb, 1*Ga+flipud(conj(-1i*Gb)));
% a2 = conv(flipud(1i*Ga), 1*Ga)+conv(flipud(1i*Ga), flipud(conj(-1i*Gb)))+conv(1*Gb, 1*Ga)+conv(1*Gb, flipud(conj(-1i*Gb)));
% a3 = 1i*conv(flipud(Ga), Ga)-conv(flipud(Ga), flipud(Gb))+conv(Gb, Ga)+1i*conv(Gb, flipud(Gb));
% figure
% subplot(2, 1, 1), plot(real(a1), '-bo'), hold on, plot(real(a2), '--rx'), plot(real(a3), ':k.'), hold off, grid minor, title('real')
% subplot(2, 1, 2), plot(imag(a1), '-bo'), hold on, plot(imag(a2), '--rx'), plot(imag(a3), ':k.'), hold off, grid minor, title('imag')





% [Ga,Gb] = wlanGolaySequence(32);
% 
% h = complex(randn(16, 1), randn(16, 1));
% y = conv(h, (1*Ga)+(-1i*Gb));
% 
% 
% 
% z = conv(y, (-1i*Ga)+(1*Gb));
% w = -1i*conv(h, conv(Ga, Ga)+conv(Gb, Gb));
% 
% hh = -20i*h;
% 
% 
% figure, subplot(2, 1, 1), plot(real(z), '-b.'), hold on, plot(real(w), ':r.'), plot(real(hh), ':ko'), hold off, grid minor, legend('practical', 'theoretical'), title('real'), subplot(2, 1, 2), plot(imag(z), '-b.'), hold on, plot(imag(w), ':r.'), plot(imag(hh), ':ko'), hold off, grid minor, legend('practical', 'theoretical', 'real channel'), title('imag')
% 
% 
% a = conv(h, flipud(Ga));
% b = xcorr(h, Ga);
% c = conv(Ga, Ga)+conv(Gb, Gb);
% d = xcorr(Ga, Ga)+xcorr(Gb, Gb);
% 
% 
% 
% figure, subplot(2, 1, 1), plot(real(a), '-b.'), hold on, plot(real(b), ':r.'), hold off, grid minor, title('real'), subplot(2, 1, 2), plot(imag(a), '-b.'), hold on, plot(imag(b), ':r.'), hold off, grid minor, title('imag')
% figure, subplot(2, 1, 1), plot(real(c), '-b.'), hold on, plot(real(d), ':r.'), hold off, grid minor, title('real'), subplot(2, 1, 2), plot(imag(c), '-b.'), hold on, plot(imag(d), ':r.'), hold off, grid minor, title('imag')
