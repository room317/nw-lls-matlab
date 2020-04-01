% effect of freq band selection

% original: original signal
% lowband: partial selection of original
% lowpass: parts of original forced to zeros

original_f = exp(-1i*2*pi*(0:255)/64);
original_t = ifft(original_f)*sqrt(256);
for i = 1 : 64
    lowband_f = original_f(1+i:256-i);
    lowpass_f = zeros(256, 1); lowpass_f(1+i:256-i) = lowband_f;
    lowband_t = ifft(lowband_f)*sqrt(224);
    lowpass_t = ifft(lowpass_f)*sqrt(256);
    idx_lowband_t = linspace(0, 255, 256-2*i);
    
    figure(1)
    plot(0:255, real(original_t), '-k.'), hold on
    plot(idx_lowband_t, real(lowband_t), '-b.')
    plot(0:255, real(lowpass_t), '-r.'), hold off, grid minor
    pause
end
