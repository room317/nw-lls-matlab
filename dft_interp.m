% dft interpolation
% xp = dft_interp(xn, p, dim, shift_pnt) interpolates xn
% xn: signal before interpolation
% p (>= size(x)): number of evenly spaced points after interpolation
% dim: dimension of interpolation
% shift_pnt: number of interpolated points to shift

function xp = dft_interp(xn, p, dim, shift_pnt)

% shift dimension
if dim ~= 2
    xn_shiftdim = xn.';
else
    xn_shiftdim = xn;
end

n = size(xn_shiftdim, 2);
yn = fft(xn_shiftdim, [], 2);
yp = zeros(size(yn, 1), p);
yp(:, 1:n/2+1) = yn(:, 1:n/2+1);
yp(:, p-n/2+2:p) = yn(:, n/2+2:n);
xp_shiftdim = ifft(yp, [], 2) * (p/n);

% circular shift
xp_shiftpnt = circshift(xp_shiftdim, shift_pnt, 2);

% restore dimension
if dim ~= 2
    xp = xp_shiftpnt.';
else
    xp = xp_shiftpnt;
end

% test result
% tn = (0 : n-1) / n;     % n evenly-spaced time points
% tp = (0 : p-1) / p;     % p evenly-spaced time points
% figure, plot(tn,real(circshift(xn, n*shift_pnt/p, dim)),'o',tp,real(xp),'-r.')

end
