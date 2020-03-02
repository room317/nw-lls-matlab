% linear interpolation
% xp = lin_interp(xn, n, p, dim) interpolates xn
% xn: signal before interpolation
% n: input signal indices (vector) on the axis after interpolation
% p: output signal indices (vector) on the axis after interpolation
% dim: dimension of interpolation

function xp = lin_interp(xn, tn, tp, dim)

% shift dimension
if dim ~= 2
    xn_shiftdim = xn.';
else
    xn_shiftdim = xn;
end

xp_shiftdim = zeros(size(xn_shiftdim, 1), length(tp));
for idx = 1 : size(xn_shiftdim, 1)
    xp_shiftdim(idx, :) = interp1(tn, xn_shiftdim(idx, :), tp, 'linear', 'extrap');
end

% restore dimension
if dim ~= 2
    xp = xp_shiftdim.';
else
    xp = xp_shiftdim;
end

% test result
% figure, plot(tn,real(xn(50, :)),'o',tp,real(xp(50, :)),'-r.')

end
