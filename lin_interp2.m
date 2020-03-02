% 2-d matrix linear interpolation
% xp = lin_interp(xn, tn, tp, dim) interpolates xn
% xn: signal before interpolation
% tn: input signal indices (vector) on the axis before interpolation
% tp: output signal indices (vector) on the axis after interpolation
% dim: dimension of interpolation

function xp = lin_interp2(xn, tn, tp, dim)

% shift dimension
if dim ~= 2
    xn_shiftdim = xn.';
else
    xn_shiftdim = xn;
end

xn_slp = (xn_shiftdim(:, end) - xn_shiftdim(:, 1)) / (tn(1, end) - tn(1, 1));
xn_ntrcpt = xn_shiftdim(:, 1) - xn_slp * tn(1, 1);
xp_shiftdim = xn_slp * tp + repmat(xn_ntrcpt, 1, length(tp));

% restore dimension
if dim ~= 2
    xp = xp_shiftdim.';
else
    xp = xp_shiftdim;
end

% test result
% for i = 1 : size(xn_shiftdim, 1)
%     subplot(2, 1, 1)
%     figure(1), plot(tn,real(xn(i, :)),'o',tp,real(xp(i, :)),'-r.')
%     subplot(2, 1, 2)
%     figure(1), plot(tn,imag(xn(i, :)),'o',tp,imag(xp(i, :)),'-r.')
%     pause
% end

end
