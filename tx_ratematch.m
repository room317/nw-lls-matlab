% rate match for turbo code

function Zo = tx_ratematch(u, rm)

x1 = [zeros(rm.nDummy, 1); squeeze(u(1 : 3 : end, 1, :))];
x2 = [zeros(rm.nDummy, 1); squeeze(u(2 : 3 : end, 1, :))];
x3 = [zeros(rm.nDummy, 1); squeeze(u(3 : 3 : end, 1, :))];

% subblock interleaving
Z1 = reshape(x1, rm.nC, [])';    
Z2 = Z1(:, rm.P + 1);
y1 = Z2(:);
Z1 = reshape(x2, rm.nC, [])';
Z2 = Z1(:, rm.P + 1);
y2 = Z2(:);
y3 = x3(rm.subblk_int_tbl + 1);
y23 = [y2 y3]';

% bit Collection and selection
Yo = [y1; y23(:)];
Ncb = rm.Kw; % for uplink
k0 = rm.nR * 2; % redundancy version number = 0
kk = 0;
Zo = zeros(rm.E(1), 1);

for kj = 0 : rm.E(1) - 1
    idx = mod(k0 + kj, Ncb) + 1;
    if rm.punc_tbl(idx) == 0
        Zo(kk + 1) = Yo(idx);
        kk = kk + 1;
    end
    if kk == rm.E(1)
        break;
    end
end

end
