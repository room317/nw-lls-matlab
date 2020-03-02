function y = rx_ratematch(Zo, rm)

% bit collection and selection
Yo = zeros(rm.Kw, 1);
Ncb = rm.Kw;
k0 = rm.nR * 2;
kk = 0;
for kj = 0 : rm.E(1) - 1
    idx = mod(k0 + kj, Ncb) + 1;
    if rm.punc_tbl(idx) == 0
        Yo(idx) = Zo(kk + 1);
%         Yo(idx) = Yo(idx) + Zo(kk + 1);
        kk = kk + 1;
    end
    if kk == rm.E(1)
        break;
    end
end

u = Yo;
u1 = u(1 : rm.Kpi);
U = reshape(u(rm.Kpi + 1 : rm.Kw), 2, [])';
u2 = U(:, 1);
u3 = U(:, 2);

Z1 = reshape(u1, rm.nR, []);
Z1_perm = zeros(size(Z1));
Z1_perm(:, rm.P + 1) = Z1;
Z1_perm_tr = Z1_perm';
x1 = Z1_perm_tr(:);

Z2 = reshape(u2, rm.nR, []);
Z2_perm = zeros(size(Z2));
Z2_perm(:, rm.P + 1) = Z2;
Z2_perm_tr = Z2_perm';
x2 = Z2_perm_tr(:);

x3 = zeros(rm.Kpi, 1);
x3(rm.subblk_int_tbl + 1) = u3;

y1 = x1(rm.nDummy + 1 : rm.nDummy + rm.D, 1);
y2 = x2(rm.nDummy + 1 : rm.nDummy + rm.D, 1);
y3 = x3(rm.nDummy + 1 : rm.nDummy + rm.D, 1);

y_tmp = [y1 y2 y3].';
y = y_tmp(:);

end
