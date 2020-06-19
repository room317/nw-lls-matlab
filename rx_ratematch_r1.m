function u = rx_ratematch_r1(Zo, rm)

num_cb = size(Zo, 2);       % number of code block

% bit collection and selection
Ncb = rm.Kw;
k0 = rm.nR * 2;

num_punc_tbl = ceil((rm.E(1)+k0)/Ncb);
len_punc_tbl_tail = (rm.E(1)+k0)-(Ncb*(num_punc_tbl-1));
punc_tbl_circ = [rm.punc_tbl(k0+1:end, 1); repmat(rm.punc_tbl, num_punc_tbl-2, 1); rm.punc_tbl(1:len_punc_tbl_tail, 1)];

Yo_circ = zeros(rm.E(1), num_cb);
Yo_circ(~punc_tbl_circ, :) = Zo(1:size(punc_tbl_circ(~punc_tbl_circ), 1), :);
Yo = Yo_circ(Ncb-k0+1:2*Ncb-k0, :);

% subblock deinterleaving
y1 = Yo(1:rm.Kpi, :);
Z12 = reshape(y1, rm.nR, [], num_cb);
Z12_perm = permute(Z12, [2 1 3]);
Z11 = zeros(size(Z12_perm));
Z11(rm.P+1, :, :) = Z12_perm;
Z11_reshape = reshape(Z11, [], 1, num_cb);
y23 = Yo(rm.Kpi+1:end, :);
y23_reshape = reshape(y23, 2, [], num_cb);
y23_perm = permute(y23_reshape, [2 1 3]);
y2 = y23_perm(:, 1, :);
Z22 = reshape(y2, rm.nR, [], num_cb);
Z22_perm = permute(Z22, [2 1 3]);
Z21 = zeros(size(Z22_perm));
Z21(rm.P+1, :, :) = Z22_perm;
Z21_reshape = reshape(Z21, [], 1, num_cb);
y3 = y23_perm(:, 2, :);
y3_intlvr = zeros(rm.Kpi, 1, num_cb);
y3_intlvr(rm.subblk_int_tbl+1, :, :) = y3;

x1 = Z11_reshape(rm.nDummy+1:end, :, :);
x2 = Z21_reshape(rm.nDummy+1:end, :, :);
x3 = y3_intlvr(rm.nDummy+1:end, :, :);

u = squeeze(reshape(permute([x1 x2 x3], [2 1 3]), [], 1, num_cb));

end
