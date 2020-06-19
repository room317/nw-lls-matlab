function u = rx_ratematch_r2(Zo, rm)

num_cb = size(Zo, 2);       % number of code block

% average
len_unpunc_tbl = size(rm.punc_tbl(~rm.punc_tbl(:, 1)), 1);
num_buff = ceil(size(Zo, 1)/len_unpunc_tbl);
len_tail = len_unpunc_tbl-mod(size(Zo, 1), len_unpunc_tbl);
Zo_tail = [Zo; zeros(len_tail, num_cb)];

Zo_sum = squeeze(sum(reshape(Zo_tail, len_unpunc_tbl, num_buff, num_cb), 2));
num_buff_vec = [num_buff*ones(len_unpunc_tbl-len_tail, num_cb); (num_buff-1)*ones(len_tail, num_cb)];
Zo_avg = Zo_sum./num_buff_vec;

% bit collection and selection
Ncb = rm.Kw;
k0 = rm.nR * 2;

Yo = zeros(Ncb, num_cb);
punc_tbl_circ = circshift(rm.punc_tbl, -k0);
Yo(~punc_tbl_circ, :) = Zo_avg;
Yo_circ = circshift(Yo, k0);

% subblock deinterleaving
y1 = Yo_circ(1:rm.Kpi, :);
Z12 = reshape(y1, rm.nR, [], num_cb);
Z12_perm = permute(Z12, [2 1 3]);
Z11 = zeros(size(Z12_perm));
Z11(rm.P+1, :, :) = Z12_perm;
Z11_reshape = reshape(Z11, [], 1, num_cb);
y23 = Yo_circ(rm.Kpi+1:end, :);
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
