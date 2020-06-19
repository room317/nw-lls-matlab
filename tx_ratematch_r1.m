% rate match for turbo code

function Zo_pad = tx_ratematch_r1(u, rm)

num_cb = size(u, 2);        % number of code block

x1 = [zeros(rm.nDummy, num_cb); u(1 : 3 : end, :)];
x2 = [zeros(rm.nDummy, num_cb); u(2 : 3 : end, :)];
x3 = [zeros(rm.nDummy, num_cb); u(3 : 3 : end, :)];

% subblock interleaving
Z11 = reshape(x1, rm.nC, [], num_cb);
Z12 = permute(Z11(rm.P + 1, :, :), [2 1 3]);
y1 = reshape(Z12, [], 1, num_cb);
Z21 = reshape(x2, rm.nC, [], num_cb);
Z22 = permute(Z21(rm.P + 1, :, :), [2 1 3]);
y2 = reshape(Z22, [], 1, num_cb);
y3 = reshape(x3(rm.subblk_int_tbl + 1, :), [], 1, num_cb);
y23 = permute([y2 y3], [2 1 3]);

% bit Collection and selection
Yo = squeeze([y1; reshape(y23, [], 1, num_cb)]);
Ncb = rm.Kw;                % for uplink
k0 = rm.nR * 2;             % redundancy version number = 0

num_punc_tbl = ceil((rm.E(1)+k0)/Ncb);
len_punc_tbl_tail = (rm.E(1)+k0)-(Ncb*(num_punc_tbl-1));
punc_tbl_circ = [rm.punc_tbl(k0+1:end, 1); repmat(rm.punc_tbl, num_punc_tbl-2, 1); rm.punc_tbl(1:len_punc_tbl_tail, 1)];
Yo_circ = [Yo(k0+1:end, :); repmat(Yo, num_punc_tbl-2, 1); Yo(1:len_punc_tbl_tail, :)];

Zo = Yo_circ(~punc_tbl_circ, :);

num_pad = rm.E(1)-size(Zo, 1);
Zo_pad = [Zo; zeros(num_pad, size(Zo, 2))];

end
