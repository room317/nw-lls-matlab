% rate match for turbo code

function Zo = tx_ratematch_r2(u, rm)

num_cb = size(u, 2);        % number of code block

x1 = [zeros(rm.nDummy, num_cb); u(1 : 3 : end, :)];
x2 = [zeros(rm.nDummy, num_cb); u(2 : 3 : end, :)];
x3 = [zeros(rm.nDummy, num_cb); u(3 : 3 : end, :)];

% subblock interleaving
y1 = x1(rm.subblk_int_tbl+1, :);
y2 = x2(rm.subblk_int_tbl+1, :);
y3 = x3(rm.subblk_int_tbl+1, :);
y23 = permute([y2 y3], [2 1 3]);

% bit Collection and selection
Yo = squeeze([y1; reshape(y23, [], 1, num_cb)]);
k0 = rm.nR * 2;             % redundancy version number = 0

len_unpunc_tbl = size(rm.punc_tbl(~rm.punc_tbl(:, 1)), 1);
num_buff = ceil(rm.E(1)/len_unpunc_tbl);

punc_tbl_circ = circshift(rm.punc_tbl, -k0);
Yo_circ = circshift(Yo, -k0);
Zo_circ = Yo_circ(~punc_tbl_circ, :);
Zo_rep = repmat(Zo_circ, num_buff, 1);

Zo = Zo_rep(1:rm.E(1), :);

end
