% rate match for turbo code
% idx_cb must be scalar
% u must be vector when idx_cb is not empty

function Zo = tx_ratematch_r2(u, idx_cb, rm)

% set the number of input code block
num_cb = size(u, 2);    % support multi-channel

x1 = [zeros(rm.nDummy, num_cb); u(1:3:end, :)];
x2 = [zeros(rm.nDummy, num_cb); u(2:3:end, :)];
x3 = [zeros(rm.nDummy, num_cb); u(3:3:end, :)];

% subblock interleaving
y1 = x1(rm.subblk_int_tbl+1, :);
y2 = x2(rm.subblk_int_tbl+1, :);
y3 = x3(rm.subblk_int_tbl+1, :);
y23 = permute([y2 y3], [2 1 3]);

% bit collection and selection
Yo = squeeze([y1; reshape(y23, [], 1, num_cb)]);
k0 = rm.nR * 2;             % redundancy version number = 0

len_unpunc_tbl = size(rm.punc_tbl(~rm.punc_tbl(:, 1)), 1);
num_buff = ceil(rm.E_true(idx_cb)/len_unpunc_tbl);       % vector (per cb)

punc_tbl_circ = circshift(rm.punc_tbl, -k0);
Yo_circ = circshift(Yo, -k0);
Zo_circ = Yo_circ(~punc_tbl_circ, :);
Zo_rep = repmat(Zo_circ, num_buff, 1);

Zo = Zo_rep(1:rm.E_true(idx_cb), :);

% % dump
% assignin('base', 'u', u)
% assignin('base', 'x1', x1)
% assignin('base', 'x2', x2)
% assignin('base', 'x3', x3)
% assignin('base', 'y1', y1)
% assignin('base', 'y2', y2)
% assignin('base', 'y3', y3)
% assignin('base', 'y23', y23)
% assignin('base', 'Yo', Yo)
% assignin('base', 'Yo_circ', Yo_circ)
% assignin('base', 'Zo_circ', Zo_circ)
% assignin('base', 'Zo_rep', Zo_rep)
% assignin('base', 'Zo', Zo)

end
