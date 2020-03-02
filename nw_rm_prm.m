function nw_rm = nw_rm_prm(len_tb_bit, mcs, num, cc)

% mcs table (Qm, code_rate_x1024)
mcs_tbl = [ ...
    2  78; 2 120; 2 193; 2 308; 2 449; 2 602; 4 378; 4 490; ...
    4 616; 6 466; 6 567; 6 666; 6 772; 6 873; 6 948];

% fixed parameters
Qm = mcs_tbl(mcs, 1);                   % bits per symbol
code_rate = mcs_tbl(mcs, 2) / 1024;     % code rate
efficiency = Qm * code_rate;            % efficiency
NL = 1;                                 % number of layers

% rate matching parameter
B = len_tb_bit;                 % length of transport block (bit)
D = cc.K + 4;
nC = 32;                        % number of column
nR = ceil(D / nC);              % number of rows
Kpi = nC * nR;                  % block interleaver length
Kw = Kpi * 3;                   % circular buffer length
nDummy = Kpi - D;
P = bitrevorder(0 : nC - 1);    % subblock permutation pattern

% rate matching output sequence length per coded block
len_tb_sym_req = ceil(B / efficiency);                              % number of req. symbols
num_ub_req = ceil(len_tb_sym_req / (num.len_rb_sym_user * cc.C));   % number of req. user blocks (1 user block per 1 cb)
G = num_ub_req * (num.len_rb_sym_user * cc.C) * Qm;                 % number of tx bits required
G_ = G / (NL * Qm);                                                 % number of tx symbols per layer required
gamma = mod(G_, cc.C);
E = NL * Qm * [floor(G_ / cc.C) * ones(cc.C - gamma, 1); ceil(G_ / cc.C) * ones(gamma, 1)];

% simulation setup
num_data_sym_per_cw = E(1) / Qm;                            % number of symbol per codeword
num_subframe = num_data_sym_per_cw / num.len_rb_sym_user;   % number of required subframes
t_tb = num_subframe * num.t_subframe;                       % transfer block length (sec)

% subblock interleaving table setup
seq_nR = (0 : nR - 1)';
tmp_blk = mod(repmat(seq_nR * nC, 1, nC) + repmat(P + 1, nR, 1), Kpi);
subblk_int_tbl = tmp_blk(:);

% puncturing table initialization
t1 = [ones(nDummy, 1); zeros(D, 1)];
t2 = [ones(nDummy, 1); zeros(D, 1)];
t3 = [ones(nDummy, 1); zeros(D, 1)];

% puncturing table subblock interleaving
Z1 = reshape(t1, nC, [])';
Z2 = Z1(:, P + 1);
s1 = Z2(:);
Z1 = reshape(t2, nC, [])';
Z2 = Z1(:, P + 1);
s2 = Z2(:);
s3 = t3(subblk_int_tbl + 1);
s23 = [s2 s3]';
punc_tbl = [s1; s23(:)];

% output
nw_rm.Qm = Qm;            % bits per symbol
nw_rm.D = D;              % rate matching input length (bit)
nw_rm.nC = nC;            % number of column
nw_rm.nR = nR;            % number of rows
nw_rm.Kpi = Kpi;          % block interleaver length
nw_rm.Kw = Kw;            % circular buffer length
nw_rm.nDummy = nDummy;
nw_rm.P = P;
nw_rm.E = E;              % rate matching output sequence length per coded block
nw_rm.num_data_sym_per_cw = num_data_sym_per_cw;
nw_rm.t_tb = t_tb;
nw_rm.subblk_int_tbl = subblk_int_tbl;    % subblock interleaving table
nw_rm.punc_tbl = punc_tbl;                % puncturing table

end
