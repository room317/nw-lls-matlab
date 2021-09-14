% ref: 3gpp ts 36.212 (36212-e40)

function nw_rm = nw_rm_prm(len_tb_bit, mcs, num, cc, sim_option, test_option)

% mcs table (Qm, code_rate_x1024)
%     mcs 16 (for test): 16QAM code rate 2/3
%     mcs 17 (for test): 16QAM code rate 3/4
mcs_tbl = [ ...
    2  78; 2 120; 2 193; 2 308; 2 449; 2 602; 4 378; 4 490; ...
    4 616; 6 466; 6 567; 6 666; 6 772; 6 873; 6 948;
    4 666; 4 772];

% primary parameters
if sim_option.override
    Qm = sim_option.Qm;
    coderate = sim_option.coderate/1024;
else
    if mcs > size(mcs_tbl, 1)
        error('MCS setting example: MCS1 ~ MCS15')
    else
        Qm = mcs_tbl(mcs, 1);                   % bits per symbol
        coderate = mcs_tbl(mcs, 2)/1024;        % code rate
    end
end

% fixed parameters
efficiency = Qm*coderate;               % efficiency (number of bits per qam symbol)
NL = 1;                                 % number of layers

% rate matching parameter
B = len_tb_bit;                 % length of transport block (bit)
D = cc.K+4;
nC = 32;                        % number of column
nR = ceil(D/nC);                % number of rows
Kpi = nC*nR;                    % block interleaver length
Kw = Kpi*3;                     % circular buffer length
nDummy = Kpi-D;
P = bitrevorder(0:nC-1);        % subblock permutation pattern

% rate matching output sequence length per coded block
%   - single physical rb doesn't have to contain only one cb data (assumption)
%   - E: last gamma cbs have residual qam symbols
%   - G_=E(1)+E(2)+E(3)+...
G = ceil(B/efficiency)*Qm;      % number of tx bits required (the total number of bits available for 1 tb)
G_ = ceil(G/(NL*Qm));           % number of tx symbols per layer required
gamma = mod(G_, cc.C);          % number of residual symbols
E = NL*Qm*[floor(G_/cc.C)*ones(cc.C-gamma, 1); ceil(G_/cc.C)*ones(gamma, 1)];       % number of total tx bits per code block (vector)

% simulation setup
num_data_sym_cb = E/Qm;                                     % number of symbol per code block (vector)
if test_option.pilot_only
    num_usrfrm_cb = 1;
else
    num_usrfrm_cb = ceil(num_data_sym_cb/num.num_qamsym_usr);   % number of required user frames (previously num_slot, number of slots, vector)
end
t_tb = sum(num_usrfrm_cb)*num.t_usrfrm;                     % transfer block length (sec)
E_true = num.num_qamsym_usr*num_usrfrm_cb*Qm;

% subblock interleaving table setup
seq_nR = (0 : nR-1)';
tmp_blk = mod(repmat(seq_nR*nC, 1, nC) + repmat(P+1, nR, 1), Kpi);
subblk_int_tbl = tmp_blk(:);

% puncturing table initialization
t1 = [ones(nDummy, 1); zeros(D, 1)];
t2 = [ones(nDummy, 1); zeros(D, 1)];
t3 = [ones(nDummy, 1); zeros(D, 1)];

% puncturing table subblock interleaving
s1 = t1(subblk_int_tbl + 1, :);
s2 = t2(subblk_int_tbl + 1, :);
s3 = t3(subblk_int_tbl + 1, :);
s23 = [s2 s3]';
punc_tbl = [s1; s23(:)];

% output
nw_rm.B = B;                % transmit block length
nw_rm.Qm = Qm;              % bits per symbol
nw_rm.D = D;                % rate matching input length (bit)
nw_rm.nC = nC;              % number of column
nw_rm.nR = nR;              % number of rows
nw_rm.Kpi = Kpi;            % block interleaver length
nw_rm.Kw = Kw;              % circular buffer length
nw_rm.nDummy = nDummy;
nw_rm.P = P;
nw_rm.E = E;                % rate matching output sequence length per coded block (number of total tx bits per code block)
nw_rm.E_true = E_true;      % rate matching output sequence length per coded block (number of total tx bits per code block)
nw_rm.len_tb_sym_req = G;   % number of req. symbols
nw_rm.num_data_sym_cb = num_data_sym_cb;
nw_rm.num_usrfrm_cb = num_usrfrm_cb;              % number of user frame
nw_rm.t_tb = t_tb;
nw_rm.subblk_int_tbl = subblk_int_tbl;      % subblock interleaving table
nw_rm.punc_tbl = punc_tbl;                  % puncturing table

% assignin('base', 'coderate', coderate)
% assignin('base', 'Qm', Qm)
% assignin('base', 'efficiency', efficiency)
% assignin('base', 'B', B)
% assignin('base', 'nC', nC)
% assignin('base', 'nR', nR)
% assignin('base', 'Kpi', Kpi)
% assignin('base', 'Kw', Kw)
% assignin('base', 'nDummy', nDummy)
% assignin('base', 'P', P)
% assignin('base', 'G', G)
% assignin('base', 'G_', G_)
% assignin('base', 'gamma', gamma)
% assignin('base', 'E', E)
% assignin('base', 'num_data_sym_cb', num_data_sym_cb)
% assignin('base', 'num_usrfrm_cb', num_usrfrm_cb)
% assignin('base', 't_tb', t_tb)
% assignin('base', 'E_true', E_true)
% pause

end
