% ref: 3gpp ts 36.212 (36212-e40)
% cell based channel coding parameters

function nw_rm = nw_rm_prm_r1(mcs, cc, num, sim_option)

% mcs table (Qm, code_rate_x1024)
%     mcs 16 (for test): 16QAM code rate 2/3
%     mcs 17 (for test): 16QAM code rate 3/4
mcs_tbl = [ ...
    2  78; 2 120; 2 193; 2 308; 2 449; 2 602; 4 378; 4 490; ...
    4 616; 6 466; 6 567; 6 666; 6 772; 6 873; 6 948;
    4 666; 4 772];

% primary parameters
if sim_option.override
    Q_m = sim_option.Qm;
    coderate = sim_option.coderate/1024;
else
    if mcs > size(mcs_tbl, 1)
        error('MCS setting example: MCS1 ~ MCS15')
    else
        Q_m = mcs_tbl(mcs, 1);                  % bits per symbol
        coderate = mcs_tbl(mcs, 2)/1024;        % code rate
    end
end

switch Q_m
    case 2
        Q_mod = 'QPSK';
        lteSymMap = 'gray';
    case 4
        Q_mod = '16QAM';
        lteSymMap = [11 10 14 15 ... 
                      9  8 12 13 ...
                      1  0  4  5 ...
                      3  2  6  7];
    case 6
        Q_mod = '64QAM';
        lteSymMap = 'gray';
    otherwise
        error('Wrong modulation order.')
end

% fixed parameters
efficiency = Q_m*coderate;           % efficiency (number of bits per qam symbol)

% sub-block interleaver
C_TC = 32;                          % number of columns of the matrix
R_TC = ceil(cc.D_r/C_TC);           % number of rows of the matrix
K_PI = R_TC*C_TC;                   % size of sub-block interleaver
N_D = K_PI-cc.D_r;                  % number of dummy bits
P = bitrevorder(0:C_TC-1)';          % inter-column permutation pattern for sub-block interleaver

Y_P = cell(1, cc.C);
PI_k = cell(1, cc.C);
for r = 1:cc.C
    Y_k = reshape(1:K_PI(r), C_TC, [])';        % sub-block interleaver matrix
    Y_P{r} = reshape(Y_k(:, P+1), [], 1);       % inter-column permuted matrix output index
    PI_k{r} = mod(P(floor((0:K_PI(r)-1)'./R_TC(r))+1)+C_TC*mod((0:K_PI(r)-1)', R_TC(r))+1, K_PI(r))+1;  % output index of the sub-block interleaver
end

% bit collection, selection and transmission
K_w = 3*K_PI;                       % circular buffer length
N_cb = K_w;                         % soft buffer size for r-th code block (for uplink)

N_RE = num.num_qamsym_usr;          % number of res per rb
N_RB = ceil(cc.A/efficiency/N_RE);  % number of required rbs
G = N_RB*N_RE*Q_m;                  % total number of bits available for the transmission of one transport block
N_L = 1;                            % number of layers a transport block is mapped onto
G_ = G/(N_L*Q_m);                   % total number of qam symbols per layer available for the transmission of one transport block
GAMMA_ = mod(G_, cc.C);             % number of remainder qam symbols for code block allocation

% rate matching output sequence length for r-th code block
E = zeros(cc.C, 1);
for r = 1:cc.C
    if r <= cc.C-GAMMA_
        E(r) = N_L*Q_m*floor(G_/cc.C);
    else
        E(r) = N_L*Q_m*ceil(G_/cc.C);
    end
end

% output
nw_rm.Q_m = Q_m;                    % number of bits per qam symbol
nw_rm.Q_mod = Q_mod;                % qam modulation ('QPSK', '16QAM', '64QAM')
nw_rm.lteSymMap = lteSymMap;        % qam symbol map
nw_rm.coderate = coderate;          % channel code rate
nw_rm.efficiency = efficiency;      % efficiency (number of information bits per qam symbol)
nw_rm.C_TC = C_TC;                  % number of columns of the matrix
nw_rm.R_TC = R_TC;                  % number of rows of the matrix
nw_rm.K_PI = K_PI;                  % size of sub-block interleaver
nw_rm.N_D = N_D;                    % number of dummy bits
nw_rm.P = P;                        % inter-column permutation pattern for sub-block interleaver
nw_rm.Y_P = Y_P;                    % output index of inter-column permuted matrix
nw_rm.PI_k = PI_k;                  % output index of the sub-block interleaver
nw_rm.K_w = K_w;                    % circular buffer length
nw_rm.N_cb = N_cb;                  % soft buffer size for r-th code block (for uplink)
nw_rm.G = G;                        % total number of bits available for the transmission of one transport block
nw_rm.N_L = N_L;                    % number of layers a transport block is mapped onto
nw_rm.G_ = G_;                      % total number of qam symbols per layer available for the transmission of one transport block
nw_rm.GAMMA_ = GAMMA_;              % number of remainder qam symbols for code block allocation
nw_rm.E = E;                        % rate matching output sequence length for r-th code block
nw_rm.N_RB = N_RB;                  % number of rbs

end
