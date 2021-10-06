% tc_prm (K, f1, f2), K: block size
tc_prm = [ ...
      40   3  10;      48   7  12;      56  19  42;      64   7  16; ...
      72   7  18;      80  11  20;      88   5  22;      96  11  24; ...
     104   7  26;     112  41  84;     120 103  90;     128  15  32; ...
     136   9  34;     144  17 108;     152   9  38;     160  21 120; ...
     168 101  84;     176  21  44;     184  57  46;     192  23  48; ...
     200  13  50;     208  27  52;     216  11  36;     224  27  56; ...
     232  85  58;     240  29  60;     248  33  62;     256  15  32; ...
     264  17 198;     272  33  68;     280 103 210;     288  19  36; ...
     296  19  74;     304  37  76;     312  19  78;     320  21 120; ...
     328  21  82;     336 115  84;     344 193  86;     352  21  44; ...
     360 133  90;     368  81  46;     376  45  94;     384  23  48; ...
     392 243  98;     400 151  40;     408 155 102;     416  25  52; ...
     424  51 106;     432  47  72;     440  91 110;     448  29 168; ...
     456  29 114;     464 247  58;     472  29 118;     480  89 180; ...
     488  91 122;     496 157  62;     504  55  84;     512  31  64; ...
     528  17  66;     544  35  68;     560 227 420;     576  65  96; ...
     592  19  74;     608  37  76;     624  41 234;     640  39  80; ...
     656 185  82;     672  43 252;     688  21  86;     704 155  44; ...
     720  79 120;     736 139  92;     752  23  94;     768 217  48; ...
     784  25  98;     800  17  80;     816 127 102;     832  25  52; ...
     848 239 106;     864  17  48;     880 137 110;     896 215 112; ...
     912  29 114;     928  15  58;     944 147 118;     960  29  60; ...
     976  59 122;     992  65 124;    1008  55  84;    1024  31  64; ...
    1056  17  66;    1088 171 204;    1120  67 140;    1152  35  72; ...
    1184  19  74;    1216  39  76;    1248  19  78;    1280 199 240; ...
    1312  21  82;    1344 211 252;    1376  21  86;    1408  43  88; ...
    1440 149  60;    1472  45  92;    1504  49 846;    1536  71  48; ...
    1568  13  28;    1600  17  80;    1632  25 102;    1664 183 104; ...
    1696  55 954;    1728 127  96;    1760  27 110;    1792  29 112; ...
    1824  29 114;    1856  57 116;    1888  45 354;    1920  31 120; ...
    1952  59 610;    1984 185 124;    2016 113 420;    2048  31  64; ...
    2112  17  66;    2176 171 136;    2240 209 420;    2304 253 216; ...
    2368 367 444;    2432 265 456;    2496 181 468;    2560  39  80; ...
    2624  27 164;    2688 127 504;    2752 143 172;    2816  43  88; ...
    2880  29 300;    2944  45  92;    3008 157 188;    3072  47  96; ...
    3136  13  28;    3200 111 240;    3264 443 204;    3328  51 104; ...
    3392  51 212;    3456 451 192;    3520 257 220;    3584  57 336; ...
    3648 313 228;    3712 271 232;    3776 179 236;    3840 331 120; ...
    3904 363 244;    3968 375 248;    4032 127 168;    4096  31  64; ...
    4160  33 130;    4224  43 264;    4288  33 134;    4352 477 408; ...
    4416  35 138;    4480 233 280;    4544 357 142;    4608 337 480; ...
    4672  37 146;    4736  71 444;    4800  71 120;    4864  37 152; ...
    4928  39 462;    4992 127 234;    5056  39 158;    5120  39  80; ...
    5184  31  96;    5248 113 902;    5312  41 166;    5376 251 336; ...
    5440  43 170;    5504  21  86;    5568  43 174;    5632  45 176; ...
    5696  45 178;    5760 161 120;    5824  89 182;    5888 323 184; ...
    5952  47 186;    6016  23  94;    6080  47 190;    6144 263 480];

% crc calculation
% caution: parity bits are padded after input bits
A = len_tb_bit;                         % size of the input sequence
gCRC24A = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];      % cyclic generator polynomial (gCRC24A)
gCRC24B = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];      % cyclic generator polynomial (gCRC24B)
L = length(gCRC24A)-1;                  % number of parity bits
B = A+L;                                % number of bits after crc attachment

% code block segmentation and code block crc attachment
Z = tc_prm(end, 1);                     % maximum code block size (Z = 6144)
Z_min = tc_prm(1, 1);                   % minimum code block size (40)

% perform segmentation of the input bit sequence
if B > Z
    % segmentation of the input bit sequence is performed
    % additional crc sequence of L = 24 bits is attached to each code block
    % caution: sequence with length K_r-L of code block number r is used to calculate the CRC parity bits with length 24.
    % caution: when crc calculation, assume filler bits have the value 0.
    % caution: additional crc sequence of L = 24 bits is attached to the end of each code block
    L_ = length(gCRC24A)-1;
    C = ceil(B/(Z-L_));     % total number of code blocks
    B_ = B+C*L_;
else
    L_ = 0;
    C = 1;                  % total number of code blocks
    B_ = B;
end

K_tmp = B_/C;
idx_K_plus = find(tc_prm(:, 1) >= K_tmp, 1, 'first');
K_plus = tc_prm(idx_K_plus, 1);                 % first segmentation size
f1_plus = tc_prm(idx_K_plus, 2);                % turbo code internal interleaver parameter
f2_plus = tc_prm(idx_K_plus, 3);                % turbo code internal interleaver parameter

% number of bits in each code block (C > 0)
if C > 1
    idx_K_minus = find(tc_prm(:, 1) < K_plus, 1, 'last');
    K_minus = tc_prm(idx_K_minus, 1);           % second segmentation size
    f1_minus = tc_prm(idx_K_minus, 2);          % turbo code internal interleaver parameter
    f2_minus = tc_prm(idx_K_minus, 3);          % turbo code internal interleaver parameter
    delta_K = K_plus - K_minus;
    C_minus = floor((C*K_plus-B_)/delta_K);     % number of segments of size K_minus
    C_plus = C-C_minus;                         % number of segments of size K_plus
else
    K_minus = 0;    % second segmentation size
    f1_minus = 0;   % turbo code internal interleaver parameter
    f2_minus = 0;   % turbo code internal interleaver parameter
    C_minus = 0;    % number of segments of size K_minus
    C_plus = 1;     % number of segments of size K_plus
end

% when F > 0, filler bits are added to the beginning of the first block
% when B < Z_min, filler bits are added to the beginning of the code block
F = C_plus*K_plus+C_minus*K_minus-B_;   % number of filler bits

% caution: C_minus first, C_plus last
K_r = [ones(C_minus, 1)*K_minus; ones(C_plus, 1)*K_plus];       % number of bits for the code block number r
f1_r = [ones(C_minus, 1)*f1_minus; ones(C_plus, 1)*f1_plus];    % number of bits for the code block number r
f2_r = [ones(C_minus, 1)*f2_minus; ones(C_plus, 1)*f2_plus];    % number of bits for the code block number r

% turbo code trellis structure
% (https://la.mathworks.com/help/comm/ref/turboencoder.html)
tc_trellis = poly2trellis(4, [13 15], 13);

D_r = K_r+4;                            % number of encoded bits per output stream

num_iter_max = 8;                       % number of maximum iteration

% turbo code internal interleaver
PI_r = cell(1, C);
for r = 1:C
    PI_r{r} = mod(f1_r(r)*(0:K_r(r)-1)'+f2_r(r)*(0:K_r(r)-1)'.^2, K_r(r))+1;
end

% output
cc.tc_trellis = tc_trellis;      % turbo code trellis structure
cc.gCRC24A = gCRC24A;            % cyclic generator polynomial (gCRC24A)
cc.gCRC24B = gCRC24B;            % cyclic generator polynomial (gCRC24B)
cc.num_iter_max = num_iter_max;  % number of maximum iteration
cc.A = A;                        % size of the input sequence
cc.L = L;                        % number of parity bits
cc.B = B;                        % number of bits after crc attachment
cc.Z = Z;                        % maximum code block size (Z = 6144)
cc.Z_min = Z_min;                % minimum code block size (40)
cc.L_ = L_;                      % additional crc sequence length attached to each code block (L_ = 24 bits)
cc.C = C;                        % total number of code blocks
cc.B_ = B_;                      % number of bits after additional crc attachment to each code block
cc.K_plus = K_plus;              % number of bits per code block (first segmentation size)
cc.f1_plus = f1_plus;            % turbo code internal interleaver parameter corresponds to K_plus
cc.f2_plus = f2_plus;            % turbo code internal interleaver parameter corresponds to K_plus
cc.C_plus = C_plus;              % number of segments of size K_plus
cc.K_minus = K_minus;            % number of bits per code block (second segmentation size)
cc.f1_minus = f1_minus;          % turbo code internal interleaver parameter corresponds to K_minus
cc.f2_minus = f2_minus;          % turbo code internal interleaver parameter corresponds to K_minus
cc.C_minus = C_minus;            % number of segments of size K_minus
cc.F = F;                        % number of filler bits
cc.K_r = K_r;                    % number of bits for the code block number r
cc.f1_r = f1_r;                  % turbo code internal interleaver parameter for the code block number r
cc.f2_r = f2_r;                  % turbo code internal interleaver parameter for the code block number r
cc.D_r = D_r;                    % number of encoded bits per output stream
cc.PI_r = PI_r;                  % turbo code internal interleaver








% generate bit stream
tx_bit = randi([0 1], A, 1);

% calculate crc
tx_crc_a = comm.CRCGenerator(gCRC24A);
tx_crc_b = comm.CRCGenerator(gCRC24B);
tx_bit_crc_ref = lteCRCEncode(tx_bit,'24A');
tx_bit_crc = tx_crc_a(tx_bit);

% segment code blocks and attach code block crcs
tx_bit_seg_ref = lteCodeBlockSegment(tx_bit_crc_ref);
tx_bit_seg = cell(1, C);
if C > 1
    for r = 1:C
        if r == 1
            % insert filler bits and calculate crc
            tx_bit_seg{r} = tx_crc_b([zeros(F, 1); tx_bit_crc(1:K_r(r)-L_-F)]);
        else
            % calculate crc
            tx_bit_seg{r} = tx_crc_b(tx_bit_crc(sum(K_r(1:r-1))-(r-1)*L_-F+1:sum(K_r(1:r))-r*L_-F));
        end
    end
else
    tx_bit_seg{1} = [zeros(F, 1); tx_bit_crc(1:K_r(r)-L_-F)];
end

% channel coding
tx_bit_enc_ref = lteTurboEncode(tx_bit_seg_ref);
tx_bit_enc = cell(1, C);
for r = 1:C
    % create turbo encoder/decoder objects
    turbo_enc = comm.TurboEncoder('TrellisStructure', tc_trellis, ...
        'InterleaverIndices', PI_r{r});
    
    % turbo encode data
    tx_bit_enc{r} = turbo_enc(tx_bit_seg{r});
end

% for r = 1:C
%     a = double(reshape(tx_bit_enc_ref{r}, [], 3));
%     b = reshape(tx_bit_enc{r}, 3, [])';
%     b(1:F, 1:2) = -1;
%     figure, plot(real(a(:, 1)), '-b.'), hold on, plot(real(b(:, 1)), ':r.'), hold off, grid minor
%     figure, plot(real(a(:, 2)), '-b.'), hold on, plot(real(b(:, 2)), ':r.'), hold off, grid minor
%     figure, plot(real(a(:, 3)), '-b.'), hold on, plot(real(b(:, 3)), ':r.'), hold off, grid minor
%     pause
% end









% mcs table (Qm, code_rate_x1024)
%     mcs 16 (for test): 16QAM code rate 2/3
%     mcs 17 (for test): 16QAM code rate 3/4
mcs_tbl = [ ...
    2  78; 2 120; 2 193; 2 308; 2 449; 2 602; 4 378; 4 490; ...
    4 616; 6 466; 6 567; 6 666; 6 772; 6 873; 6 948;
    4 666; 4 772];

mcs = 17;

% primary parameters
if mcs > size(mcs_tbl, 1)
    error('MCS setting example: MCS1 ~ MCS15')
else
    Q_m = mcs_tbl(mcs, 1);                   % bits per symbol
    coderate = mcs_tbl(mcs, 2)/1024;        % code rate
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
    Y_k = reshape(1:K_PI(r), C_TC, [])';    % sub-block interleaver matrix
    Y_P{r} = reshape(Y_k(:, P+1), [], 1);      % inter-column permuted matrix output index
    PI_k{r} = mod(P(floor((0:K_PI(r)-1)'./R_TC(r))+1)+C_TC*mod((0:K_PI(r)-1)', R_TC(r))+1, K_PI(r))+1;    % output index of the sub-block interleaver
end

% bit collection, selection and transmission
K_w = 3*K_PI;                       % circular buffer length
N_cb = K_w;                         % soft buffer size for r-th code block (for uplink)

G = ceil(A/efficiency)*Q_m;         % total number of bits available for the transmission of one transport block
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
rm.Qm = Q_m;                    % number of bits per qam symbol
rm.coderate = coderate;         % channel code rate
rm.efficiency = efficiency;     % efficiency (number of information bits per qam symbol)
rm.C_TC = C_TC;                 % number of columns of the matrix
rm.R_TC = R_TC;                 % number of rows of the matrix
rm.K_PI = K_PI;                 % size of sub-block interleaver
rm.N_D = N_D;                   % number of dummy bits
rm.P = P;                       % inter-column permutation pattern for sub-block interleaver
rm.Y_P = Y_P;                   % output index of inter-column permuted matrix
rm.PI_k = PI_k;                 % output index of the sub-block interleaver
rm.K_w = K_w;                   % circular buffer length
rm.N_cb = N_cb;                 % soft buffer size for r-th code block (for uplink)
rm.G = G;                       % total number of bits available for the transmission of one transport block
rm.N_L = N_L;                   % number of layers a transport block is mapped onto
rm.G_ = G_;                     % total number of qam symbols per layer available for the transmission of one transport block
rm.GAMMA_ = GAMMA_;             % number of remainder qam symbols for code block allocation
rm.E = E;                       % rate matching output sequence length for r-th code block





% tx_bit_enc_ref{1} = [(-1)*ones(F, 1); (F+N_D(1)+1:F+N_D(1)+D_r(1)-F)'; (-1)*ones(F, 1); 100+(F+N_D(1)+1:F+N_D(1)+D_r(1)-F)'; 200+(N_D(1)+1:N_D(1)+D_r(1))'];
% tx_bit_enc_ref{2} = [(-1)*ones(F, 1); 1000+(F+N_D(2)+1:F+N_D(2)+D_r(2)-F)'; (-1)*ones(F, 1); 1100+(F+N_D(2)+1:F+N_D(2)+D_r(2)-F)'; 1200+(N_D(2)+1:N_D(2)+D_r(2))'];
% % tx_bit_enc_ref{3} = [(-1)*ones(F, 1); 2000+(F+N_D(3)+1:F+N_D(3)+D_r(3)-F)'; (-1)*ones(F, 1); 2100+(F+N_D(3)+1:F+N_D(3)+D_r(3)-F)'; 2200+(N_D(3)+1:N_D(3)+D_r(3))'];
% % tx_bit_enc_ref{4} = [(-1)*ones(F, 1); 3000+(F+N_D(4)+1:F+N_D(4)+D_r(4)-F)'; (-1)*ones(F, 1); 3100+(F+N_D(4)+1:F+N_D(4)+D_r(4)-F)'; 3200+(N_D(4)+1:N_D(4)+D_r(4))'];
% 
% tx_bit_enc{1} = reshape([[(-1)*ones(1, F) (F+N_D(1)+1:F+N_D(1)+D_r(1)-F)]; [(-1)*ones(1, F) 100+(F+N_D(1)+1:F+N_D(1)+D_r(1)-F)]; 200+(N_D(1)+1:N_D(1)+D_r(1))], [], 1);
% tx_bit_enc{2} = reshape([[(-1)*ones(1, F) 1000+(F+N_D(2)+1:F+N_D(2)+D_r(2)-F)]; [(-1)*ones(1, F) 1100+(F+N_D(2)+1:F+N_D(2)+D_r(2)-F)]; 1200+(N_D(2)+1:N_D(2)+D_r(2))], [], 1);
% % tx_bit_enc{3} = reshape([[(-1)*ones(1, F) 2000+(F+N_D(3)+1:F+N_D(3)+D_r(3)-F)]; [(-1)*ones(1, F) 2100+(F+N_D(3)+1:F+N_D(3)+D_r(3)-F)]; 2200+(N_D(3)+1:N_D(3)+D_r(3))], [], 1);
% % tx_bit_enc{4} = reshape([[(-1)*ones(1, F) 3000+(F+N_D(4)+1:F+N_D(4)+D_r(4)-F)]; [(-1)*ones(1, F) 3100+(F+N_D(4)+1:F+N_D(4)+D_r(4)-F)]; 3200+(N_D(4)+1:N_D(4)+D_r(4))], [], 1);



% A1 = randi(2, D_r(1)-F, 1)-1;
% B1 = randi(2, D_r(1)-F, 1)-1;
% C1 = randi(2, D_r(1), 1)-1;
% A2 = randi(2, D_r(1)-F, 1)-1;
% B2 = randi(2, D_r(1)-F, 1)-1;
% C2 = randi(2, D_r(1), 1)-1;
% 
% tx_bit_enc_ref{1} = [(-1)*ones(F, 1); A1; (-1)*ones(F, 1); B1; C1];
% tx_bit_enc_ref{2} = [(-1)*ones(F, 1); A2; (-1)*ones(F, 1); B2; C2];
% % tx_bit_enc_ref{3} = [(-1)*ones(F, 1); 2000+(F+N_D(3)+1:F+N_D(3)+D_r(3)-F)'; (-1)*ones(F, 1); 2100+(F+N_D(3)+1:F+N_D(3)+D_r(3)-F)'; 2200+(N_D(3)+1:N_D(3)+D_r(3))'];
% % tx_bit_enc_ref{4} = [(-1)*ones(F, 1); 3000+(F+N_D(4)+1:F+N_D(4)+D_r(4)-F)'; (-1)*ones(F, 1); 3100+(F+N_D(4)+1:F+N_D(4)+D_r(4)-F)'; 3200+(N_D(4)+1:N_D(4)+D_r(4))'];
% 
% tx_bit_enc{1} = reshape([[(-1)*ones(1, F) A1']; [(-1)*ones(1, F) B1']; C1'], [], 1);
% tx_bit_enc{2} = reshape([[(-1)*ones(1, F) A2']; [(-1)*ones(1, F) B2']; C2'], [], 1);
% % tx_bit_enc{3} = reshape([[(-1)*ones(1, F) 2000+(F+N_D(3)+1:F+N_D(3)+D_r(3)-F)]; [(-1)*ones(1, F) 2100+(F+N_D(3)+1:F+N_D(3)+D_r(3)-F)]; 2200+(N_D(3)+1:N_D(3)+D_r(3))], [], 1);
% % tx_bit_enc{4} = reshape([[(-1)*ones(1, F) 3000+(F+N_D(4)+1:F+N_D(4)+D_r(4)-F)]; [(-1)*ones(1, F) 3100+(F+N_D(4)+1:F+N_D(4)+D_r(4)-F)]; 3200+(N_D(4)+1:N_D(4)+D_r(4))], [], 1);





% rate matching
rv_idx = 0;                                                 % redundancy version number
k_0 = rm.R_TC.*(2*ceil(rm.N_cb./(8*rm.R_TC))*rv_idx+2);     % start index of circular buffer

tx_bit_ratematch_ref = lteRateMatchTurbo(tx_bit_enc_ref, rm.G, rv_idx);

e_k = cell(1, cc.C);
e_k_ = cell(1, cc.C);
f_k = zeros(sum(rm.E), 1);
for r = 1:cc.C
    % rate matching input
    d_k = reshape(tx_bit_enc{r}, 3, [])';
    if r == 1
        d_k(1:cc.F, 1:2) = -1;      % null for filler bits
    end
    
    % sub-block interleaver
    y_k = [(-1)*ones(rm.N_D(r), 3); d_k];
    v_k = [y_k(rm.Y_P{r}, 1:2) y_k(rm.PI_k{r}, 3)];
    w_k = [v_k(:, 1); reshape(v_k(:, 2:3)', [], 1)];
    
    % bit collection, selection and transmission
    w_k_ = w_k;
    w_k_idx = 0:rm.K_w(r)-1;
    w_k_(w_k == -1) = [];
    w_k_idx(w_k == -1) = [];
    k_0_ = find(w_k_idx >= k_0(r), 1)-1;
    N_cb_ = length(w_k_);
%     e_k{r} = w_k_(mod(rm.k_0(r):rm.k_0(r)+rm.E(r)-1, N_cb_)+1);
    e_k{r} = w_k_(mod(k_0_:k_0_+rm.E(r)-1, N_cb_)+1);
    
    e_k_{r} = zeros(rm.E(r), 1);
    k = 0;
    j = 0;
    while k < rm.E(r)
        w_k_tmp = w_k(mod(k_0(r)+j, rm.N_cb(r))+1);
        if w_k_tmp ~= -1
            e_k_{r}(k+1) = w_k_tmp;
            k = k+1;
        end
        j = j+1;
    end
    
    
    
    
    
    % code block concatenation
    f_k(sum(rm.E(1:r-1))+1:sum(rm.E(1:r))) = e_k{r};
    f_k_(sum(rm.E(1:r-1))+1:sum(rm.E(1:r))) = e_k_{r};
end

% tx_bit_ratematch = [0; 0; f_k(1:end-2)];
tx_bit_ratematch = f_k;
tx_bit_ratematch_ = f_k_;




for r = 1:C
    a = double(reshape(tx_bit_enc_ref{r}, [], 3));
    b = reshape(tx_bit_enc{r}, 3, [])';
    if r == 1
        b(1:F, 1:2) = -1;
    end
    figure
    subplot(3, 1, 1), plot(xor(a(:, 1), b(:, 1)), '-b.'), grid minor
    subplot(3, 1, 2), plot(xor(a(:, 2), b(:, 2)), '-b.'), grid minor
    subplot(3, 1, 3), plot(xor(a(:, 3), b(:, 3)), '-b.'), grid minor
    pause
end

% figure, plot(tx_bit_ratematch_ref(:), '-bo'), hold on, plot(tx_bit_ratematch(:), ':rx'), hold off, grid minor
% figure, plot((tx_bit_ratematch_ref(:)-tx_bit_ratematch(:)), '-b.'), grid minor
% figure, plot(xor(tx_bit_ratematch_ref(:), tx_bit_ratematch(:)), '-b.'), grid minor
figure, plot(tx_bit_ratematch_ref(:), '-bo'), hold on, plot(tx_bit_ratematch_(:), ':rx'), hold off, grid minor
figure, plot((double(tx_bit_ratematch_ref(:))-tx_bit_ratematch_(:)), '-b.'), grid minor
figure, plot(xor(tx_bit_ratematch_ref(:), tx_bit_ratematch_(:)), '-b.'), grid minor











