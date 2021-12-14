% ref: 3gpp ts 36.212 (36212-e40)
% cell based channel coding parameters

function nw_cc = nw_cc_prm_r1(len_tb_bit)

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
nw_cc.tc_prm = tc_prm;              % turbo code table
nw_cc.tc_trellis = tc_trellis;      % turbo code trellis structure
nw_cc.gCRC24A = gCRC24A;            % cyclic generator polynomial (gCRC24A)
nw_cc.gCRC24B = gCRC24B;            % cyclic generator polynomial (gCRC24B)
nw_cc.num_iter_max = num_iter_max;  % number of maximum iteration
nw_cc.A = A;                        % size of the input sequence
nw_cc.L = L;                        % number of parity bits
nw_cc.B = B;                        % number of bits after crc attachment
nw_cc.Z = Z;                        % maximum code block size (Z = 6144)
nw_cc.Z_min = Z_min;                % minimum code block size (40)
nw_cc.L_ = L_;                      % additional crc sequence length attached to each code block (L_ = 24 bits)
nw_cc.C = C;                        % total number of code blocks
nw_cc.B_ = B_;                      % number of bits after additional crc attachment to each code block
nw_cc.K_plus = K_plus;              % number of bits per code block (first segmentation size)
nw_cc.f1_plus = f1_plus;            % turbo code internal interleaver parameter corresponds to K_plus
nw_cc.f2_plus = f2_plus;            % turbo code internal interleaver parameter corresponds to K_plus
nw_cc.C_plus = C_plus;              % number of segments of size K_plus
nw_cc.K_minus = K_minus;            % number of bits per code block (second segmentation size)
nw_cc.f1_minus = f1_minus;          % turbo code internal interleaver parameter corresponds to K_minus
nw_cc.f2_minus = f2_minus;          % turbo code internal interleaver parameter corresponds to K_minus
nw_cc.C_minus = C_minus;            % number of segments of size K_minus
nw_cc.F = F;                        % number of filler bits
nw_cc.K_r = K_r;                    % number of bits for the code block number r
nw_cc.f1_r = f1_r;                  % turbo code internal interleaver parameter for the code block number r
nw_cc.f2_r = f2_r;                  % turbo code internal interleaver parameter for the code block number r
nw_cc.D_r = D_r;                    % number of encoded bits per output stream
nw_cc.PI_r = PI_r;                  % turbo code internal interleaver

end
