N_RE = 12*14-24;
R = 3/4;
Q_m = 4;
v = 1;
N_SC = 12;
N_symbol = 14;
N_DMRS = 24;
N_oh = 0;

N_info = N_RE*R*Q_m*v;
N_RE_ = N_SC*N_symbol-N_DMRS-N_oh;

if N_info <= 3824
    n = max(3, floor(log2(N_info))-6);
    N_info_ = max(24, (2^n)*floor(N_info/(2^n)));
else
    n = floor(log2(N_info))-5;
    N_info_ = (2^n)*round((N_info-24)/(2^n));
    
    if R <= 1/4
        C = ceil((N_info_+24)/3816);
        TBS = 8*C*ceil((N_info_+24)/(8*C))-24;
    else
        if N_info_ >= 8424
            C = ceil((N_info_+24)/8424);
            TBS = 8*C*ceil((N_info_+24)/(8*C))-24;
        else
            TBS = 8*ceil((N_info_+24)/8)-24;
        end
    end
end


