% rate match for turbo code (cell-based signal processing)
%   - test_option.rm_version: 1(do-while operation), 2(vector operation)
%   - d must shall be cell-type
%   - rv_idx: redundancy version number

function f_k = tx_ratematch_r3(d, cc, rm, test_option)

% start index of circular buffer
k_0 = rm.R_TC.*(2*ceil(rm.N_cb./(8*rm.R_TC))*test_option.rv_idx+2);

e_k = cell(1, cc.C);
f_k = zeros(sum(rm.E), 1);
for r = 1:cc.C
    % rate matching input
    d_k = reshape(d{r}, 3, [])';
    if r == 1
        d_k(1:cc.F, 1:2) = -1;          % null for filler bits
    end
    
    % sub-block interleaver
    y_k = [(-1)*ones(rm.N_D(r), 3); d_k];
    v_k = [y_k(rm.Y_P{r}, 1:2) y_k(rm.PI_k{r}, 3)];
    w_k = [v_k(:, 1); reshape(v_k(:, 2:3)', [], 1)];
    
    % bit collection, selection and transmission
    if test_option.rm_version == 1      % follows standard doc
        e_k{r} = zeros(rm.E(r), 1);
        k = 0;
        j = 0;
        while k < rm.E(r)
            w_k_tmp = w_k(mod(k_0(r)+j, rm.N_cb(r))+1);
            if w_k_tmp ~= -1
                e_k{r}(k+1) = w_k_tmp;
                k = k+1;
            end
            j = j+1;
        end
    else                                % vectorized equations
        w_k_idx = 0:rm.K_w(r)-1;
        w_k_idx(w_k == -1) = [];
        w_k(w_k == -1) = [];
        k_0_ = find(w_k_idx >= k_0(r), 1)-1;
        N_cb_ = length(w_k);
        e_k{r} = w_k(mod(k_0_:k_0_+rm.E(r)-1, N_cb_)+1);
    end
    
    % code block concatenation
    f_k(sum(rm.E(1:r-1))+1:sum(rm.E(1:r))) = e_k{r};
end

% dump
assignin('base', 'k_0', k_0)
assignin('base', 'd_k', d_k)
assignin('base', 'y_k', y_k)
assignin('base', 'v_k', v_k)
assignin('base', 'w_k', w_k)
assignin('base', 'e_k', e_k)
assignin('base', 'f_k', f_k)

end
