% demap_real_ch demaps user channel from real time-frequency channel
% matrix.
%     demap_real_ch(ch_mat_tf, list_subc_usr, list_sym_usr, cheq_option, test_option)
%       - ch_onetap_usr_tf: time-frequency channel matrix (diagonal elements only)
%       - ch_eff_usr_tf: time-frequency effective user channel matrix
%       - ch_eff_usr_dd: delay-doppler effective user channel matrix
%       - ch_mat_tf: time-frequency channel matrix (full matrix)
%       - list_subc_usr: index vector of user subcarriers
%       - list_ofdmsym_usr: index vector of user ofdm(otfs) symbols
%       - chest_option: channel estimation option
%       - cheq_option: channel equalization option
%       - test_option: full-tap/one-tap equalization option

function [ch_onetap_usr_tf, ch_eff_usr_tf, ch_eff_usr_dd] = demap_real_ch(ch_mat_tf, list_subc_usr, list_ofdmsym_usr, chest_option, cheq_option, test_option)

% demap user channel matrix
nsubc_usr = length(list_subc_usr);
nsym_usr = length(list_ofdmsym_usr);
ch_mat_usr_tf = ch_mat_tf(list_subc_usr, list_subc_usr, list_ofdmsym_usr);

% generate real full-tap time-frequency channel matrix
if (strcmp(chest_option, 'real') && test_option.fulltap_eq) || test_option.ch_mse
    % generate effective time-frequency channel matrix
    ch_eff_usr_tf = zeros(nsubc_usr*nsym_usr);
    for idx_sym = 1:nsym_usr
        ch_eff_usr_tf((idx_sym-1)*nsubc_usr+1:idx_sym*nsubc_usr, (idx_sym-1)*nsubc_usr+1:idx_sym*nsubc_usr) = ch_mat_usr_tf(:, :, idx_sym);
    end
else
    ch_eff_usr_tf = [];
end

% generate real full-tap delay-doppler channel matrix
if strcmp(chest_option, 'real') && test_option.fulltap_eq && (strcmp(cheq_option, 'ddeq_zf') || strcmp(cheq_option, 'ddeq_mmse'))
    % generate sfft matrix
    % idft_column = kron(eye(nsym), conj(dftmtx(nbw))/sqrt(nbw));
    % dft_row = kron(dftmtx(nsym)/sqrt(nsym), eye(nbw));
    % sfft_mtx = dft_row*idft_column;
    sfft_mtx = kron(dftmtx(nsym_usr), conj(dftmtx(nsubc_usr))/nsubc_usr);       % kron(A, B)*kron(C, D) = kron(AC, BD)
    isfft_mtx = kron(conj(dftmtx(nsym_usr))/nsym_usr, dftmtx(nsubc_usr));     % inv(kron(A, B)) = kron(inv(A), inv(B))
    
    % generate effective delay-doppler channel matrix
    % ch_eff_dd = sfft_mtx*ch_eff_tf/sfft_mtx;
    ch_eff_usr_dd = sfft_mtx*ch_eff_usr_tf*isfft_mtx;
else
    ch_eff_usr_dd = [];
end

% generate real one-tap time-frequency channel
if (strcmp(chest_option, 'real') && ~test_option.fulltap_eq) || test_option.ch_mse
    % extract diagonal elements of real channel
    ch_onetap_usr_tf = zeros(nsubc_usr, nsym_usr);
    for idx_sym = 1:nsym_usr
        ch_onetap_usr_tf(:, idx_sym) = diag(ch_mat_usr_tf(:, :, idx_sym));  % diagonal term only
    end
else
    ch_onetap_usr_tf = [];
end

end
