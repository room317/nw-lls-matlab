% gen_real_ch reproduce real time-frequency channel matrix from channel (time saving)
% object.
%     [ch_onetap_tf, ch_mat_tf, ch_eff_tf, ch_eff_dd] = gen_real_ch(ch_obj, ch_pg)
%     ch_onetap_tf: time-frequency channel matrix (diagonal elements only)
%     ch_mat_tf: time-frequency channel matrix (full matrix)
%     ch_eff_tf: time-frequency effective channel matrix
%     ch_eff_dd: delay-doppler effective channel matrix
%     ch_obj: channel object (eg. comm.RayleighChannel)
%     ch_pg: channel path gain
%     nfft: fft size
%     ncp: number of cp samples
%     nbw: number of inband subcarriers
%     nsym: number of ofdm(otfs) symbols
%     list_subc_usr: index vector of user subcarriers
%     list_sym_usr: index vector of user ofdm(otfs) symbols
%     demap_usr: user demapping option (true: generate all output, false: generate 'ch_mat_tf' only)

function [ch_mat_tf, ch_onetap_usr_tf, ch_eff_usr_tf, ch_eff_usr_dd] = gen_real_ch_r1(ch_obj, ch_pg, nfft, ncp, nbw, nsym, list_subc_usr, list_sym_usr, demap_usr)

% get channel object info
ch_info = info(ch_obj);
ch_filter = ch_info.ChannelFilterCoefficients;    % ch_filter_coeff: channel filter coefficient per path
npath = length(ch_obj.PathDelays);                % num_path: number of path

% initialize
nfilter = size(ch_filter, 2);
ch_coeff_mat = zeros(nfft, nfft, npath);
ch_mat_t = zeros(nfft, nfft, nsym);
ch_halfmap_mat_tf = zeros(nfft, nfft, nsym);

% reshape path gain matrix
ch_pg_reshape = reshape(ch_pg, nfft+ncp, nsym, npath);
ch_pg_reshape(1:ncp, :, :) = [];

% generate channel coeff. matrix
for idx_path = 1:npath
    ch_coeff_row = circshift([ch_filter(idx_path, end:-1:1) zeros(1, nfft-nfilter)], -nfilter+1);
    ch_coeff_col = [ch_filter(idx_path, :) zeros(1, nfft-nfilter)];
    ch_coeff_mat(:, :, idx_path) = toeplitz(ch_coeff_col, ch_coeff_row);
end

% generate tf domain channel (circulant matrix svd)
for idx_sym = 1:nsym
    pg_mat = ch_pg_reshape(:, idx_sym, :);
    ch_mat_per_path = pg_mat.*ch_coeff_mat;
    ch_mat_t(:, :, idx_sym) = sum(ch_mat_per_path, 3);
    ch_halfmap_mat_tf(:, :, idx_sym) = ifft(fft(ch_mat_t(:, :, idx_sym), [], 1), [], 2);    % circulant matrix svd
end

% demap bandwidth
ch_mat_shift_tf = fftshift(fftshift(ch_halfmap_mat_tf, 2), 1);
ch_mat_tf = ch_mat_shift_tf((nfft/2)-(nbw/2)+1:(nfft/2)+(nbw/2), (nfft/2)-(nbw/2)+1:(nfft/2)+(nbw/2), :);

% generate output
if demap_usr
    % demap user block
    nsubc_usr = length(list_subc_usr);
    nsym_usr = length(list_sym_usr);
    ch_mat_usr_tf = ch_mat_tf(list_subc_usr, list_subc_usr, list_sym_usr);
    
    % extract diagonal elements of real channel
    ch_onetap_usr_tf = zeros(nsubc_usr, nsym_usr);
    for idx_sym = list_sym_usr
        ch_onetap_usr_tf(:, idx_sym) = diag(ch_mat_usr_tf(:, :, idx_sym));  % diagonal term only
    end
    
    % generate effective time-frequency channel matrix
    ch_eff_usr_tf = zeros(nsubc_usr*nsym_usr);
    for idx_sym = list_sym_usr
        ch_eff_usr_tf((idx_sym-1)*nsubc_usr+1:idx_sym*nsubc_usr, (idx_sym-1)*nsubc_usr+1:idx_sym*nsubc_usr) = ch_mat_tf(:, :, idx_sym);
    end
    
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
    ch_onetap_usr_tf = [];
    ch_eff_usr_tf = [];
    ch_eff_usr_dd = [];
end

% % dump
% assignin('base', 'ch_mat_tf', ch_mat_tf);
% assignin('base', 'ch_onetap_tf', ch_onetap_tf);
% assignin('base', 'idft_column', idft_column);
% assignin('base', 'dft_row', dft_row);
% assignin('base', 'sfft_mtx', sfft_mtx);
% assignin('base', 'ch_eff_tf', ch_eff_tf);
% assignin('base', 'ch_eff_dd', ch_eff_dd);
% pause

end