% gen_real_ch reproduce real time-frequency channel matrix from channel (time saving)
% object.
%     [ch_onetap_tf, ch_mat_tf, ch_eff_tf, ch_eff_dd] = gen_real_ch(ch_obj, ch_pg)
%     ch_onetap_tf: time-frequency channel matrix (diagonal elements only)
%     ch_fulltap_tf: time-frequency channel matrix (full matrix)
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
%     gpu_flag: use gpu when true

function [ch_mat_t, ch_fulltap_tf, ch_onetap_usr_tf, ch_eff_usr_tf, ch_eff_usr_dd] = gen_real_ch_r1(ch_obj, ch_pg, num, list_subc_usr, list_ofdmsym_usr, demap_usr, test_option)

% variables  
nfft = num.num_fft;
ncp = num.num_cp;
nbw = num.num_subc_bw;
nsym = num.num_ofdmsym;
sfft_mtx = num.sfft_mtx;
isfft_mtx = num.isfft_mtx;

% get channel object info
% ch_info = info(ch_obj);
% ch_filter = ch_info.ChannelFilterCoefficients;    % ch_filter_coeff: channel filter coefficient per path
% npath = length(ch_obj.PathDelays);                % num_path: number of path
ch_filter = getPathFilters(ch_obj);    % ch_filter_coeff: channel filter coefficient per path
ch_filter = ch_filter.';
npath = size(ch_pg, 2);

% initialize
nfilter = size(ch_filter, 2);
if test_option.gpu_flag
    ch_coeff_mat = gpuArray(zeros(nfft, nfft, npath));
    ch_mat_t = gpuArray(zeros(nfft, nfft, nsym));
    ch_mat_fftshift_tf = gpuArray(zeros(nfft, nfft, nsym));
else
    ch_coeff_mat = zeros(nfft, nfft, npath);
    ch_mat_t = zeros(nfft, nfft, nsym);
    ch_mat_fftshift_tf = zeros(nfft, nfft, nsym);
end

% reshape path gain matrix
if test_option.gpu_flag
    ch_pg_reshape = gpuArray(reshape(ch_pg, nfft+ncp, nsym, npath));
else
    ch_pg_reshape = reshape(ch_pg, nfft+ncp, nsym, npath);
end
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
    ch_mat_fftshift_tf(:, :, idx_sym) = ifft(fft(ch_mat_t(:, :, idx_sym), [], 1), [], 2);    % circulant matrix svd
end

% convert gpu arrays to normal arrays
if test_option.gpu_flag
    ch_mat_fftshift_tf = gather(ch_mat_fftshift_tf);
    clear ch_mat_t
    clear ch_mat_per_path
    clear pg_mat
    clear ch_coeff_mat
    clear ch_pg_reshape
end

% demap bandwidth
ch_mat_tf = fftshift(fftshift(ch_mat_fftshift_tf, 2), 1);
ch_fulltap_tf = ch_mat_tf((nfft/2)-(nbw/2)+1:(nfft/2)+(nbw/2), (nfft/2)-(nbw/2)+1:(nfft/2)+(nbw/2), :);

% generate output
if demap_usr
    % demap user block
    nsubc_usr = length(list_subc_usr);
    nsym_usr = length(list_ofdmsym_usr);
    ch_mat_usr_tf = ch_fulltap_tf(list_subc_usr, list_subc_usr, list_ofdmsym_usr);
    
    % extract diagonal elements of real channel
    ch_onetap_usr_tf = zeros(nsubc_usr, nsym_usr);
    for idx_sym = list_ofdmsym_usr
        ch_onetap_usr_tf(:, idx_sym) = diag(ch_mat_usr_tf(:, :, idx_sym));  % diagonal term only
    end
    
    % generate effective time-frequency channel matrix
    ch_eff_usr_tf = zeros(nsubc_usr*nsym_usr);
    for idx_sym = 1:nsym_usr
        ch_eff_usr_tf((idx_sym-1)*nsubc_usr+1:idx_sym*nsubc_usr, (idx_sym-1)*nsubc_usr+1:idx_sym*nsubc_usr) = ch_mat_usr_tf(:, :, idx_sym);
    end
    
    % generate effective delay-doppler channel matrix
    % ch_eff_dd = sfft_mtx*ch_eff_tf/sfft_mtx;
    ch_eff_usr_dd = sfft_mtx*ch_eff_usr_tf*isfft_mtx;
else
    ch_onetap_usr_tf = [];
    ch_eff_usr_tf = [];
    ch_eff_usr_dd = [];
end

% % dump
% assignin('base', 'ch_fulltap_tf', ch_fulltap_tf);
% assignin('base', 'ch_onetap_tf', ch_onetap_tf);
% assignin('base', 'idft_column', idft_column);
% assignin('base', 'dft_row', dft_row);
% assignin('base', 'ch_eff_tf', ch_eff_tf);
% assignin('base', 'ch_eff_dd', ch_eff_dd);
% pause

end
