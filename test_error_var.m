% purpose
% 0. to test the dd-domain input-output relation
% 1. to check per-qam noise variance before/after equalization (dd domain)
% 2. to decide which noise variance should be used for llr calculation
% 
% how to use
% 0. set channel estimation option to 'perfect channel'
% 1. dump these variables at the otfs top:
%   - 'tx_sym_rbs_dd'
%   - 'rx_sym_rbs_usr_tf'
%   - 'ch_real_eff_usr_dd'
%   - 'ch_est_rbs_usr_dd'
%   - 'noise_var'
%   - 'rx_sym_rbs_usr_eq_tf'
%   - 'rx_sym_rbs_usr_eq_dd'
% 2. run

tic

% check variables
if ~exist('tx_sym_rbs_dd', 'var') || ...
        ~exist('rx_sym_rbs_usr_tf', 'var') || ...
        ~exist('ch_real_eff_usr_tf', 'var') || isempty(ch_real_eff_usr_tf) || ...
        ~exist('ch_real_eff_usr_dd', 'var') || isempty(ch_real_eff_usr_dd) || ...
        ~exist('ch_est_rbs_usr_dd', 'var') || isempty(ch_est_rbs_usr_dd) || ...
        ~exist('noise_var', 'var') || ...
        ~exist('rx_sym_rbs_usr_eq_tf', 'var') || ...
        ~exist('rx_sym_rbs_usr_eq_dd', 'var')
    error('Dump variables first!')
end

toc

% set signals
x_d = reshape(tx_sym_rbs_dd(:, :, 1, 1), [], 1);        % tx in dd-eff domain
r_d = reshape(sqrt(nw_num.num_subc_usr/nw_num.num_ofdmsym_usr)*fft(ifft(rx_sym_rbs_usr_tf, [], 1), [], 2), [], 1);          % rx in dd domain
h_d = ch_real_eff_usr_dd;                               % ch_real in dd-eff domain
h_t = ch_real_eff_usr_tf;                               % ch_real in tf-eff domain
hh0_dd = ch_est_rbs_usr_dd;                             % ch_est in dd domain
x_eq_d = rx_sym_rbs_usr_eq_dd(:);                       % tx eq in dd-eff domain
x_eq_t = rx_sym_rbs_usr_eq_tf(:);                       % tx eq in tf-eff domain

toc

% find basic signals
n_d = r_d-h_d*x_d;                                      % noise in dd-eff domain
h_dd = gen_conv_ch(h_d, nw_num.num_subc_usr);           % full real matrix in dd domain

toc

% channel estimation 0
% % % % hh0_d = gen_eff_ch(hh0_dd);                             % ch_est in dd-eff domain
hh0_tf = sqrt(nw_num.num_doppler_usr/nw_num.num_delay_usr)*fft(ifft(hh0_dd, [], 2), [], 1);     % ch_est in tf domain
hh0_t = diag(hh0_tf(:));                                   % ch_est in tf-eff domain
ihh0_t = diag(hh0_tf(:).^-1);                              % inv(ch_est) in tf-eff domain
ihh0_d = nw_num.sfft_mtx*ihh0_t*nw_num.isfft_mtx;       % inv(ch_est) in dd-eff domain
hh0_rmse = sqrt(mean(abs(h_t-hh0_t).^2, 'all'));        % channel estimation error (rmse)
fprintf('ch. est. error(rmse) with pilot in dd        :%10.6f\n', hh0_rmse)

toc

% channel estimation 1 (perfect in tf domain: diagonal in tf domain)
hh1_tf = reshape(diag(h_t), nw_num.num_subc_usr, nw_num.num_ofdmsym_usr);                       % ch_est in tf domain
% % % % hh1_dd = sqrt(nw_num.num_subc_usr/nw_num.num_ofdmsym_usr)*fft(ifft(hh1_tf, [], 1), [], 2);      % ch_est in dd domain
hh1_t = diag(hh1_tf(:));                                   % ch_est in tf-eff domain
ihh1_t = diag(hh1_tf(:).^-1);                           % inv(ch_est) in tf-eff domain
ihh1_d = nw_num.sfft_mtx*ihh1_t*nw_num.isfft_mtx;       % inv(ch_est) in dd-eff domain
hh1_rmse = sqrt(mean(abs(h_t-hh1_t).^2, 'all'));        % channel estimation error (rmse)
fprintf('ch. est. error(rmse) with diag element in tf :%10.6f\n', hh1_rmse)

toc

% channel estimation 2 (perfect in dd domain: center pilot with full resource)
hh2_dd = reshape(h_dd(:, :, nw_num.num_subc_usr*nw_num.num_ofdmsym_usr/2), nw_num.num_subc_usr, nw_num.num_ofdmsym_usr);    % ch_est in dd domain
hh2_tf = sqrt(nw_num.num_doppler_usr/nw_num.num_delay_usr)*fft(ifft(hh2_dd, [], 2), [], 1);                                 % ch_est in tf domain
hh2_t = diag(hh2_tf(:));                                   % ch_est in tf-eff domain
ihh2_t = diag(hh2_tf(:).^-1);                           % inv(ch_est) in tf-eff domain
ihh2_d = nw_num.sfft_mtx*ihh2_t*nw_num.isfft_mtx;       % inv(ch_est) in dd-eff domain
hh2_rmse = sqrt(mean(abs(h_t-hh2_t).^2, 'all'));        % channel estimation error (rmse)
fprintf('ch. est. error(rmse) with perfect pilot in dd:%10.6f\n', hh2_rmse)

toc

% equalization 1 (with channel estimation 1)
y0_d = ihh0_d*r_d;                                      % eq in dd-eff domain
% % % % y0_d_tmp = ihh0_d*h_d*x_d+ihh0_d*n_d;                   % eq in dd-eff domain for comparison
e0_d = y0_d-x_d;                                            % est error in dd domain
% % % % e0_d_tmp = (ihh0_d*h_d-eye(length(x_d)))*x_d+ihh0_d*n_d;    % est error in dd domain for comparison
e0_rmse = sqrt(mean(abs(x_d-y0_d).^2, 'all'));          % equalized symbol error (rmse)
fprintf('sym. error(rmse) with pilot in dd        :%10.6f\n', e0_rmse)

toc

% equalization 2 (with channel estimation 2)
y1_d = ihh1_d*r_d;                                      % eq in dd-eff domain
% % % % y1_d_tmp = ihh1_d*h_d*x_d+ihh1_d*n_d;                   % eq in dd-eff domain for comparison
e1_d = y1_d-x_d;                                            % est error in dd domain
% % % % e1_d_tmp = (ihh1_d*h_d-eye(length(x_d)))*x_d+ihh1_d*n_d;    % est error in dd domain for comparison
e1_rmse = sqrt(mean(abs(x_d-y1_d).^2, 'all'));          % equalized symbol error (rmse)
fprintf('sym. error(rmse) with diag element in tf :%10.6f\n', e1_rmse)

toc

% equalization 3 (with channel estimation 3)
y2_d = ihh2_d*r_d;                                      % eq in dd-eff domain
% % % % y2_d_tmp = ihh2_d*h_d*x_d+ihh2_d*n_d;                   % eq in dd-eff domain for comparison
e2_d = y2_d-x_d;                                            % est error in dd domain
% % % % e2_d_tmp = (ihh2_d*h_d-eye(length(x_d)))*x_d+ihh2_d*n_d;    % est error in dd domain for comparison
e2_rmse = sqrt(mean(abs(x_d-y2_d).^2, 'all'));          % equalized symbol error (rmse)
fprintf('sym. error(rmse) with perfect pilot in dd:%10.6f\n', e2_rmse)

toc

% plot
figure, plot(abs(n_d).^2, '-b.'), hold on, plot(noise_var*ones(size(n_d)), '-r.'), hold off, legend('real noise', 'given noise variance'), grid minor
figure, plot(real(r_d), '-b.'), hold on, plot(real(h_d*x_d), ':r.'), hold off, legend('rx in dd', 'rx in dd without noise'), grid minor
% % % % figure, plot(real(y0_d), '-b.'), hold on, plot(real(y0_d_tmp), ':r.'), hold off, legend('eq0 rx in dd', 'for comparison'), grid minor
% % % % figure, plot(real(y1_d), '-b.'), hold on, plot(real(y1_d_tmp), ':r.'), hold off, legend('eq1 rx in dd', 'for comparison'), grid minor
% % % % figure, plot(real(y2_d), '-b.'), hold on, plot(real(y2_d_tmp), ':r.'), hold off, legend('eq2 rx in dd', 'for comparison'), grid minor
% % % % figure, plot(real(e0_d), '-b.'), hold on, plot(real(e0_d_tmp), ':r.'), hold off, legend('sym0 err in dd', 'for comparison'), grid minor
% % % % figure, plot(real(e1_d), '-b.'), hold on, plot(real(e1_d_tmp), ':r.'), hold off, legend('sym1 err in dd', 'for comparison'), grid minor
% % % % figure, plot(real(e2_d), '-b.'), hold on, plot(real(e2_d_tmp), ':r.'), hold off, legend('sym2 err in dd', 'for comparison'), grid minor
figure, plot(real(x_eq_d), '-b.'), hold on, plot(real(y0_d), ':r.'), hold off, legend('eq dump in dd', 'eq0 rx cal'), grid minor
figure, plot(abs(e0_d).^2, '-b'), hold on, plot(abs(e1_d).^2, '-r'), plot(abs(e2_d).^2, '-k'), hold off, legend('sym0 err pwr', 'sym1 err pwr', 'sym2 err pwr'), grid minor

% check fast fading in dd domain
% % % % for i = 1:10:nw_num.num_subc_usr*nw_num.num_ofdmsym_usr
% % % %     figure(21), mesh(1:size(h_dd, 2), 1:size(h_dd, 1), fftshift(fftshift(abs(h_dd(:, :, i)), 1), 2))
% % % % end

toc
