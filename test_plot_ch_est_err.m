% plot one-tap channel estimation error

% load sum_ch_pwr first!
if ~exist('ch_real_eff_usr_dd', 'var') || ...
        ~exist('ch_real_eff_usr_tf', 'var')
    error('Load variables first!')
end

% set parameters (manually)
df = 60e3;
nfft = 256;
nsubc = 11*12;
nsym = 4*14;
if ~all(size(ch_real_eff_usr_dd) == size(ch_real_eff_usr_tf)) || ...
        ~(size(ch_real_eff_usr_dd, 1) == size(ch_real_eff_usr_dd, 2)) || ...
        ~(nsubc*nsym == size(ch_real_eff_usr_dd, 1))
    error('Bad parameters!')
end

% set parameters (automatically)
fs = df*nfft;
Ts = 1/fs;
ncp = nfft*72/1024;
Tsym = Ts*(nfft+ncp);
fbw = df*nsubc;
T = Tsym*nsym;
ddelay = 1/fbw;
ddoppler = 1/T;
rdelay = ddelay*nsubc;
rdoppler = ddoppler*nsym;

% sfft transform matrix
sfft_mtx = kron(dftmtx(nsym), conj(dftmtx(nsubc))/nsubc);
isfft_mtx = kron(conj(dftmtx(nsym))/nsym, dftmtx(nsubc));

% set index (scaled)
x_axis = ((0:nsym-1)-(nsym/2))*ddoppler;
y_axis = ((0:nsubc-1)-(nsubc/2))*ddelay;

% % set index (not scaled)
% x_axis = (0:nsym-1)-(nsym/2);
% y_axis = (0:nsubc-1)-(nsubc/2);

% estimate 1-tap channel
ch_est_onetap_tf = diag(diag(ch_real_eff_usr_tf));
ch_est_onetap_dd = gen_eff_ch(sqrt(nsym/nsubc)*fft(ifft(reshape(diag(ch_real_eff_usr_tf), nsubc, nsym), [], 2), [], 1));
ch_intf_eff_usr_tf = ch_real_eff_usr_tf - ch_est_onetap_tf;
ch_intf_eff_usr_dd = sfft_mtx*ch_intf_eff_usr_tf*isfft_mtx;
tmp_ch = ch_real_eff_usr_dd - ch_est_onetap_dd;

% % test
% tx_dd = complex(randn(nsubc*nsym, 1), randn(nsubc*nsym, 1));
% rx_dd = ch_real_eff_usr_dd*tx_dd;
% rx_onetap_dd = ch_est_onetap_dd*tx_dd;
% tx_tf = reshape(sqrt(nsym/nsubc)*fft(ifft(reshape(tx_dd, nsubc, nsym), [], 2), [], 1), [], 1);
% rx_tf = ch_real_eff_usr_tf*tx_tf;
% rx_onetap_tf = ch_est_onetap_tf*tx_tf;

% plot
ntest = 1000;

figure
mesh(1:ntest, 1:ntest, abs(ch_real_eff_usr_tf(1:ntest, 1:ntest)))
zlabel('Channel Amplitude'), title('Real Channel in TF'), axis tight, grid minor

figure
mesh(1:ntest, 1:ntest, abs(ch_est_onetap_tf(1:ntest, 1:ntest)))
zlabel('Channel Amplitude'), title('One-tap Channel Estimation in TF'), axis tight, grid minor

figure
mesh(1:ntest, 1:ntest, abs(ch_real_eff_usr_tf(1:ntest, 1:ntest)-ch_est_onetap_tf(1:ntest, 1:ntest)))
zlabel('Channel Amplitude'), title('One-tap Channel Estimation in TF'), axis tight, grid minor

figure
mesh(1:ntest, 1:ntest, abs(ch_real_eff_usr_dd(1:ntest, 1:ntest)))
zlabel('Channel Amplitude'), title('Real Channel in DD'), axis tight, grid minor

figure
mesh(1:ntest, 1:ntest, abs(ch_est_onetap_dd(1:ntest, 1:ntest)))
zlabel('Channel Amplitude'), title('One-tap Channel Estimation in DD'), axis tight, grid minor

figure
mesh(1:ntest, 1:ntest, abs(ch_real_eff_usr_dd(1:ntest, 1:ntest)-ch_est_onetap_dd(1:ntest, 1:ntest)))
zlabel('Channel Amplitude'), title('One-tap Channel Estimation in DD'), axis tight, grid minor

% use 'xlim' and 'ylim' to save axes
