% how to use
% 0. set channel estimation option to 'perfect channel'
% 1. dump these variables at the otfs top:
%   - 'ch_est_rbs_usr_tf'
%   - 'rx_sym_rbs_usr_tf'
%   - 'noise_var'
%   - 'tx_sym_rbs_tf'
%   - 'tx_sym_rbs_dd'
% 2. run

% memo for perfect channel
%   - ch_est_rbs_usr_tf : exist
%   - ch_est_rbs_usr_dd : empty
%   - rx_sym_rbs_usr_tf : exist
%   - rx_sym_rbs_usr_dd : empty

% check variables
if ~exist('ch_est_rbs_usr_tf', 'var') || ...
        ~exist('rx_sym_rbs_usr_tf', 'var') || ...
        ~exist('noise_var', 'var') || ...
        ~exist('tx_sym_rbs_tf', 'var') || ...
        ~exist('tx_sym_rbs_dd', 'var')
    error('Dump variables first!')
end

% sfft transform matrix
sfft_mtx = kron(dftmtx(nw_num.num_doppler_usr), conj(dftmtx(nw_num.num_delay_usr))/nw_num.num_delay_usr);       % kron(A, B)*kron(C, D) = kron(AC, BD)
isfft_mtx = kron(conj(dftmtx(nw_num.num_doppler_usr))/nw_num.num_doppler_usr, dftmtx(nw_num.num_delay_usr));     % inv(kron(A, B)) = kron(inv(A), inv(B))

% full-tap noise in dd
ch_inv_fulltap_dd = reshape(sum(inv(ch_real_eff_usr_dd), 2), nw_num.num_delay_usr, nw_num.num_doppler_usr);
noise_var_mat_fulltap_dd = noise_var*(abs(ch_inv_fulltap_dd).^2);

% full-tap noise in tf
ch_inv_fulltap_tf = reshape(sum(inv(ch_real_eff_usr_tf), 2), nw_num.num_subc_usr, nw_num.num_ofdmsym_usr);
noise_var_mat_fulltap_tf = noise_var*(abs(ch_inv_fulltap_tf).^2);

% 1. zf-eq in dd
rx_sym_rbs_usr_dd = sqrt(nw_num.num_subc_usr/nw_num.num_ofdmsym_usr)*fft(ifft(rx_sym_rbs_usr_tf, [], 1), [], 2);
ch_est_rbs_usr_dd = sqrt(nw_num.num_subc_usr/nw_num.num_ofdmsym_usr)*fft(ifft(ch_est_rbs_usr_tf, [], 1), [], 2);
ch_est_eff_usr_dd = gen_eff_ch(ch_est_rbs_usr_dd);
rx_sym_rbs_usr_eq_vec1_dd = ch_est_eff_usr_dd\rx_sym_rbs_usr_dd(:);
rx_sym_rbs_usr_eq1_dd = reshape(rx_sym_rbs_usr_eq_vec1_dd, nw_num.num_delay_usr, nw_num.num_doppler_usr);
rx_sym_rbs_usr_eq1_tf = sqrt(nw_num.num_doppler_usr/nw_num.num_delay_usr)*fft(ifft(rx_sym_rbs_usr_eq1_dd, [], 2), [], 1);

% 1.1. one-tap noise in dd
ch_inv_dd = reshape(sum(inv(ch_est_eff_usr_dd), 2), nw_num.num_delay_usr, nw_num.num_doppler_usr);
noise_var_mat1_dd = noise_var*(abs(ch_inv_dd).^2);

% 1.2. one-tap noise in tf
ch_inv_tf = 1./ch_est_rbs_usr_tf;
noise_var_mat1_tf = noise_var*(abs(ch_inv_tf).^2);

% 1.3 one-tap noise in tf alternative (not good)
ch_inv_alt_tf = sqrt(nw_num.num_doppler_usr/nw_num.num_delay_usr)*fft(ifft(ch_inv_dd, [], 2), [], 1);
noise_var_mat1_alt_tf = noise_var*(abs(ch_inv_alt_tf).^2);

% 2. mmse-eq in dd
ch_est_eff_usr_mmse_dd = (ch_est_eff_usr_dd'*ch_est_eff_usr_dd+noise_var*eye(nw_num.num_delay_usr*nw_num.num_doppler_usr))\ch_est_eff_usr_dd';
% ch_est_eff_usr_mmse_dd = (ch_est_eff_usr_dd'*ch_est_eff_usr_dd)\ch_est_eff_usr_dd';
rx_sym_rbs_usr_eq_vec2_dd = ch_est_eff_usr_mmse_dd*rx_sym_rbs_usr_dd(:);
rx_sym_rbs_usr_eq2_dd = reshape(rx_sym_rbs_usr_eq_vec2_dd, nw_num.num_delay_usr, nw_num.num_doppler_usr);
rx_sym_rbs_usr_eq2_tf = sqrt(nw_num.num_doppler_usr/nw_num.num_delay_usr)*fft(ifft(rx_sym_rbs_usr_eq2_dd, [], 2), [], 1);

% 2.1. one-tap noise in dd
ch_inv_dd = reshape(sum(ch_est_eff_usr_mmse_dd, 2), nw_num.num_delay_usr, nw_num.num_doppler_usr);
noise_var_mat2_dd = noise_var*(abs(ch_inv_dd).^2);

% 2.2. one-tap noise in tf
ch_est_eff_usr_mmse_tf = isfft_mtx*ch_est_eff_usr_mmse_dd*sfft_mtx;
ch_inv_tf = reshape(diag(ch_est_eff_usr_mmse_tf), nw_num.num_subc_usr, nw_num.num_ofdmsym_usr);
noise_var_mat2_tf = noise_var*(abs(ch_inv_tf).^2);

% 2.3 one-tap noise in tf alternative (not good)
ch_inv_alt_tf = sqrt(nw_num.num_doppler_usr/nw_num.num_delay_usr)*fft(ifft(ch_inv_dd, [], 2), [], 1);
noise_var_mat2_alt_tf = noise_var*(abs(ch_inv_alt_tf).^2);

% plot zf vs. mmse in dd (equalized in dd)
figure
subplot(3, 1, 1)
plot(real(tx_sym_rbs_dd(:)), '-k.'), hold on
plot(real(rx_sym_rbs_usr_eq1_dd(:)), '-b.')
plot(real(rx_sym_rbs_usr_eq2_dd(:)), ':r.'), hold off
title('real'), legend('tx', 'zf-eq rx', 'mmse-eq rx')
subplot(3, 1, 2)
plot(imag(tx_sym_rbs_dd(:)), '-k.'), hold on
plot(imag(rx_sym_rbs_usr_eq1_dd(:)), '-b.')
plot(imag(rx_sym_rbs_usr_eq2_dd(:)), ':r.'), hold off
title('imag'), legend('tx', 'zf-eq rx', 'mmse-eq rx')
subplot(3, 1, 3)
plot(abs(tx_sym_rbs_dd(:)-rx_sym_rbs_usr_eq1_dd(:)).^2, '-b.'), hold on
plot(abs(tx_sym_rbs_dd(:)-rx_sym_rbs_usr_eq2_dd(:)).^2, ':r.'), hold off
title('tx est. error'), legend('zf-eq', 'mmse-eq')

% plot zf vs. mmse in tf (equalized in dd)
figure
subplot(3, 1, 1)
plot(real(tx_sym_rbs_tf(:)), '-k.'), hold on
plot(real(rx_sym_rbs_usr_eq1_tf(:)), '-b.')
plot(real(rx_sym_rbs_usr_eq2_tf(:)), ':r.'), hold off
title('real'), legend('tx', 'zf-eq rx', 'mmse-eq rx')
subplot(3, 1, 2)
plot(imag(tx_sym_rbs_tf(:)), '-k.'), hold on
plot(imag(rx_sym_rbs_usr_eq1_tf(:)), '-b.')
plot(imag(rx_sym_rbs_usr_eq2_tf(:)), ':r.'), hold off
title('imag'), legend('tx', 'zf-eq rx', 'mmse-eq rx')
subplot(3, 1, 3)
plot(abs(tx_sym_rbs_tf(:)-rx_sym_rbs_usr_eq1_tf(:)).^2, '-b.'), hold on
plot(abs(tx_sym_rbs_tf(:)-rx_sym_rbs_usr_eq2_tf(:)).^2, ':r.'), hold off
title('tx est. error'), legend('zf-eq', 'mmse-eq')

% plot constellation
figure
scatter(real(tx_sym_rbs_dd(:)), imag(tx_sym_rbs_dd(:)), 'bo'), hold on
scatter(real(rx_sym_rbs_usr_eq1_dd(:)), imag(rx_sym_rbs_usr_eq1_dd(:)), 'rx'), hold off
title('equalized in dd'), legend('tx', 'rx (zf-eq)')

% plot noise variance
figure
subplot(2, 1, 1)
plot(abs(tx_sym_rbs_dd(:)-rx_sym_rbs_usr_eq1_dd(:)).^2, '-k.'), hold on
plot(noise_var_mat1_dd(:), '-b.')
plot(noise_var_mat_fulltap_dd(:), '-r.'), hold off
title('dd domain (zf-eq in dd domain)'), legend('real noise', 'est. noise (one-tap)', 'est. noise (full-tap)')
subplot(2, 1, 2)
plot(abs(tx_sym_rbs_tf(:)-rx_sym_rbs_usr_eq1_tf(:)).^2, '-k.'), hold on
plot(noise_var_mat1_tf(:), '-b.')
plot(noise_var_mat_fulltap_tf(:), '-r.')
plot(noise_var_mat1_alt_tf(:), '-m.'), hold off
title('tf domain (zf-eq in dd domain)'), legend('real noise', 'est. noise (one-tap)', 'est. noise (full-tap)', 'est. noise (alt)')

% plot noise variance
figure
subplot(2, 1, 1)
plot(abs(tx_sym_rbs_dd(:)-rx_sym_rbs_usr_eq2_dd(:)).^2, '-k.'), hold on
plot(noise_var_mat2_dd(:), '-b.')
plot(noise_var_mat_fulltap_dd(:), '-r.'), hold off
title('dd domain (mmse-eq in dd domain)'), legend('real noise', 'est. noise (one-tap)', 'est. noise (full-tap)')
subplot(2, 1, 2)
plot(abs(tx_sym_rbs_tf(:)-rx_sym_rbs_usr_eq2_tf(:)).^2, '-k.'), hold on
plot(noise_var_mat2_tf(:), '-b.')
plot(noise_var_mat_fulltap_tf(:), '-r.')
plot(noise_var_mat2_alt_tf(:), '-m.'), hold off
title('tf domain (mmse-eq in dd domain)'), legend('real noise', 'est. noise (one-tap)', 'est. noise (full-tap)', 'est. noise (alt)')

% 3. zf-eq in tf
rx_sym_rbs_usr_eq3_tf = rx_sym_rbs_usr_tf./ch_est_rbs_usr_tf;
rx_sym_rbs_usr_eq3_dd = sqrt(nw_num.num_subc_usr/nw_num.num_ofdmsym_usr)*fft(ifft(rx_sym_rbs_usr_eq3_tf, [], 1), [], 2);

% 3.1. one-tap noise in tf
noise_var_mat3_tf = noise_var_mat1_tf;

% 3.2. one-tap noise in dd
noise_var_mat3_dd = noise_var_mat1_dd;

% 3.3 one-tap noise in tf alternative
ch_inv_dd = sqrt(nw_num.num_subc_usr/nw_num.num_ofdmsym_usr)*fft(ifft(ch_inv_tf, [], 1), [], 2);
ch_inv_alt_dd = reshape(sum(gen_eff_ch(ch_inv_dd), 2), nw_num.num_delay_usr, nw_num.num_doppler_usr);
% ch_inv_alt_dd = sum(ch_inv_dd, 'all')*ones(nw_num.num_delay_usr, nw_num.num_doppler_usr)/sqrt(numel(ch_inv_dd));    % simplified
noise_var_mat3_alt_dd = noise_var*(abs(ch_inv_alt_dd).^2);

% 4. mmse-eq in tf
ch_est_rbs_usr_mmse_tf = conj(ch_est_rbs_usr_tf)./(noise_var+abs(ch_est_rbs_usr_tf).^2);
rx_sym_rbs_usr_eq4_tf = ch_est_rbs_usr_mmse_tf.*rx_sym_rbs_usr_tf;
rx_sym_rbs_usr_eq4_dd = sqrt(nw_num.num_subc_usr/nw_num.num_ofdmsym_usr)*fft(ifft(rx_sym_rbs_usr_eq4_tf, [], 1), [], 2);

% 4.1. one-tap noise in tf
ch_inv_tf = ch_est_rbs_usr_mmse_tf;
noise_var_mat4_tf = noise_var*(abs(ch_inv_tf).^2);

% 4.2. one-tap noise in dd
ch_inv_dd = reshape(sum(sfft_mtx*diag(ch_est_rbs_usr_mmse_tf(:))*isfft_mtx, 2), nw_num.num_delay_usr, nw_num.num_doppler_usr);
noise_var_mat4_dd = noise_var*(abs(ch_inv_dd).^2);

% 4.3. one-tap noise in dd alternative (same with 4.1)
ch_mmse_inv_dd = sqrt(nw_num.num_subc_usr/nw_num.num_ofdmsym_usr)*fft(ifft(ch_inv_tf, [], 1), [], 2);
ch_inv_alt_dd = reshape(sum(gen_eff_ch(ch_mmse_inv_dd), 2), nw_num.num_delay_usr, nw_num.num_doppler_usr);
noise_var_mat4_alt_dd = noise_var*(abs(ch_inv_alt_dd).^2);

% plot zf vs. mmse in tf (equalized in tf)
figure
subplot(3, 1, 1)
plot(real(tx_sym_rbs_tf(:)), '-k.'), hold on
plot(real(rx_sym_rbs_usr_eq3_tf(:)), '-b.')
plot(real(rx_sym_rbs_usr_eq4_tf(:)), ':r.'), hold off
title('real'), legend('tx', 'zf-eq', 'mmse-eq rx')
subplot(3, 1, 2)
plot(imag(tx_sym_rbs_tf(:)), '-k.'), hold on
plot(imag(rx_sym_rbs_usr_eq3_tf(:)), '-b.')
plot(imag(rx_sym_rbs_usr_eq4_tf(:)), ':r.'), hold off
title('imag'), legend('tx', 'zf-eq', 'mmse-eq rx')
subplot(3, 1, 3)
plot(abs(tx_sym_rbs_tf(:)-rx_sym_rbs_usr_eq3_tf(:)).^2, '-b.'), hold on
plot(abs(tx_sym_rbs_tf(:)-rx_sym_rbs_usr_eq4_tf(:)).^2, ':r.'), hold off
title('tx est. error'), legend('zf-eq', 'mmse-eq')

% plot zf vs. mmse in dd (equalized in tf)
figure
subplot(3, 1, 1)
plot(real(tx_sym_rbs_dd(:)), '-k.'), hold on
plot(real(rx_sym_rbs_usr_eq3_dd(:)), '-b.')
plot(real(rx_sym_rbs_usr_eq4_dd(:)), ':r.'), hold off
title('real'), legend('tx', 'zf-eq rx', 'mmse-eq rx')
subplot(3, 1, 2)
plot(imag(tx_sym_rbs_dd(:)), '-k.'), hold on
plot(imag(rx_sym_rbs_usr_eq3_dd(:)), '-b.')
plot(imag(rx_sym_rbs_usr_eq4_dd(:)), ':r.'), hold off
title('imag'), legend('tx', 'zf-eq rx', 'mmse-eq rx')
subplot(3, 1, 3)
plot(abs(tx_sym_rbs_dd(:)-rx_sym_rbs_usr_eq3_dd(:)).^2, '-b.'), hold on
plot(abs(tx_sym_rbs_dd(:)-rx_sym_rbs_usr_eq4_dd(:)).^2, ':r.'), hold off
title('tx est. error'), legend('zf-eq', 'mmse-eq')

% plot constellation
figure
scatter(real(tx_sym_rbs_dd(:)), imag(tx_sym_rbs_dd(:)), 'bo'), hold on
scatter(real(rx_sym_rbs_usr_eq3_dd(:)), imag(rx_sym_rbs_usr_eq3_dd(:)), 'rx'), hold off
title('equalized in tf'), legend('tx', 'rx (zf-eq)')

% plot noise variance
figure
subplot(2, 1, 1)
plot(abs(tx_sym_rbs_dd(:)-rx_sym_rbs_usr_eq3_dd(:)).^2, '-k.'), hold on
plot(noise_var_mat3_dd(:), '-b.')
plot(noise_var_mat_fulltap_dd(:), '-r.')
plot(noise_var_mat3_alt_dd(:), '-m.'), hold off
title('dd domain (zf-eq in tf domain)'), legend('real noise', 'est. noise (one-tap)', 'est. noise (full-tap)', 'est. noise (one-tap alt.)')
subplot(2, 1, 2)
plot(abs(tx_sym_rbs_tf(:)-rx_sym_rbs_usr_eq3_tf(:)).^2, '-k.'), hold on
plot(noise_var_mat3_tf(:), '-b.')
plot(noise_var_mat_fulltap_tf(:), '-r.'), hold off
title('tf domain (zf-eq in tf domain)'), legend('real noise', 'est. noise (one-tap)', 'est. noise (full-tap)')

% plot noise variance
figure
subplot(2, 1, 1)
plot(abs(tx_sym_rbs_dd(:)-rx_sym_rbs_usr_eq4_dd(:)).^2, '-k.'), hold on
plot(noise_var_mat4_dd(:), '-b.')
plot(noise_var_mat_fulltap_dd(:), '-r.')
plot(noise_var_mat4_alt_dd(:), '-m.'), hold off
title('dd domain (mmse-eq in tf domain)'), legend('real noise', 'est. noise (one-tap)', 'est. noise (full-tap)', 'est. noise (one-tap alt.)')
subplot(2, 1, 2)
plot(abs(tx_sym_rbs_tf(:)-rx_sym_rbs_usr_eq4_tf(:)).^2, '-k.'), hold on
plot(noise_var_mat4_tf(:), '-b.'),
plot(noise_var_mat_fulltap_tf(:), '-r.'), hold off
title('tf domain (mmse-eq in ft domain)'), legend('real noise', 'est. noise (one-tap)', 'est. noise (full-tap)')
