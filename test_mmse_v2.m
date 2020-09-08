% simulation for channel estimation techniequs using LS, LMMMSE, and
% computationally efficient LMMMSE methods. 
% Prepared by: Hiren Gami

% Ref: J J Van de Beek, "Synchronization and Channel Estimation in OFDM 
% systems", Ph.D thesis,Sept. 1998 

nsim = 1500;
ncp = 8;
nfft = 64;
esn0_db = 0:5:40;
snr = 10.^(esn0_db/10);
beta = 17/9;
m = 16;
ltap = 5;-

h_est_ls = zeros(1,length(esn0_db));
h_est_lmmse = zeros(1,length(esn0_db));
ht_est_lmmse = zeros(1,length(esn0_db));
htd_est_lmmse = zeros(1,length(esn0_db));
htq_est_lmmse = zeros(1,length(esn0_db));

for isnr = 1:length(esn0_db)
    disp(['EsN0dB is :', num2str(esn0_db(isnr))]);
    tic
    
    h_mse_ls = 0;
    h_mse_lmmse=0; 
    ht_mse_lmmse =0;
    htd_mse_lmmse=0;
    htq_mse_lmmse =0;
    
    for isim = 1:nsim
        
        % random channel taps
        g = complex(randn(ltap, 1), randn(ltap, 1));
        g_norm = g/norm(g);
        h = fft(g_norm, nfft);
        
        % generation of symbol
        x_bit = randi([0 m-1], nfft, 1);
%         xf = qammod(x_bit, m, 'UnitAveragePower', true); % normalizing symbol power
        xf = qammod(x_bit, m)/sqrt(10); % normalizing symbol power (for time saving)
        xt = ifft(xf)*sqrt(nfft);
        xt_cp = xt([nfft-ncp+1:nfft, 1:nfft], :);
        
        % channel convolution and awgn
        xt_fading = conv(xt_cp, g_norm);
        nt = complex(randn(nfft+ncp+ltap-1, 1), randn(nfft+ncp+ltap-1, 1));
        n0 = 10^(-esn0_db(isnr)/10);
        yt_cp =  xt_fading + sqrt(n0/2)*nt;
        
        % receiver processing
        yt = yt_cp(ncp+1:ncp+nfft);
        yf = fft(yt)/sqrt(nfft);
        
        % frequency doimain ls channel estimation 
        h_ls = yf./xf;
        h_mse_ls = h_mse_ls+((h-h_ls)'*(h-h_ls))/nfft;
        
        % frequency domain lmmse estimation
        r_hh = h*h';
        wh = r_hh/(r_hh+(beta/snr(isnr))*eye(nfft));
        h_lmmse = wh*h_ls;
        h_mse_lmmse = h_mse_lmmse+((h-h_lmmse)'*(h-h_lmmse))/nfft;
        
        % time domain lmmse estimation
        g_ls = ifft(h_ls, nfft);
        r_gg = g_norm*g_norm';
        wg = r_gg/(r_gg+(beta/snr(isnr))*eye(ltap));
        g_lmmse = wg*g_ls(1:ltap);
        ht_lmmse = fft(g_lmmse, nfft);
        ht_mse_lmmse = ht_mse_lmmse+((h-ht_lmmse)'*(h-ht_lmmse))/nfft;

        % time domain lmmse estimation - ignoring channel covariance
        rd_gg = diag(g_norm.*conj(g_norm));
        wgd = rd_gg/(rd_gg+(beta/snr(isnr))*eye(ltap));
        gd_lmmse = wgd*g_ls(1:ltap);
        htd_lmmse = fft(gd_lmmse, nfft);
        htd_mse_lmmse = htd_mse_lmmse+((h-htd_lmmse)'*(h-htd_lmmse))/nfft;
        
        % time domain lmmse estimation - ignoring smoothing matrix
        htq_lmmse = fft(g_ls, nfft);
        htq_mse_lmmse = htq_mse_lmmse+((h-htq_lmmse)'*(h-htq_lmmse))/nfft;
         
    end
    
    h_est_ls(isnr) = h_mse_ls/nsim;
    h_est_lmmse(isnr) = h_mse_lmmse/nsim;
    ht_est_lmmse(isnr) = ht_mse_lmmse/nsim;
    htd_est_lmmse(isnr) = htd_mse_lmmse/nsim;
    htq_est_lmmse(isnr) = htq_mse_lmmse/nsim;
    
    toc
end

% channel estimation 
figure
semilogy(esn0_db, h_est_ls, 'r', 'LineWidth', 2), hold on
semilogy(esn0_db, h_est_lmmse, 'k', 'LineWidth', 2);
semilogy(esn0_db, ht_est_lmmse, 'g', 'LineWidth', 2);
semilogy(esn0_db, htd_est_lmmse, 'm', 'LineWidth', 2);
semilogy(esn0_db, htq_est_lmmse, 'b', 'LineWidth', 2);

% theoratical bound calculation
semilogy(esn0_db, beta./snr, '-.r*', 'LineWidth', 2);
theory_lmmse = (1/nfft)*(beta./snr).*(1./(1+(beta./snr)));
semilogy(esn0_db, theory_lmmse, '-.k*', 'LineWidth', 2);

hold off
xlabel('esn0 (dB)'); ylabel('channel mse');
legend('ls','mmse', 'td lmmse','tdd lmmse','td qabs lmmse','theory-ls', 'theory-lmmse');
grid minor
