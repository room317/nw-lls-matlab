% ofdm system with channel estimation
% based on mimo-ofdm wireless communications with matlab, 2010 john wiley & sons

% set numerical parameters 
num_ofdmsym = 14;
num_fft = 1024;
qam_size = 16;
subc_pilot = 6:6:num_fft;
pwr_pilot = 1;
num_cp = 72;

% set simulation parameter
snr_list = 0:4:40;
num_ch_tap = 10;

% generate bits
tx_bit = randi([0 qam_size-1], num_fft, num_ofdmsym);

% map qam symbols (baseband modulation)
tx_sym = qammod(tx_bit, qam_size);

% map data and pilots
subc_data = setxor(1:num_fft, subc_pilot);
tx_sym(subc_pilot, :) = pwr_pilot*tx_sym(subc_pilot, :);

% ifft
tx_ofdmsym = ifft(tx_sym);

% insert cp
tx_ofdmsym_cp = tx_ofdmsym([end-num_cp+1:end, 1:end], :);

% generate random fading channel (block fading)
ch = complex(randn(num_ch_tap, 1), randn(num_ch_tap, 1));
ch_norm = ch./norm(ch);
tx_ofdmsym_cp_fading = filter(ch_norm, 1, tx_ofdmsym_cp);
ch_freq = fft(ch_norm, num_fft);        % frequency-domain channel
ch_power_db = 20*log10(abs(ch_freq));   % true channel power in dB

% calculate rms delay spread
t0 = (0:num_ch_tap-1)';
t1_delay = ((t0.*ch_norm)'*ch_norm)/(ch_norm'*ch_norm);
t2_delay = ((t0.*ch_norm)'*(t0.*ch_norm))/(ch_norm'*ch_norm);
t_rms = sqrt(t2_delay-t1_delay^2);      % rms delay

% calculate common channel parameters
denom = 1i*2*pi*t_rms/num_fft;          % denominator of eq(6.16) p.192, mimo-ofdm wireless communications with matlab, 2010 john wiley & sons

% initialize variables
ch_ls_interp = zeros(num_fft, num_ofdmsym);
ser_awgn = zeros(1, length(snr_list));
ser_ls = zeros(1, length(snr_list));
ser_mmse = zeros(1, length(snr_list));

% receiver
for idx_snr = 1:length(snr_list)
    
    disp(['simulation @ ', num2str(snr_list(idx_snr)),' dB (', num2str(round(idx_snr/length(snr_list)*100)), ' %)'])
    
    % add noise
    snr_db = snr_list(idx_snr);
    snr_ebn0 = snr_db + 10*log10(log2(qam_size));          % eb/n0
    noise_var = 10^(snr_ebn0*0.1);
    rx_ofdmsym_awgn = awgn(tx_ofdmsym_cp, snr_ebn0, 'measured');
    rx_ofdmsym_fading = awgn(tx_ofdmsym_cp_fading, snr_ebn0, 'measured');
    
    % remove cyclic prefix
    rx_ofdmsym_awgn = rx_ofdmsym_awgn(num_cp+1:num_cp+num_fft, :);
    rx_ofdmsym_fading = rx_ofdmsym_fading(num_cp+1:num_cp+num_fft, :);
    
    % fft
    rx_sym_awgn = fft(rx_ofdmsym_awgn);
    rx_sym_fading = fft(rx_ofdmsym_fading);
    
    % extract pilots
    tx_sym_pilot = tx_sym(subc_pilot, :);               % tx pilots
    rx_sym_pilot = rx_sym_fading(subc_pilot, :);        % rx pilots
    
    % estimate channel (ls)
    ch_ls = rx_sym_pilot ./ tx_sym_pilot;
    
    % calculate channel/pilot cross-correlation
    k1 = repmat((0:num_fft-1)', 1, length(subc_pilot));
    k2 = repmat(subc_pilot, num_fft, 1);
    rf1 = 1./(1+denom*(k1-k2));
    rhp = rf1;
    
    % calculate pilot auto-correlation
    k3 = repmat(subc_pilot', 1, length(subc_pilot));
    k4 = repmat(subc_pilot, length(subc_pilot), 1);
    rf2 = 1./(1+denom*(k3-k4));
    rpp = rf2+eye(length(subc_pilot))/noise_var;
    
    % estimate channel (mmse)
    ch_mmse = rhp/rpp*ch_ls;
    
    % interpolate ls channel estimation
    for idx_ofdmsym = 1 : num_ofdmsym
        ch_ls_interp(:,idx_ofdmsym) = interp1(subc_pilot, ch_ls(:,idx_ofdmsym), 1:num_fft, 'linear', 'extrap');
    end
    
    % demapping 
    rx_bit_awgn = qamdemod(rx_sym_awgn, qam_size);                  % no fading
    rx_bit_ls = qamdemod(rx_sym_fading./ch_ls_interp, qam_size);    % ls estimation
    rx_bit_mmse = qamdemod(rx_sym_fading./ch_mmse, qam_size);       % mmse estimation
    
    % extract data bits
    tx_bit_data = tx_bit(subc_data, :);
    rx_bit_data_awgn = rx_bit_awgn(subc_data, :);
    rx_bit_data_ls = rx_bit_ls(subc_data, :);
    rx_bit_data_mmse = rx_bit_mmse(subc_data, :);
    
    % calculate ser
    [~, ser_awgn(idx_snr)] = symerr(tx_bit_data, rx_bit_data_awgn);
    [~, ser_ls(idx_snr)] = symerr(tx_bit_data, rx_bit_data_ls);
    [~, ser_mmse(idx_snr)] = symerr(tx_bit_data, rx_bit_data_mmse);
end

figure
semilogy(snr_list, ser_awgn, '-+'), hold on
semilogy(snr_list, ser_ls, '-o')
semilogy(snr_list, ser_mmse, '-s'), hold off
legend('awgn', 'ls ce', 'mmse ce');
grid minor

ch_pwr_db_ls = 20*log10(abs(ch_ls_interp));
ch_pwr_db_mmse = 20*log10(abs(ch_mmse));

figure
plot(ch_power_db(1:8:end), '+k', 'LineWidth', 3), hold on
plot(ch_pwr_db_ls((1:8:end), 1), 'or', 'LineWidth', 3)
plot(ch_pwr_db_mmse((1:8:end), 1), 'Sb', 'LineWidth', 1), hold off
title('channel power')
xlabel('samples')
ylabel('magnitudes')
legend('actual','ls estimate','mmse estimate')
grid minor

