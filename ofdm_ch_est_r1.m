% ofdm_ch_est_r1 estimates time-frequency domain channel.
% ch_est_rbs = ofdm_ch_est_r1(tx_sym_rbs, rx_sym_rbs, ch_out_time_serial, list_subc_usr, list_ofdmsym_usr, num, chest_option)
%   - ch_est_rbs: time-frequency channel estimation
%   - tx_sym_rbs: transmit signal for perfect channel estimation
%   - ch_out_time_serial: channel output signal for for perfert channel estimation
%   - list_subc_usr: subcarrier index vector of user resource 
%   - list_ofdmsym_usr: ofdm symbol index vector of user resource
%   - num: numerology of ofdm system
%   - chest_option: channel estimation options

function ch_est_rbs = ofdm_ch_est_r1(tx_sym_rbs, rx_sym_rbs, ch_out_time_serial, list_subc_usr, list_ofdmsym_usr, num, chest_option)

if strcmp(chest_option, 'tf_ltedown') || strcmp(chest_option, 'tf_lteup') || strcmp(chest_option, 'tf_nr')
    
    % set parameters
    ofdmsym_avg_win_size = 1;   % ofdm symbol average window size for channel estimation
    
    rx_sym_pilot_rbs = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
    rx_ofdmsym_pilot_interp_1d = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
    rx_ofdmsym_pilot_interp = zeros(num.num_subc_usr, num.num_ofdmsym_usr);
    for i = 1:num.num_ofdmsym_pilot_usr
        % demap pilots
        rx_sym_pilot_rbs(num.idx_subc_pilot_usr(:, i), num.idx_ofdmsym_pilot_usr(i)) = ...
            rx_sym_rbs(num.idx_subc_pilot_usr(:, i), num.idx_ofdmsym_pilot_usr(i)) .* ...
            conj(tx_sym_rbs(num.idx_subc_pilot_usr(:, i), num.idx_ofdmsym_pilot_usr(i)));
        
        % find symbol index within averaging window
        idx_ofdmsym_pilot_avg = ...
            num.idx_ofdmsym_pilot_usr( ...
            num.idx_ofdmsym_pilot_usr >= max(1, num.idx_ofdmsym_pilot_usr(i)-ofdmsym_avg_win_size) & ...
            num.idx_ofdmsym_pilot_usr <= num.idx_ofdmsym_pilot_usr(i));
        num_ofdmsym_pilot_avg = length(idx_ofdmsym_pilot_avg);
        
        % average pilots across symbols
        cnt_ofdmsym_pilot_avg = ones(num.num_subc_usr, 1);
        cnt_ofdmsym_pilot_avg(num.idx_subc_pilot_usr(:, i), 1) = ...
            sum(double(ismember(num.idx_subc_pilot_usr(:, i-num_ofdmsym_pilot_avg+1:i), num.idx_subc_pilot_usr(:, i))), 2);
        mask_subc_avg = zeros(num.num_subc_usr, 1);
        mask_subc_avg(num.idx_subc_pilot_usr(:, i), 1) = ones(num.num_subc_pilot_usr, 1);
        rx_ofdmsym_pilot_avg = sum(rx_sym_pilot_rbs(:, idx_ofdmsym_pilot_avg), 2) .* mask_subc_avg ./ cnt_ofdmsym_pilot_avg;
        
        % generate virtual pilots (ignored)
        rx_ofdmsym_pilot_vp = rx_ofdmsym_pilot_avg;
        
        % interpolate 1st dimension (along subcarriers)
        rx_ofdmsym_pilot_interp_1d(:, num.idx_ofdmsym_pilot_usr(i)) = ...
            interp1(num.idx_subc_pilot_usr(:, i), rx_ofdmsym_pilot_vp(num.idx_subc_pilot_usr(:, i), 1), ...
            1:num.num_subc_usr, 'linear', 'extrap');
        
    end
    
    % interpolate 2nd dimension (along ofdm symbols)
    for i = 1:num.num_subc_usr
        rx_ofdmsym_pilot_interp(i, :) = ...
            interp1(num.idx_ofdmsym_pilot_usr, rx_ofdmsym_pilot_interp_1d(i, num.idx_ofdmsym_pilot_usr), ...
            1:num.num_ofdmsym_usr, 'linear', 'extrap');
    end
    
    % output
    ch_est_rbs = rx_ofdmsym_pilot_interp;
    
elseif strcmp(chest_option, 'perfect')
    
    % demap channel input
    ch_in_sym_rbs = tx_sym_rbs;
    
    % demap channel output
    ch_out_ofdmsym = reshape(ch_out_time_serial, num.num_fft+num.num_cp, []);
    ch_out_sym_nfft = (1/sqrt(num.num_fft))*fft(ch_out_ofdmsym(num.num_cp+1:end, :), [], 1);
    ch_out_sym_nfft_shift = fftshift(ch_out_sym_nfft, 1);
    ch_out_sym_bw = ch_out_sym_nfft_shift(num.num_fft/2-num.num_subc_bw/2+1:num.num_fft/2+num.num_subc_bw/2, :);
    ch_out_sym_rbs = ch_out_sym_bw(list_subc_usr, list_ofdmsym_usr);
    
    % estimate channels
    ch_perfect_rbs_zf = ch_out_sym_rbs./ch_in_sym_rbs;  % estimate channel
    ch_sat = 100;   % ch_est saturation (zero-division prevention)
    ch_perfect_rbs_sat = complex(sign(real(ch_perfect_rbs_zf)).*min(abs(real(ch_perfect_rbs_zf)),ch_sat), ...
        sign(imag(ch_perfect_rbs_zf)).*min(abs(imag(ch_perfect_rbs_zf)),ch_sat));
    
    % output
    ch_est_rbs = ch_perfect_rbs_sat;
    
else
    error('This function supports following channel estimation options: {''tf_ltedown'', ''tf_lteup'', ''tf_nr'', ''perfect''}')
end

end
