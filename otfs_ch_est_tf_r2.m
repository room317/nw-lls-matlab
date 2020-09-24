function ch_est_rbs_tf = otfs_ch_est_tf_r2(tx_sym_rbs_tf, ch_out_time_serial, list_subc_usr, list_ofdmsym_usr, num, chest_option)

% perfect channel estimation including inter-subcarrier interference
if strcmp(chest_option, 'perfect')
    
    % demap channel input
    ch_in_sym_rbs_sfft = tx_sym_rbs_tf;
    
    % demap channel output
    ch_out_ofdmsym = reshape(ch_out_time_serial, num.num_fft+num.num_cp, num.num_ofdmsym);
    ch_out_sym_nfft = (1/sqrt(num.num_fft)) * fft(ch_out_ofdmsym(num.num_cp+1:end, :), [], 1);
    ch_out_sym_nfft_shift = fftshift(ch_out_sym_nfft, 1);
    ch_out_sym_bw = ch_out_sym_nfft_shift(num.num_fft/2-num.num_subc_bw/2+1:num.num_fft/2+num.num_subc_bw/2, :);
    ch_out_sym_rbs_sfft = ch_out_sym_bw(list_subc_usr, list_ofdmsym_usr);      % temporary
    
    % estimate channels
    ch_perfect_rbs_zf = ch_out_sym_rbs_sfft ./ ch_in_sym_rbs_sfft;  % estimate channel
    ch_sat = 100;   % ch_est saturation (zero-division prevention)
    ch_perfect_rbs_sat = complex(sign(real(ch_perfect_rbs_zf)).*min(abs(real(ch_perfect_rbs_zf)),ch_sat), ...
        sign(imag(ch_perfect_rbs_zf)).*min(abs(imag(ch_perfect_rbs_zf)),ch_sat));
    
    % output
    ch_est_rbs_tf = ch_perfect_rbs_sat;
    
else
    error('This function only supports ''perfect channel'' estimation.')
end

% % dump variables
% assignin('base', 'rx_sym_dd_ndft', rx_sym_dd_ndft);
% assignin('base', 'pilot_sym', pilot_sym);
% assignin('base', 'ch_est_dd_ndft_raw', ch_est_dd_ndft_raw);
% assignin('base', 'ch_est_dd_ndft_smooth', ch_est_dd_ndft_smooth);
% 
% % plot
% figure
% subplot(1, 2, 1), mesh(1:num.num_ofdmsym_per_usrfrm, 1:num.ndft, real(ch_est_dd_ndft_raw)), title('real(raw channel estimation) in dd-domain')
% subplot(1, 2, 2), mesh(1:num.num_ofdmsym_per_usrfrm, 1:num.ndft, imag(ch_est_dd_ndft_raw)), title('imag(raw channel estimation) in dd-domain')
% figure
% subplot(1, 2, 1), mesh(1:num.num_ofdmsym_per_usrfrm, 1:num.ndft, real(ch_est_dd_ndft_smooth)), title('real(smoothed channel estimation) in dd-domain')
% subplot(1, 2, 2), mesh(1:num.num_ofdmsym_per_usrfrm, 1:num.ndft, imag(ch_est_dd_ndft_smooth)), title('imag(smoothed channel estimation) in dd-domain')
% pause

end
