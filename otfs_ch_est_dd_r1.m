% len_buff: 2-dimensional

% map plan example (d: data, p: pilot, -: pilot, .: guard)
%   (11)data       :  d d d d d d d d d d d d d d d d
%   (10)data       :  d d d d d d d d d d d d d d d d
%    (9)data&pilot :  d d d d . . . . . . . . d d d d
%    (8)data&pilot :  d d d d . . . . . . . . d d d d
%    (7)data&pilot :  d d d d . . - - - - . . d d d d
%    (6)data&pilot :  d d d d . . - - p - . . d d d d
%    (5)data&guard :  d d d d . . - - - - . . d d d d
%    (4)data&guard :  d d d d . . - - - - . . d d d d
%    (3)data&guard :  d d d d . . . . . . . . d d d d
%    (2)data&guard :  d d d d . . . . . . . . d d d d
%    (1)data       :  d d d d d d d d d d d d d d d d
%    (0)data       :  d d d d d d d d d d d d d d d d

function ch_est_rbs_dd = otfs_ch_est_dd_r1(rx_sym_rbs_dd, num, test_option)

% demap pilot symbols
[~, rx_sym_pilot_subfrm] = otfs_sym_demap_r2(rx_sym_rbs_dd, num, test_option);

% set resource block center index
idx_delay_usr = floor(num.num_delay_usr/2)+1;
idx_doppler_usr = floor(num.num_doppler_usr/2)+1;

% despread pilot symbols
if test_option.otfs_map_plan == 1 || test_option.otfs_map_plan == 2 || test_option.otfs_map_plan == 3   % use impulse pilot
    
    % remove guard
    rx_sym_pilot_guard_removed = rx_sym_pilot_subfrm(num.num_delay_guard_usr+1:end-num.num_delay_guard_usr, num.num_doppler_guard_usr+1:end-num.num_doppler_guard_usr);
    
    % find dd-domain channel impulse response
    rx_sym_pilot_imp_resp = rx_sym_pilot_guard_removed;
    
    % set pilot resource center index
    idx_delay_pilot_usr = floor((num.num_delay_pilot_usr-2*num.num_delay_guard_usr)/2)+1;
    idx_doppler_pilot_usr = floor((num.num_doppler_pilot_usr-2*num.num_doppler_guard_usr)/2)+1;
    
    % set scaling factor
    scale_ch_est_rbs = sqrt((num.num_delay_usr*num.num_doppler_usr)/(num.num_delay_pilot_usr*num.num_doppler_pilot_usr));
    
elseif test_option.otfs_map_plan == 4       % use pilot sequence (guard removing -> despreading -> interpolating)
    
    % remove guard
    rx_sym_pilot_guard_removed = rx_sym_pilot_subfrm(num.num_delay_guard_usr+1:end-num.num_delay_guard_usr, num.num_doppler_guard_usr+1:end-num.num_doppler_guard_usr);
    
    % despread pilot
    rx_sym_pilot_despread = zeros(size(rx_sym_pilot_guard_removed)+[length(test_option.otfs_pilot_spread_seq)-1, 0]);
    for idx_doppler = 1:size(rx_sym_pilot_guard_removed, 2)
        rx_sym_pilot_despread(:, idx_doppler) = conv(rx_sym_pilot_guard_removed(:, idx_doppler), flipud(conj(test_option.otfs_pilot_spread_seq)));
    end
    
    % find dd-domain channel impulse response
    rx_sym_pilot_imp_resp = rx_sym_pilot_despread;
    
%     assignin('base', 'rx_sym_pilot_subfrm', rx_sym_pilot_subfrm)
%     assignin('base', 'rx_sym_pilot_guard_removed', rx_sym_pilot_guard_removed)
%     assignin('base', 'rx_sym_pilot_despread', rx_sym_pilot_despread)
%     assignin('base', 'rx_sym_pilot_imp_resp', rx_sym_pilot_imp_resp)
%     pause
    
%     % demap pilot
%     idx_pilot_seq = [ceil((length(test_option.pilot_spread_seq)-1)/2), num.num_doppler_guard_usr];
% %     rx_sym_pilot_subfrm = rx_sym_pilot_despread(idx_pilot_seqlength(test_option.otfs_pilot_spread_seq):length(test_option.otfs_pilot_spread_seq)+size(rx_sym_pilot_spread, 1)-1, :);
%     rx_sym_pilot_subfrm = sqrt(pwr_pilot)*rx_sym_pilot_despread(idx_pilot_seq(1)+1:idx_pilot_seq(1)+size(rx_sym_pilot_spread, 1), :);
    
    % set pilot resource center index
    idx_delay_pilot_usr = floor(length(test_option.otfs_pilot_spread_seq)/2)+floor((num.num_delay_pilot_usr-2*num.num_delay_guard_usr)/2)+1;
    idx_doppler_pilot_usr = floor((num.num_doppler_pilot_usr-2*num.num_doppler_guard_usr)/2)+1;
    
    % set scaling factor
    scale_ch_est_rbs = sqrt((num.num_delay_usr*num.num_doppler_usr)/(num.num_delay_pilot_usr*num.num_doppler_pilot_usr)* ...
        (length(test_option.otfs_pilot_spread_seq)/(2*length(test_option.otfs_pilot_spread_seq)-1))/ ...
        mean(abs(xcorr(test_option.otfs_pilot_spread_seq, test_option.otfs_pilot_spread_seq)).^2, 'all'));
    
elseif test_option.otfs_map_plan == 5       % use pilot sequence (guard removing -> despreading -> pruning -> interpolating)
    
    % remove guard
    rx_sym_pilot_guard_removed = rx_sym_pilot_subfrm(num.num_delay_guard_usr+1:end-num.num_delay_guard_usr, num.num_doppler_guard_usr+1:end-num.num_doppler_guard_usr);
    
    % despread pilot
    rx_sym_pilot_despread = zeros(size(rx_sym_pilot_guard_removed)+[length(test_option.otfs_pilot_spread_seq)-1, 0]);
    for idx_doppler = 1:size(rx_sym_pilot_guard_removed, 2)
        rx_sym_pilot_despread(:, idx_doppler) = conv(rx_sym_pilot_guard_removed(:, idx_doppler), flipud(conj(test_option.otfs_pilot_spread_seq)));
    end
    
    % find dd-domain channel impulse response (pruning)
    rx_sym_pilot_imp_resp = rx_sym_pilot_despread(floor(length(test_option.otfs_pilot_spread_seq)/2)+1: ...
        floor(length(test_option.otfs_pilot_spread_seq)/2)+(num.num_delay_pilot_usr-2*num.num_delay_guard_usr), :);
    
%     % demap pilot
%     idx_pilot_seq = [ceil((length(test_option.pilot_spread_seq)-1)/2), num.num_doppler_guard_usr];
% %     rx_sym_pilot_subfrm = rx_sym_pilot_despread(idx_pilot_seqlength(test_option.otfs_pilot_spread_seq):length(test_option.otfs_pilot_spread_seq)+size(rx_sym_pilot_spread, 1)-1, :);
%     rx_sym_pilot_subfrm = sqrt(pwr_pilot)*rx_sym_pilot_despread(idx_pilot_seq(1)+1:idx_pilot_seq(1)+size(rx_sym_pilot_spread, 1), :);
    
    % set pilot resource center index
    idx_delay_pilot_usr = floor((num.num_delay_pilot_usr-2*num.num_delay_guard_usr)/2)+1;
    idx_doppler_pilot_usr = floor((num.num_doppler_pilot_usr-2*num.num_doppler_guard_usr)/2)+1;
    
    % set scaling factor
    scale_ch_est_rbs = sqrt((num.num_delay_usr*num.num_doppler_usr)/(num.num_delay_pilot_usr*num.num_doppler_pilot_usr)* ...
        (length(test_option.otfs_pilot_spread_seq)/(2*length(test_option.otfs_pilot_spread_seq)-1))/ ...
        mean(abs(xcorr(test_option.otfs_pilot_spread_seq, test_option.otfs_pilot_spread_seq)).^2, 'all'));
    
elseif test_option.otfs_map_plan == 6       % use pilot sequence (despreading -> guard removing & pruning -> interpolating)
    
    % despread pilot
    rx_sym_pilot_despread = zeros(size(rx_sym_pilot_subfrm)+[length(test_option.otfs_pilot_spread_seq)-1, 0]);
    for idx_doppler = 1:size(rx_sym_pilot_subfrm, 2)
        rx_sym_pilot_despread(:, idx_doppler) = conv(rx_sym_pilot_subfrm(:, idx_doppler), flipud(conj(test_option.otfs_pilot_spread_seq)));
    end
    
    % find dd-domain channel impulse response (pruning)
    rx_sym_pilot_imp_resp = rx_sym_pilot_despread(num.num_delay_guard_usr+floor(length(test_option.otfs_pilot_spread_seq)/2)+1: ...
        num.num_delay_guard_usr+floor(length(test_option.otfs_pilot_spread_seq)/2)+(num.num_delay_pilot_usr-2*num.num_delay_guard_usr), :);
    
%     % demap pilot
%     idx_pilot_seq = [ceil((length(test_option.pilot_spread_seq)-1)/2), num.num_doppler_guard_usr];
% %     rx_sym_pilot_subfrm = rx_sym_pilot_despread(idx_pilot_seqlength(test_option.otfs_pilot_spread_seq):length(test_option.otfs_pilot_spread_seq)+size(rx_sym_pilot_spread, 1)-1, :);
%     rx_sym_pilot_subfrm = sqrt(pwr_pilot)*rx_sym_pilot_despread(idx_pilot_seq(1)+1:idx_pilot_seq(1)+size(rx_sym_pilot_spread, 1), :);
    
    % set pilot resource center index
    idx_delay_pilot_usr = floor((num.num_delay_pilot_usr-2*num.num_delay_guard_usr)/2)+1;
    idx_doppler_pilot_usr = floor((num.num_doppler_pilot_usr-2*num.num_doppler_guard_usr)/2)+1;
    
    % set scaling factor
    scale_ch_est_rbs = sqrt((num.num_delay_usr*num.num_doppler_usr)/(num.num_delay_pilot_usr*num.num_doppler_pilot_usr)* ...
        (length(test_option.otfs_pilot_spread_seq)/(2*length(test_option.otfs_pilot_spread_seq)-1))/ ...
        mean(abs(xcorr(test_option.otfs_pilot_spread_seq, test_option.otfs_pilot_spread_seq)).^2, 'all'));
    
else
    error('''otfs_map_plan'' must be one of these: {1, 2, 3, 4, 5, 6}')
end

% estimate channels
ch_est_rbs_raw = zeros(num.num_delay_usr, num.num_doppler_usr);
ch_est_rbs_raw(1:size(rx_sym_pilot_imp_resp, 1), 1:size(rx_sym_pilot_imp_resp, 2)) = rx_sym_pilot_imp_resp*scale_ch_est_rbs;
ch_est_rbs_raw_cntr = circshift(ch_est_rbs_raw, [idx_delay_usr-idx_delay_pilot_usr, idx_doppler_usr-idx_doppler_pilot_usr]);

% set pilot resource index
idx_ch_est_delay = (1:size(rx_sym_pilot_imp_resp, 1))+(idx_delay_usr-idx_delay_pilot_usr);
idx_ch_est_doppler = (1:size(rx_sym_pilot_imp_resp, 2))+(idx_doppler_usr-idx_doppler_pilot_usr);

% interpolate channels (linear)
if num.num_delay_data_usr > 0 || num.num_delay_guard_usr > 0
    idx_ch_est_delay_interp = [1, idx_ch_est_delay, num.num_delay_usr];
else
    idx_ch_est_delay_interp = idx_ch_est_delay;
end
if num.num_doppler_data_usr > 0 || num.num_doppler_guard_usr > 0
    idx_ch_est_doppler_interp = [1; idx_ch_est_doppler'; num.num_doppler_usr];
else
    idx_ch_est_doppler_interp = idx_ch_est_doppler';
end

ch_est_rbs_interp = interp2(idx_ch_est_doppler_interp, idx_ch_est_delay_interp, ...
    ch_est_rbs_raw_cntr(idx_ch_est_delay_interp, idx_ch_est_doppler_interp), (1:num.num_doppler_usr)', 1:num.num_delay_usr);

% circular shift channel to edge
ch_est_rbs_interp_edge = circshift(ch_est_rbs_interp, [-ceil(num.num_delay_usr/2), -ceil(num.num_doppler_usr/2)]);
% ch_est_rbs_interpl_edge = circshift(ch_est_rbs_raw_cntr, [-ceil(num.num_delay_usr/2), -ceil(num.num_doppler_usr/2)]);

% output
if test_option.ch_edge_interp
    
    % 2d inverse sfft for channel transformation
    ch_est_rbs_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(ch_est_rbs_interp_edge, [], 2), [], 1);
    
    % set guard parameter
    len_ch_edge_delay = 8;
    len_ch_edge_doppler = 0;
    
    % set edge interpolation index
    ch_est_rbs_tf_edge_interp = interp2((1+len_ch_edge_doppler:num.num_doppler_usr-len_ch_edge_doppler)', [1, (1+len_ch_edge_delay:num.num_delay_usr-len_ch_edge_delay), num.num_delay_usr], ...
        ch_est_rbs_tf(1+len_ch_edge_delay-1:num.num_delay_usr-len_ch_edge_delay+1, 1+len_ch_edge_doppler:num.num_doppler_usr-len_ch_edge_doppler), ...
        (1:num.num_doppler_usr)', 1:num.num_delay_usr);
    
    % 2d sfft for channel transformation
    ch_est_rbs_dd = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(ch_est_rbs_tf_edge_interp, [], 1), [], 2);
    
%     assignin('base', 'ch_est_rbs_tf_edge_interp', ch_est_rbs_tf_edge_interp)
%     assignin('base', 'ch_est_rbs_dd', ch_est_rbs_dd)
    
else
    ch_est_rbs_dd = ch_est_rbs_interp_edge;
end

% % dump variables
% assignin('base', 'idx_delay_pilot_usr', idx_delay_pilot_usr)
% assignin('base', 'idx_doppler_pilot_usr', idx_doppler_pilot_usr)
% assignin('base', 'rx_sym_rbs_dd', rx_sym_rbs_dd)
% assignin('base', 'rx_sym_pilot_subfrm', rx_sym_pilot_subfrm)
% assignin('base', 'rx_sym_pilot_guard_removed', rx_sym_pilot_guard_removed)
% assignin('base', 'rx_sym_pilot_despread', rx_sym_pilot_despread)
% assignin('base', 'rx_sym_pilot_imp_resp', rx_sym_pilot_imp_resp)
% assignin('base', 'ch_est_rbs_raw', ch_est_rbs_raw)
% assignin('base', 'ch_est_rbs_raw_cntr', ch_est_rbs_raw_cntr)
% assignin('base', 'idx_ch_est_delay', idx_ch_est_delay);
% assignin('base', 'idx_ch_est_doppler', idx_ch_est_doppler);
% assignin('base', 'idx_ch_est_delay_interp', idx_ch_est_delay_interp);
% assignin('base', 'idx_ch_est_doppler_interp', idx_ch_est_doppler_interp);
% assignin('base', 'ch_est_rbs_interpl', ch_est_rbs_interpl);
% assignin('base', 'ch_est_rbs_dd', ch_est_rbs_dd);
% pause
 
% % plot
% figure
% subplot(1, 2, 1), mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft, real(ch_est_dd_ndft_raw)), title('real(raw channel estimation) in dd-domain')
% subplot(1, 2, 2), mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft, imag(ch_est_dd_ndft_raw)), title('imag(raw channel estimation) in dd-domain')
% figure
% subplot(1, 2, 1), mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft, real(ch_est_dd_ndft_smooth)), title('real(smoothed channel estimation) in dd-domain')
% subplot(1, 2, 2), mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft, imag(ch_est_dd_ndft_smooth)), title('imag(smoothed channel estimation) in dd-domain')
% pause

end
