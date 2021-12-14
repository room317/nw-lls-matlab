% pilot map plan
%   1. pilot plan: impulse (d: data, p: pilot, -: null, .: guard in pilot)
%       (11)data       :  d d d d d d d d d d d d d d d d
%       (10)data       :  d d d d d d d d d d d d d d d d
%        (9)data&pilot :  d d d d . . . . . . . . d d d d
%        (8)data&pilot :  d d d d . . . . . . . . d d d d
%        (7)data&pilot :  d d d d . . - - - - . . d d d d
%        (6)data&pilot :  d d d d . . - - p - . . d d d d
%        (5)data&guard :  d d d d . . - - - - . . d d d d
%        (4)data&guard :  d d d d . . - - - - . . d d d d
%        (3)data&guard :  d d d d . . . . . . . . d d d d
%        (2)data&guard :  d d d d . . . . . . . . d d d d
%        (1)data       :  d d d d d d d d d d d d d d d d
%        (0)data       :  d d d d d d d d d d d d d d d d
% 
%   2. map plan: random, zc (d: data, p: pilot (sequence), -: null, .: guard in pilot)
%       (11)data       :  d d d d d d d d d d d d d d d d
%       (10)data       :  d d d d d d d d d d d d d d d d
%        (9)data&pilot :  d d d d . . . . . . . . d d d d
%        (8)data&pilot :  d d d d . - - - p - - . d d d d
%        (7)data&pilot :  d d d d . - - - p - - . d d d d
%        (6)data&pilot :  d d d d . - - - p - - . d d d d
%        (5)data&guard :  d d d d . - - - p - - . d d d d
%        (4)data&guard :  d d d d . - - - p - - . d d d d
%        (3)data&guard :  d d d d . - - - - - - . d d d d
%        (2)data&guard :  d d d d . . . . . . . . d d d d
%        (1)data       :  d d d d d d d d d d d d d d d d
%        (0)data       :  d d d d d d d d d d d d d d d d
% 
%   3. map plan: golay_serial, (d: data, p/q: pilot (sequence), -: null, .: guard in pilot)
%       (11)data       :  d d d d d d d d d d d d d d d d
%       (10)data       :  d d d d d d d d d d d d d d d d
%        (9)data&pilot :  d d d d . . . . . . . . d d d d
%        (8)data&pilot :  d d d d . - - p - - - . d d d d
%        (7)data&pilot :  d d d d . - - p - - - . d d d d
%        (6)data&pilot :  d d d d . . . . . . . . d d d d
%        (5)data&guard :  d d d d . . . . . . . . d d d d
%        (4)data&guard :  d d d d . - - - q - - . d d d d
%        (3)data&guard :  d d d d . - - - q - - . d d d d
%        (2)data&guard :  d d d d . . . . . . . . d d d d
%        (1)data       :  d d d d d d d d d d d d d d d d
%        (0)data       :  d d d d d d d d d d d d d d d d
% 
%   4. map plan: golay_parallel (d: data, p/q: pilot (sequence), -: null, .: guard in pilot)
%       (11)data       :  d d d d d d d d d d d d d d d d
%       (10)data       :  d d d d d d d d d d d d d d d d
%        (9)data&pilot :  d d d d . . . . . . . . d d d d
%        (8)data&pilot :  d d d d . - p . . - q . d d d d
%        (7)data&pilot :  d d d d . - p . . - q . d d d d
%        (6)data&pilot :  d d d d . - p . . - q . d d d d
%        (5)data&guard :  d d d d . - p . . - q . d d d d
%        (4)data&guard :  d d d d . - p . . - q . d d d d
%        (3)data&guard :  d d d d . - p . . - q . d d d d
%        (2)data&guard :  d d d d . . . . . . . . d d d d
%        (1)data       :  d d d d d d d d d d d d d d d d
%        (0)data       :  d d d d d d d d d d d d d d d d
% 
%   5. map plan: golay_diag (d: data, p/q: pilot (sequence), -: null, .: guard in pilot)
%       (11)data       :  d d d d d d d d d d d d d d d d
%       (10)data       :  d d d d d d d d d d d d d d d d
%        (9)data&pilot :  d d d d d d d d . . . . d d d d
%        (8)data&pilot :  d d d d d d d d . - q . d d d d
%        (7)data&pilot :  d d d d d d d d . - q . d d d d
%        (6)data&pilot :  d d d d d d d d . . . . d d d d
%        (5)data&guard :  d d d d . . . . d d d d d d d d
%        (4)data&guard :  d d d d . - p . d d d d d d d d
%        (3)data&guard :  d d d d . - p . d d d d d d d d
%        (2)data&guard :  d d d d . . . . d d d d d d d d
%        (1)data       :  d d d d d d d d d d d d d d d d
%        (0)data       :  d d d d d d d d d d d d d d d d

function [ch_est_rbs_dd, ch_est_rbs_tf] = otfs_ch_est_dd_r3(rx_sym_rbs_dd, ch_real_rbs_usr_tf, noise_var, num, chest_option, test_option)

% demap pilot symbols
[~, rx_sym_pilot1_usrfrm, rx_sym_pilot2_usrfrm] = otfs_sym_demap_r3(rx_sym_rbs_dd, num, chest_option, test_option);

% set resource block center index
idx_delay_usr = floor(num.num_delay_usr/2)+1;
idx_doppler_usr = floor(num.num_doppler_usr/2)+1;

% remove guard
if strncmp(chest_option, 'dd_golay_', 9)        % use golay sequence pilot
    rx_sym_pilot1_guard_removed = rx_sym_pilot1_usrfrm(num.num_delay_guard_a_usr+1:end-num.num_delay_guard_b_usr, ...
        num.num_doppler_guard_usr+1:end-num.num_doppler_guard_usr);
    rx_sym_pilot2_guard_removed = rx_sym_pilot2_usrfrm(num.num_delay_guard_a_usr+1:end-num.num_delay_guard_b_usr, ...
        num.num_doppler_guard_usr+1:end-num.num_doppler_guard_usr);
else
    rx_sym_pilot1_guard_removed = rx_sym_pilot1_usrfrm(num.num_delay_guard_a_usr+1:end-num.num_delay_guard_b_usr, ...
        num.num_doppler_guard_usr+1:end-num.num_doppler_guard_usr);
    rx_sym_pilot2_guard_removed = [];
end

% extract channel info from rx symbols
if strcmp(chest_option, 'dd_impulse')       % use impulse pilot
    % find dd-domain channel impulse response
    rx_sym_pilot_imp_resp = rx_sym_pilot1_guard_removed;
    
    % set pilot resource center index
    idx_delay_pilot_usr = num.num_delay_pilot_half_usr-num.num_delay_guard_a_usr+1;           % pilot: even
    idx_doppler_pilot_usr = num.num_doppler_pilot_half_usr-num.num_doppler_guard_usr+1;     % pilot: even
    
    % set scaling factor
    scale_rx_ch_info = sqrt((num.num_delay_usr*num.num_doppler_usr)/(num.num_delay_pilot_usr*num.num_doppler_pilot_usr));
    
    % find channel info from rx signal
    rx_ch_info = rx_sym_pilot_imp_resp;
elseif strcmp(chest_option, 'dd_zc')        % use zadoff-chu sequence pilot
    % despread pilot
    rx_sym_pilot1_despread = zeros(size(rx_sym_pilot1_guard_removed)+[test_option.zc_seq_len-1, 0]);
    for idx_doppler = 1:size(rx_sym_pilot1_guard_removed, 2)
        rx_sym_pilot1_despread(:, idx_doppler) = conv(rx_sym_pilot1_guard_removed(:, idx_doppler), flipud(conj(test_option.zc_seq)));
    end
    
    % prune edges and set indices
    if test_option.ch_est_xcorr_prune
        % find dd-domain channel impulse response (pruning)
        rx_sym_pilot_imp_resp = rx_sym_pilot1_despread(floor(test_option.zc_seq_len/2)+1: ...
            floor(test_option.zc_seq_len/2)+(num.num_delay_pilot_usr-test_option.zc_seq_len), :);
        
        % set pilot resource center index
        idx_delay_pilot_usr = num.num_delay_pilot_half_usr-num.num_delay_guard_a_usr+1;           % pilot: even
        idx_doppler_pilot_usr = num.num_doppler_pilot_half_usr-num.num_doppler_guard_usr+1;     % pilot: even
    elseif size(rx_sym_pilot1_despread, 1) > num.num_delay_usr
        % find dd-domain channel impulse response (pruning)
        max_num_pilot_delay = min(num.num_delay_usr, size(rx_sym_pilot1_despread, 1)-floor(test_option.zc_seq_len/2));
        rx_sym_pilot_imp_resp = rx_sym_pilot1_despread(floor(test_option.zc_seq_len/2)+1: ...
            floor(test_option.zc_seq_len/2)+max_num_pilot_delay, :);
        
        % set pilot resource center index
        idx_delay_pilot_usr = num.num_delay_pilot_half_usr-num.num_delay_guard_a_usr+1;           % pilot: even
        idx_doppler_pilot_usr = num.num_doppler_pilot_half_usr-num.num_doppler_guard_usr+1;     % pilot: even
    else
        % find dd-domain channel impulse response
        rx_sym_pilot_imp_resp = rx_sym_pilot1_despread;
        
        % set pilot resource center index
        idx_delay_pilot_usr = floor(test_option.zc_seq_len/2)+num.num_delay_pilot_half_usr-num.num_delay_guard_a_usr+1;
        idx_doppler_pilot_usr = num.num_doppler_pilot_half_usr-num.num_doppler_guard_usr+1;
    end
    
    % check error
    if size(rx_sym_pilot_imp_resp, 1) > num.num_delay_usr
        error('Adjust spread sequence length and resource size.')
    end
    
    % set scaling factor
    scale_rx_ch_info = sqrt((num.num_delay_usr*num.num_doppler_usr)/(num.num_delay_pilot_usr*num.num_doppler_pilot_usr)* ...
        (test_option.zc_seq_len/(2*test_option.zc_seq_len-1))/mean(abs(xcorr(test_option.zc_seq, test_option.zc_seq)).^2, 'all'));
    
    % find channel info from rx signal
    rx_ch_info = rx_sym_pilot_imp_resp;
elseif strcmp(chest_option, 'dd_random')    % use random sequence pilot (generate toeplitz pilot matrix -> pseudo-inverse)
    % set parameters for tx pilot
    pwr_pilot = num.num_delay_pilot_usr*num.num_doppler_pilot_usr/test_option.rand_seq_len;
    idx_delay_pilot_usr = num.num_delay_pilot_half_usr+1;
    idx_doppler_pilot_usr = num.num_doppler_pilot_half_usr+1;
    
    % generate tx pilot
    tx_sym_pilot_dd = zeros(num.num_delay_usr, num.num_doppler_usr);
    tx_sym_pilot_dd(idx_delay_pilot_usr-ceil(test_option.rand_seq_len/2)+1: ...
        idx_delay_pilot_usr+floor(test_option.rand_seq_len/2), ...
        idx_doppler_pilot_usr) = ...
        sqrt(pwr_pilot)*test_option.rand_seq;
    
    % circshift tx pilot
    map_shift = [idx_delay_usr-idx_delay_pilot_usr, idx_doppler_usr-idx_doppler_pilot_usr];
    tx_sym_pilot_ctr_dd = circshift(tx_sym_pilot_dd, map_shift);
    
    % gen tx pilot in tf domain
    tx_sym_pilot_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(tx_sym_pilot_ctr_dd, [], 2), [], 1);
    
    % set pilot resource center index
    idx_delay_pilot_usr = num.num_delay_pilot_half_usr-num.num_delay_guard_a_usr+1;           % pilot: even
    idx_doppler_pilot_usr = num.num_doppler_pilot_half_usr-num.num_doppler_guard_usr+1;     % pilot: even
    
    % set scaling factor
    scale_rx_ch_info = 1;
%     scale_rx_ch_info = sqrt((num.num_delay_usr*num.num_doppler_usr)/(num.num_delay_pilot_usr*num.num_doppler_pilot_usr));
%     scale_rx_ch_info = sqrt(num.num_delay_usr*num.num_doppler_usr);
    
    % find channel info from rx signal
    rx_ch_info = rx_sym_pilot1_guard_removed;
elseif strncmp(chest_option, 'dd_golay_', 9)        % use zadoff-chu sequence pilot
    % despread pilot
    rx_sym_pilot1_despread = zeros(size(rx_sym_pilot1_guard_removed)+[test_option.golay_seq_len-1, 0]);
    rx_sym_pilot2_despread = zeros(size(rx_sym_pilot2_guard_removed)+[test_option.golay_seq_len-1, 0]);
    for idx_doppler = 1:size(rx_sym_pilot1_guard_removed, 2)
        rx_sym_pilot1_despread(:, idx_doppler) = conv(rx_sym_pilot1_guard_removed(:, idx_doppler), flipud(conj(test_option.golay_seq_a)));
        rx_sym_pilot2_despread(:, idx_doppler) = conv(rx_sym_pilot2_guard_removed(:, idx_doppler), flipud(conj(test_option.golay_seq_b)));
    end
    rx_sym_pilot_despread = rx_sym_pilot1_despread+rx_sym_pilot2_despread;
    
    % prune edges and set indices
    if test_option.ch_est_xcorr_prune
        % find dd-domain channel impulse response (pruning)
        rx_sym_pilot_imp_resp = rx_sym_pilot_despread(floor(test_option.golay_seq_len/2)+1: ...
            floor(test_option.golay_seq_len/2)+(size(rx_sym_pilot_despread, 1)-test_option.golay_seq_len), :);
        
        % set pilot resource center index
        idx_delay_pilot_usr = floor(size(rx_sym_pilot1_usrfrm, 1)/2)-num.num_delay_guard_a_usr+1;
        idx_doppler_pilot_usr = floor(size(rx_sym_pilot1_usrfrm, 2)/2)-num.num_doppler_guard_usr+1;     % pilot: even
    elseif size(rx_sym_pilot_despread, 1) > num.num_delay_usr
        % find dd-domain channel impulse response (pruning)
        max_num_pilot_delay = min(num.num_delay_usr, size(rx_sym_pilot_despread, 1)-floor(test_option.golay_seq_len/2));
        rx_sym_pilot_imp_resp = rx_sym_pilot_despread(floor(test_option.golay_seq_len/2)+1: ...
            floor(test_option.golay_seq_len/2)+max_num_pilot_delay, :);
        
        % set pilot resource center index
        idx_delay_pilot_usr = floor(size(rx_sym_pilot1_usrfrm, 1)/2)-num.num_delay_guard_a_usr+1;
        idx_doppler_pilot_usr = floor(size(rx_sym_pilot1_usrfrm, 2)/2)-num.num_doppler_guard_usr+1;     % pilot: even
    else
        % find dd-domain channel impulse response
        rx_sym_pilot_imp_resp = rx_sym_pilot_despread;
        
        % set pilot resource center index
        idx_delay_pilot_usr = floor(test_option.golay_seq_len/2)+floor(size(rx_sym_pilot1_usrfrm, 1)/2)-num.num_delay_guard_a_usr+1;
        idx_doppler_pilot_usr = floor(size(rx_sym_pilot1_usrfrm, 2)/2)-num.num_doppler_guard_usr+1;
    end
    
    % check error
    if size(rx_sym_pilot_imp_resp, 1) > num.num_delay_usr
        error('Adjust spread sequence length and resource size.')
    end
    
    % set scaling factor
    scale_rx_ch_info = sqrt((num.num_delay_usr*num.num_doppler_usr)/(num.num_delay_pilot_usr*num.num_doppler_pilot_usr)* ...
        (test_option.golay_seq_len/(2*test_option.golay_seq_len-1))/mean(abs(xcorr(test_option.golay_seq_a, test_option.golay_seq_a)).^2, 'all'));
    
    % find channel info from rx signal
    rx_ch_info = rx_sym_pilot_imp_resp;
else
    error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
end

% map channel info
rx_ch_info_rbs = zeros(num.num_delay_usr, num.num_doppler_usr);
rx_ch_info_rbs(1:size(rx_ch_info, 1), 1:size(rx_ch_info, 2)) = rx_ch_info*scale_rx_ch_info;
rx_ch_info_rbs_cntr = circshift(rx_ch_info_rbs, [idx_delay_usr-idx_delay_pilot_usr, idx_doppler_usr-idx_doppler_pilot_usr]);

% interpolate guard area
if test_option.ch_dd_guard_interp
    % set pilot resource index
    idx_ch_info_delay = (1:size(rx_ch_info, 1))+(idx_delay_usr-idx_delay_pilot_usr);
    idx_ch_info_doppler = (1:size(rx_ch_info, 2))+(idx_doppler_usr-idx_doppler_pilot_usr);
    
    % set interpolation index
    if idx_delay_usr > idx_delay_pilot_usr
        idx_ch_info_delay_interp = [1, idx_ch_info_delay, num.num_delay_usr];
    else
        idx_ch_info_delay_interp = idx_ch_info_delay;
    end
    if idx_doppler_usr > idx_doppler_pilot_usr
        idx_ch_info_doppler_interp = [1; idx_ch_info_doppler'; num.num_doppler_usr];
    else
        idx_ch_info_doppler_interp = idx_ch_info_doppler';
    end
    
    % interpolate channels (linear)
    rx_ch_info_rbs_interp_guard = interp2(idx_ch_info_doppler_interp, idx_ch_info_delay_interp, ...
        rx_ch_info_rbs_cntr(idx_ch_info_delay_interp, idx_ch_info_doppler_interp), (1:num.num_doppler_usr)', 1:num.num_delay_usr, 'linear');
else
    rx_ch_info_rbs_interp_guard = rx_ch_info_rbs_cntr;
end

% estimate and shift channel
if strcmp(chest_option, 'dd_random')    % use random sequence pilot
    % generate tf domain rx ch info
    rx_ch_info_rbs_interp_guard_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(rx_ch_info_rbs_interp_guard, [], 2), [], 1);
    
%     % estimate channnel with zf
%     ch_est_rbs_raw_tf = rx_ch_info_rbs_interp_guard_tf./tx_sym_pilot_tf;
    
    % estimate channel with mmse
    tx_sym_pilot_mmse_tf = conj(tx_sym_pilot_tf)./(noise_var+abs(tx_sym_pilot_tf).^2);
    ch_est_rbs_raw_tf = rx_ch_info_rbs_interp_guard_tf.*tx_sym_pilot_mmse_tf;
    
    % generate dd domain rx ch info
    ch_est_rbs_raw_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(ch_est_rbs_raw_tf, [], 1), [], 2);
else
    % circular shift channel to edge
    ch_est_rbs_raw_dd = circshift(rx_ch_info_rbs_interp_guard, [-idx_delay_usr+1, -idx_doppler_usr+1]);
    ch_est_rbs_raw_tf = [];
end

% output (time domain delay edge interpolation: no use)
if test_option.ch_tf_edge_interp
    % 2d inverse sfft for channel transformation
    if isempty(ch_est_rbs_raw_tf)
        ch_est_rbs_tf = sqrt(num.num_doppler_usr/num.num_delay_usr)*fft(ifft(ch_est_rbs_raw_dd, [], 2), [], 1);
    else
        ch_est_rbs_tf = ch_est_rbs_raw_tf;
    end
    
    % set new guard parameters
    len_ch_edge_delay_tf = 8;
    len_ch_edge_doppler_tf = 0;
    
    % set edge interpolation index
    ch_est_rbs_tf_interp_edge = interp2((1+len_ch_edge_doppler_tf:num.num_doppler_usr-len_ch_edge_doppler_tf)', ...
        [1, (1+len_ch_edge_delay_tf:num.num_delay_usr-len_ch_edge_delay_tf), num.num_delay_usr], ...
        ch_est_rbs_tf(1+len_ch_edge_delay_tf-1:num.num_delay_usr-len_ch_edge_delay_tf+1, ...
        1+len_ch_edge_doppler_tf:num.num_doppler_usr-len_ch_edge_doppler_tf), ...
        (1:num.num_doppler_usr)', 1:num.num_delay_usr);
    
    % 2d sfft for channel transformation
    ch_est_rbs_dd = sqrt(num.num_delay_usr/num.num_doppler_usr)*fft(ifft(ch_est_rbs_tf_interp_edge, [], 1), [], 2);
    ch_est_rbs_tf = ch_est_rbs_tf_interp_edge;
else
    ch_est_rbs_dd = ch_est_rbs_raw_dd;
    ch_est_rbs_tf = ch_est_rbs_raw_tf;
end

% % test
% if test_option.ch_mse
%     ch_real_rbs_usr_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(ch_real_rbs_usr_tf, [], 1), [], 2);
%     
%     fprintf('pilot power: %8.4f\n', mean(abs(ch_est_rbs_dd).^2, 'all'))
%     fprintf('ch est mse : %8.4f\n', mean(abs(ch_real_rbs_usr_dd-ch_est_rbs_dd).^2, 'all'))
%     figure
%     subplot(1, 3, 1), mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, abs(fftshift(fftshift(ch_real_rbs_usr_dd, 1), 2)))   % dd real channel
%     xlabel('doppler'), ylabel('delay'), zlabel('amplitude'), title('perfect channel'), axis_ref = axis;
%     subplot(1, 3, 2), mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, abs(fftshift(fftshift(ch_est_rbs_dd, 1), 2)))   % dd real channel
%     xlabel('doppler'), ylabel('delay'), zlabel('amplitude'), title('estimated channel'), axis(axis_ref)
%     subplot(1, 3, 3), mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, abs(fftshift(fftshift(ch_real_rbs_usr_dd-ch_est_rbs_dd, 1), 2)))   % dd real channel
%     xlabel('doppler'), ylabel('delay'), zlabel('rmse'), title('channel estimation error'), axis(axis_ref)
% end

% % dump
% assignin('base', 'rx_sym_rbs_dd', rx_sym_rbs_dd)
% assignin('base', 'rx_sym_pilot1_usrfrm', rx_sym_pilot1_usrfrm)
% assignin('base', 'rx_sym_pilot2_usrfrm', rx_sym_pilot2_usrfrm)
% assignin('base', 'idx_delay_usr', idx_delay_usr)
% assignin('base', 'idx_doppler_usr', idx_doppler_usr)
% assignin('base', 'rx_sym_pilot1_guard_removed', rx_sym_pilot1_guard_removed)
% assignin('base', 'rx_sym_pilot2_guard_removed', rx_sym_pilot2_guard_removed)
% assignin('base', 'idx_delay_pilot_usr', idx_delay_pilot_usr)
% assignin('base', 'idx_doppler_pilot_usr', idx_doppler_pilot_usr)
% assignin('base', 'rx_ch_info', rx_ch_info)
% if strcmp(chest_option, 'dd_random')
%     assignin('base', 'tx_sym_pilot_dd', tx_sym_pilot_dd)
%     assignin('base', 'tx_sym_pilot_tf', tx_sym_pilot_tf)
%     assignin('base', 'tx_sym_pilot_ctr_dd', tx_sym_pilot_ctr_dd)
%     assignin('base', 'rx_ch_info_rbs_interp_guard_tf', rx_ch_info_rbs_interp_guard_tf)
% elseif strncmp(chest_option, 'dd_golay_', 9)
%     assignin('base', 'rx_sym_pilot1_despread', rx_sym_pilot1_despread)
%     assignin('base', 'rx_sym_pilot2_despread', rx_sym_pilot2_despread)
%     assignin('base', 'rx_sym_pilot_despread', rx_sym_pilot_despread)
%     assignin('base', 'rx_sym_pilot_imp_resp', rx_sym_pilot_imp_resp)
% elseif strcmp(chest_option, 'dd_zc')
%     assignin('base', 'rx_sym_pilot1_despread', rx_sym_pilot1_despread)
%     assignin('base', 'rx_sym_pilot_imp_resp', rx_sym_pilot_imp_resp)
% elseif strcmp(chest_option, 'dd_impulse')
%     assignin('base', 'rx_sym_pilot_imp_resp', rx_sym_pilot_imp_resp)
% end
% assignin('base', 'rx_ch_info_rbs', rx_ch_info_rbs)
% assignin('base', 'rx_ch_info_rbs_cntr', rx_ch_info_rbs_cntr)
% if test_option.ch_dd_guard_interp
%     assignin('base', 'idx_ch_info_delay', idx_ch_info_delay)
%     assignin('base', 'idx_ch_info_doppler', idx_ch_info_doppler)
%     assignin('base', 'idx_ch_info_delay_interp', idx_ch_info_delay_interp)
%     assignin('base', 'idx_ch_info_doppler_interp', idx_ch_info_doppler_interp)
% end
% assignin('base', 'rx_ch_info_rbs_interp_guard', rx_ch_info_rbs_interp_guard)
% assignin('base', 'ch_est_rbs_raw_tf', ch_est_rbs_raw_tf)
% assignin('base', 'ch_est_rbs_raw_dd', ch_est_rbs_raw_dd)
% if test_option.ch_tf_edge_interp
%     assignin('base', 'ch_est_rbs_tf', ch_est_rbs_tf)
%     assignin('base', 'ch_est_rbs_tf_interp_edge', ch_est_rbs_tf_interp_edge)
% end
% assignin('base', 'ch_est_rbs_dd', ch_est_rbs_dd);
% if test_option.ch_mse
%     ch_real_rbs_usr_dd = sqrt(num.num_subc_usr/num.num_ofdmsym_usr)*fft(ifft(ch_real_rbs_usr_tf, [], 1), [], 2);
%     assignin('base', 'ch_real_rbs_usr_dd', ch_real_rbs_usr_dd)
%     assignin('base', 'ch_real_rbs_usr_tf', ch_real_rbs_usr_tf)
% end
% pause

end
