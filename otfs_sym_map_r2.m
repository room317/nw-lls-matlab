% map otfs data and pilot qam symbols
% map plan 1 (d: data, p: pilot, -: pilot, .: guard in pilot)
%   (11)guard      :  . . . . d d d d d d d d . . . .
%   (10)guard      :  . . . . d d d d d d d d . . . .
%    (9)data       :  d d d d d d d d d d d d d d d d
%    (8)data       :  d d d d d d d d d d d d d d d d
%    (7)data       :  d d d d d d d d d d d d d d d d
%    (6)data       :  d d d d d d d d d d d d d d d d
%    (5)data&pilot :  . . . . d d d d d d d d . . . .
%    (4)data&pilot :  . . . . d d d d d d d d . . . .
%    (3)data&pilot :  - - . . d d d d d d d d . . - -
%    (2)data&pilot :  p - . . d d d d d d d d . . - -
%    (1)guard      :  - - . . d d d d d d d d . . - -
%    (0)guard      :  - - . . d d d d d d d d . . - -
% 
% map plan 2 (d: data, p: pilot, -: pilot, .: guard in pilot)
%   (11)guard      :  - - . . d d d d d d d d . . - -
%   (10)guard      :  - - . . d d d d d d d d . . - -
%    (9)data&guard :  . . . . d d d d d d d d . . . .
%    (8)data&guard :  . . . . d d d d d d d d . . . .
%    (7)data       :  d d d d d d d d d d d d d d d d
%    (6)data       :  d d d d d d d d d d d d d d d d
%    (5)data       :  d d d d d d d d d d d d d d d d
%    (4)data       :  d d d d d d d d d d d d d d d d
%    (3)data&guard :  . . . . d d d d d d d d . . . .
%    (2)data&guard :  . . . . d d d d d d d d . . . .
%    (1)data&pilot :  - - . . d d d d d d d d . . - -
%    (0)data&pilot :  p - . . d d d d d d d d . . - -
% 
% map plan 3 (d: data, p: pilot, -: pilot, .: guard in pilot)
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
% 
% map plan 4/5/6 (d: data, p: pilot (sequence), -: pilot, .: guard in pilot)
%   (11)data       :  d d d d d d d d d d d d d d d d
%   (10)data       :  d d d d d d d d d d d d d d d d
%    (9)data&pilot :  d d d d . . . . . . . . d d d d
%    (8)data&pilot :  d d d d . - - - p - - . d d d d
%    (7)data&pilot :  d d d d . - - - p - - . d d d d
%    (6)data&pilot :  d d d d . - - - p - - . d d d d
%    (5)data&guard :  d d d d . - - - p - - . d d d d
%    (4)data&guard :  d d d d . - - - p - - . d d d d
%    (3)data&guard :  d d d d . - - - - - - . d d d d
%    (2)data&guard :  d d d d . . . . . . . . d d d d
%    (1)data       :  d d d d d d d d d d d d d d d d
%    (0)data       :  d d d d d d d d d d d d d d d d

function tx_sym_rbs_dd = otfs_sym_map_r2(tx_sym_data_subfrm, num, test_option)

% set normalized data and pilot symbol power
if test_option.otfs_map_plan == 1 || test_option.otfs_map_plan == 2 || test_option.otfs_map_plan == 3   % use impulse pilot
    
    % set power for papr reduction
    if test_option.otfs_pilot_impulse_pwr_reduction
        % reduce 2d-impulse pilot symbol power and boost data symbol power
        % set power as 'time domain' pilot impulse not to exceed average data power for papr reduction
        pwr_data = (num.num_delay_usr*num.num_doppler_usr)/ ...
            ((num.num_delay_usr*num.num_doppler_usr)-(num.num_delay_pilot_usr*num.num_doppler_pilot_usr)+(num.num_doppler_pilot_usr.^2)/2);
        pwr_pilot = pwr_data*(num.num_doppler_pilot_usr.^2)/2;
    else
        % normal 2d-impulse pilot
        pwr_data = 1;
        pwr_pilot = num.num_delay_pilot_usr*num.num_doppler_pilot_usr;
    end
    
elseif test_option.otfs_map_plan == 4 || test_option.otfs_map_plan == 5   % use pilot sequence (spreaded by zadoff-chu)
    
    % set power for sequence spreading
    pwr_data = 1;
    pwr_pilot = num.num_delay_pilot_usr*num.num_doppler_pilot_usr/length(test_option.otfs_pilot_spread_seq);
    
elseif test_option.otfs_map_plan == 6   % use pilot sequence (all ones)
    
    % set power for sequence spreading
    pwr_data = 1;
    pwr_pilot = num.num_delay_pilot_usr*num.num_doppler_pilot_usr/length(test_option.otfs_pilot_spread_seq);
    
else
    error('''otfs_map_plan'' must be one of these: {1, 2, 3, 4, 5, 6}')
end

% set index
idx_delay_usr = floor(num.num_delay_usr/2)+1;
idx_doppler_usr = floor(num.num_doppler_usr/2)+1;
idx_delay_pilot_usr = floor(num.num_delay_pilot_usr/2)+1;
idx_doppler_pilot_usr = floor(num.num_doppler_pilot_usr/2)+1;

% set pilot resource block
if test_option.otfs_map_plan == 1 || test_option.otfs_map_plan == 2 || test_option.otfs_map_plan == 3   % use impulse pilot
    
    % generate pilot symbols
    tx_sym_pilot_subfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_usr);
    tx_sym_pilot_subfrm(idx_delay_pilot_usr, idx_doppler_pilot_usr) = sqrt(pwr_pilot);     % boosts pilots
    
elseif test_option.otfs_map_plan == 4 || test_option.otfs_map_plan == 5   % use pilot sequence (spreaded by zadoff-chu)
    
    % map spread sequence
%     idx_pilot_seq = [ceil((num.num_delay_pilot_usr+num.num_delay_guard_usr-length(test_option.pilot_spread_seq))/2), num.num_doppler_guard_usr];
%     tx_sym_pilot_guard_subfrm = zeros(num.num_delay_pilot_usr+num.num_delay_guard_usr, num.num_doppler_pilot_usr+num.num_doppler_guard_usr);
%     tx_sym_pilot_guard_subfrm(idx_pilot_seq(1)+1:idx_pilot_seq(1)+length(test_option.pilot_spread_seq), idx_pilot_seq(2)+1) = sqrt(pwr_pilot)*test_option.pilot_spread_seq;
    tx_sym_pilot_subfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_usr);
    tx_sym_pilot_subfrm(idx_delay_pilot_usr-ceil(length(test_option.otfs_pilot_spread_seq)/2)+1: ...
        idx_delay_pilot_usr-ceil(length(test_option.otfs_pilot_spread_seq)/2)+length(test_option.otfs_pilot_spread_seq), ...
        idx_doppler_pilot_usr) = ...
        sqrt(pwr_pilot)*test_option.otfs_pilot_spread_seq;
    
elseif test_option.otfs_map_plan == 6   % use pilot sequence (all ones)
    
    % map spread sequence
%     idx_pilot_seq = [ceil((num.num_delay_pilot_usr+num.num_delay_guard_usr-length(test_option.pilot_spread_seq))/2), num.num_doppler_guard_usr];
%     tx_sym_pilot_guard_subfrm = zeros(num.num_delay_pilot_usr+num.num_delay_guard_usr, num.num_doppler_pilot_usr+num.num_doppler_guard_usr);
%     tx_sym_pilot_guard_subfrm(idx_pilot_seq(1)+1:idx_pilot_seq(1)+length(test_option.pilot_spread_seq), idx_pilot_seq(2)+1) = sqrt(pwr_pilot)*test_option.pilot_spread_seq;
    tx_sym_pilot_subfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_usr);
    tx_sym_pilot_subfrm(idx_delay_pilot_usr-ceil(length(test_option.otfs_pilot_seq_ones)/2)+1: ...
        idx_delay_pilot_usr-ceil(length(test_option.otfs_pilot_seq_ones)/2)+length(test_option.otfs_pilot_seq_ones), ...
        idx_doppler_pilot_usr) = ...
        sqrt(pwr_pilot)*test_option.otfs_pilot_seq_ones;
    
else
    error('''otfs_map_plan'' must be one of these: {1, 2, 3, 4, 5, 6}')
end

% fprintf('pilot:%6.3f\n', sqrt(mean(abs(tx_sym_pilot_guard_subfrm).^2, 'all')))

% reshape data symbols
tx_sym_data1_subfrm = sqrt(pwr_data)*reshape(tx_sym_data_subfrm(1:num.num_delay_data_usr*num.num_doppler_usr), num.num_delay_data_usr, []);
tx_sym_data2_subfrm = sqrt(pwr_data)*reshape(tx_sym_data_subfrm(num.num_delay_data_usr*num.num_doppler_usr+1:end), num.num_delay_pilot_usr, []);

% fprintf('data1:%6.3f\n', sqrt(mean(abs(tx_sym_data1_subfrm).^2, 'all')))
% fprintf('data2:%6.3f\n', sqrt(mean(abs(tx_sym_data2_subfrm).^2, 'all')))

% map resource block
tx_sym_map_basic = [ ...
    tx_sym_pilot_subfrm, tx_sym_data2_subfrm;
    tx_sym_data1_subfrm];

% fprintf('total:%6.3f\n', sqrt(mean(abs(tx_sym_map_basic).^2, 'all')))
% pause

% relocate symbols according to map plans
if test_option.otfs_map_plan == 1
    map_shift = [-ceil((idx_delay_pilot_usr-1)/2), -(idx_doppler_pilot_usr-1)];
elseif test_option.otfs_map_plan == 2
    map_shift = [-(idx_delay_pilot_usr-1), -(idx_doppler_pilot_usr-1)];
elseif test_option.otfs_map_plan == 3 || test_option.otfs_map_plan == 4 || test_option.otfs_map_plan == 5 || test_option.otfs_map_plan == 6
    map_shift = [idx_delay_usr-idx_delay_pilot_usr, idx_doppler_usr-idx_doppler_pilot_usr];
else
    error('''otfs_map_plan'' must be one of these: {1, 2, 3, 4, 5, 6}.')
end
tx_sym_rbs_dd = circshift(tx_sym_map_basic, map_shift);

% % calculate single-tone position
% tone_pos = [num.num_delay_guard_usr, num.num_doppler_guard_usr]+map_shift;

% % dump variables
% assignin('base', 'tx_sym_data_subfrm', tx_sym_data_subfrm);
% assignin('base', 'tx_sym_pilot_subfrm', tx_sym_pilot_subfrm);
% assignin('base', 'tx_sym_data1_subfrm', tx_sym_data1_subfrm);
% assignin('base', 'tx_sym_data2_subfrm', tx_sym_data2_subfrm);
% assignin('base', 'tx_sym_map_basic', tx_sym_map_basic);
% assignin('base', 'tx_sym_rbs_dd', tx_sym_rbs_dd);
% pause

% % plot results
% figure, mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, abs(tx_sym_rbs_dd))   % dd real channel
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Pilot Plan')
% fprintf('tone position: (%d, %d)\n', tone_pos)
% pause

end
