% demap otfs data and pilot qam symbols
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

function [rx_sym_data_usrfrm, rx_sym_pilot_usrfrm] = otfs_sym_demap_r2(rx_sym_rbs_dd, num, test_option)

% set power for normalization
if test_option.otfs_pilot_impulse_pwr_reduction && (test_option.otfs_map_plan == 1 || test_option.otfs_map_plan == 2 || test_option.otfs_map_plan == 3)   % use impulse pilot
    
    % set power for papr reduction
    % 2d-impulse pilot with reduced power
    % set power as 'time domain' pilot impulse not to exceed average data power for papr reduction
    pwr_data = ((num.num_delay_usr*num.num_doppler_usr)-((num.num_delay_pilot_usr-1)*num.num_doppler_pilot_usr))/ ...
        (num.num_delay_usr*num.num_doppler_usr);
    pwr_pilot = pwr_data*num.num_delay_pilot_usr;
else
    % set power for normal 2d-impulse pilot
    % set power for pilot sequence spreading
    pwr_data = 1;
    pwr_pilot = 1;      % 1/length(test_option.otfs_pilot_spread_seq);
end

% set index
idx_delay_usr = floor(num.num_delay_usr/2)+1;
idx_doppler_usr = floor(num.num_doppler_usr/2)+1;
idx_delay_pilot_usr = floor(num.num_delay_pilot_usr/2)+1;
idx_doppler_pilot_usr = floor(num.num_doppler_pilot_usr/2)+1;

% relocate symbols according to map plans
if test_option.otfs_map_plan == 1
    map_shift = [ceil((idx_delay_pilot_usr-1)/2), idx_doppler_pilot_usr-1];
elseif test_option.otfs_map_plan == 2
    map_shift = [idx_delay_pilot_usr-1, idx_doppler_pilot_usr-1];
elseif test_option.otfs_map_plan == 3 || test_option.otfs_map_plan == 4 || test_option.otfs_map_plan == 5 || test_option.otfs_map_plan == 6
    map_shift = [-(idx_delay_usr-idx_delay_pilot_usr), -(idx_doppler_usr-idx_doppler_pilot_usr)];
else
    error('''otfs_map_plan'' must be one of these: {1, 2, 3, 4}.')
end
rx_sym_map_basic = circshift(rx_sym_rbs_dd, map_shift);

% demap pilot from user resource block (guard included)
rx_sym_pilot_usrfrm = sqrt(pwr_pilot)*rx_sym_map_basic(1:num.num_delay_pilot_usr, 1:num.num_doppler_pilot_usr);

% % demap pilot from user resource block (guard included)
% if test_option.otfs_pilot_plan == 3
%     % demap spreaded pilot from user resource block (guard included)
%     rx_sym_pilot_spread = rx_sym_map_basic(1:num.num_delay_guard_usr+num.num_delay_pilot_usr, ...
%         1:num.num_doppler_guard_usr+num.num_doppler_pilot_usr);
%     
%     % despread pilot
%     rx_sym_pilot_despread = zeros(size(rx_sym_pilot_spread, 1)+length(test_option.otfs_pilot_spread_seq)-1, size(rx_sym_pilot_spread, 2));
%     for idx_doppler = 1:size(rx_sym_pilot_spread, 2)
%         rx_sym_pilot_despread(:, idx_doppler) = conv(rx_sym_pilot_spread(:, idx_doppler), flipud(conj(test_option.otfs_pilot_spread_seq)));
%     end
%     
%     % demap pilot
%     idx_pilot_seq = [ceil((length(test_option.pilot_spread_seq)-1)/2), num.num_doppler_guard_usr];
% %     rx_sym_pilot_usrfrm = rx_sym_pilot_despread(idx_pilot_seqlength(test_option.otfs_pilot_spread_seq):length(test_option.otfs_pilot_spread_seq)+size(rx_sym_pilot_spread, 1)-1, :);
%     rx_sym_pilot_usrfrm = sqrt(pwr_pilot)*rx_sym_pilot_despread(idx_pilot_seq(1)+1:idx_pilot_seq(1)+size(rx_sym_pilot_spread, 1), :);
% else
%     rx_sym_pilot_usrfrm = sqrt(pwr_pilot)*rx_sym_map_basic(1:num.num_delay_guard_usr+num.num_delay_pilot_usr, ...
%         1:num.num_doppler_guard_usr+num.num_doppler_pilot_usr);
% end

% assignin('base', 'rx_sym_rbs_dd', rx_sym_rbs_dd)
% assignin('base', 'rx_sym_map_basic', rx_sym_map_basic)
% assignin('base', 'rx_sym_pilot_usrfrm', rx_sym_pilot_usrfrm)
% pause

% demap data from user resource block
rx_sym_data1_usrfrm = rx_sym_map_basic(num.num_delay_pilot_usr+1:end, :);
rx_sym_data2_usrfrm = rx_sym_map_basic(1:num.num_delay_pilot_usr, ...
    num.num_doppler_pilot_usr+1:end);

% serialize data symbols
rx_sym_data_usrfrm = sqrt(pwr_data)*[rx_sym_data1_usrfrm(:); rx_sym_data2_usrfrm(:)];

% % dump variables
% assignin('base', 'rx_sym_rbs', rx_sym_rbs);
% assignin('base', 'rx_sym_map_basic', rx_sym_map_basic);
% assignin('base', 'rx_sym_data1_usrfrm', rx_sym_data1_usrfrm);
% assignin('base', 'rx_sym_data2_usrfrm', rx_sym_data2_usrfrm);
% assignin('base', 'rx_sym_data_usrfrm', rx_sym_data_usrfrm);

% % plot results
% figure, mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, abs(rx_sym_map_basic))   % dd real channel
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Pilot Plan')
% figure, plot(abs(rx_sym_data_usrfrm))
% xlabel('symbols'), ylabel('magnitude'), title('data symbols')
% pause

end
