% demap otfs data and pilot qam symbols
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

function [rx_sym_data_usrfrm, rx_sym_pilot1_usrfrm, rx_sym_pilot2_usrfrm] = otfs_sym_demap_r3(rx_sym_rbs_dd, num, chest_option, test_option)

% set power for normalization
if strcmp(chest_option, 'dd_impulse')       % use impulse pilot
    % set power for papr reduction
    % 2d-impulse pilot with reduced power
    % set power as 'time domain' pilot impulse not to exceed average data power for papr reduction
    if test_option.otfs_pilot_pwr_set
        pwr_data = ((num.num_delay_usr*num.num_doppler_usr)-((num.num_delay_pilot_usr-1)*num.num_doppler_pilot_usr))/ ...
            (num.num_delay_usr*num.num_doppler_usr);
        pwr_pilot = pwr_data*num.num_delay_pilot_usr;
    else
        % set power for normal 2d-impulse pilot
        pwr_data = 1;
        pwr_pilot = 1;
    end
else
    % set power for sequence spreading
    pwr_data = 1;
    pwr_pilot = 1;
end

% set index
idx_delay_usr = floor(num.num_delay_usr/2)+1;
idx_doppler_usr = floor(num.num_doppler_usr/2)+1;
idx_delay_pilot_usr = num.num_delay_pilot_half_usr+1;       % pilot: even
idx_doppler_pilot_usr = num.num_doppler_pilot_half_usr+1;   % pilot: even

% relocate symbols according to map plans
if strcmp(chest_option, 'dd_impulse') || ...
        strcmp(chest_option, 'dd_zc') || strcmp(chest_option, 'dd_random') || ...
        strcmp(chest_option, 'dd_golay_serial') || strcmp(chest_option, 'dd_golay_parallel')
    map_shift = (-1)*[idx_delay_usr-idx_delay_pilot_usr, idx_doppler_usr-idx_doppler_pilot_usr];
elseif strcmp(chest_option, 'dd_golay_diag')        % use golay sequence
    map_shift = (-1)*[idx_delay_usr-num.num_delay_pilot_usr, idx_doppler_usr-idx_doppler_pilot_usr];
else
    error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
end
rx_sym_map_base = circshift(rx_sym_rbs_dd, map_shift);

% demap data and pilot symbols from resource block
if strcmp(chest_option, 'dd_impulse') || ...
        strcmp(chest_option, 'dd_zc') || strcmp(chest_option, 'dd_random')
    % demap pilot from user resource block (guard included)
    rx_sym_pilot1_usrfrm = sqrt(pwr_pilot)*rx_sym_map_base(1:num.num_delay_pilot_usr, 1:num.num_doppler_pilot_usr);
    rx_sym_pilot2_usrfrm = [];
    
    % demap data from user resource block
    rx_sym_data1_usrfrm = rx_sym_map_base(1:num.num_delay_pilot_usr, num.num_doppler_pilot_usr+1:end);
    rx_sym_data2_usrfrm = rx_sym_map_base(num.num_delay_pilot_usr+1:end, :);
    rx_sym_data3_usrfrm = [];
elseif strcmp(chest_option, 'dd_golay_serial')
    % demap pilot from user resource block (guard included)
    if num.num_doppler_usr == num.num_doppler_pilot_usr             % full doppler spreading
        rx_sym_pilot1_usrfrm = sqrt(pwr_pilot)*circshift(rx_sym_map_base(1:num.num_delay_pilot_half_usr, 1:num.num_doppler_pilot_usr), [0, -floor(num.num_doppler_usr/4)-1+idx_doppler_usr]);
        rx_sym_pilot2_usrfrm = sqrt(pwr_pilot)*circshift(rx_sym_map_base(num.num_delay_pilot_half_usr+1:num.num_delay_pilot_usr, 1:num.num_doppler_pilot_usr), [0, -floor(num.num_doppler_usr/4)]);
    else
        rx_sym_pilot1_usrfrm = sqrt(pwr_pilot)*rx_sym_map_base(1:num.num_delay_pilot_half_usr, 1:num.num_doppler_pilot_usr);
        rx_sym_pilot2_usrfrm = sqrt(pwr_pilot)*rx_sym_map_base(num.num_delay_pilot_half_usr+1:num.num_delay_pilot_usr, 1:num.num_doppler_pilot_usr);
    end
    
    % demap data from user resource block
    rx_sym_data1_usrfrm = rx_sym_map_base(1:num.num_delay_pilot_usr, num.num_doppler_pilot_usr+1:end);
    rx_sym_data2_usrfrm = rx_sym_map_base(num.num_delay_pilot_usr+1:end, :);
    rx_sym_data3_usrfrm = [];
elseif strcmp(chest_option, 'dd_golay_parallel')
    % demap pilot from user resource block (guard included)
    rx_sym_pilot1_usrfrm = sqrt(pwr_pilot)*rx_sym_map_base(1:num.num_delay_pilot_usr, 1:num.num_doppler_pilot_half_usr);
    rx_sym_pilot2_usrfrm = sqrt(pwr_pilot)*rx_sym_map_base(1:num.num_delay_pilot_usr, num.num_doppler_pilot_half_usr+1:num.num_doppler_pilot_usr);
    
    % demap data from user resource block
    rx_sym_data1_usrfrm = rx_sym_map_base(1:num.num_delay_pilot_usr, num.num_doppler_pilot_usr+1:end);
    rx_sym_data2_usrfrm = rx_sym_map_base(num.num_delay_pilot_usr+1:end, :);
    rx_sym_data3_usrfrm = [];
elseif strcmp(chest_option, 'dd_golay_diag')        % use golay sequence
    % demap pilot from user resource block (guard included)
    rx_sym_pilot1_usrfrm = sqrt(pwr_pilot)*rx_sym_map_base(1:num.num_delay_pilot_usr, 1:num.num_doppler_pilot_half_usr);
    rx_sym_pilot2_usrfrm = sqrt(pwr_pilot)*rx_sym_map_base(num.num_delay_pilot_usr+1:2*num.num_delay_pilot_usr, num.num_doppler_pilot_half_usr+1:num.num_doppler_pilot_usr);
    
    % demap data from user resource block
    rx_sym_data1_usrfrm = rx_sym_map_base(1:num.num_delay_pilot_usr, num.num_doppler_pilot_half_usr+1:num.num_doppler_pilot_usr);
    rx_sym_data2_usrfrm = rx_sym_map_base(num.num_delay_pilot_usr+1:2*num.num_delay_pilot_usr, 1:num.num_doppler_pilot_half_usr);
    rx_sym_data3_usrfrm = rx_sym_map_base(2*num.num_delay_pilot_usr+1:end, :);
else
    error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
end

% serialize data symbols
rx_sym_data_usrfrm = sqrt(pwr_data)*[rx_sym_data1_usrfrm(:); rx_sym_data2_usrfrm(:); rx_sym_data3_usrfrm(:)];

% % test
% fprintf('total rx avg pwr:%6.3f\n', sqrt(mean(abs(rx_sym_map_base).^2, 'all')))
% fprintf('rx data avg pwr:%6.3f\n', sqrt(mean(abs(rx_sym_data_usrfrm).^2, 'all')))
% if strcmp(chest_option, 'dd_impulse') || ...
%         strcmp(chest_option, 'dd_zc') || ...
%         strcmp(chest_option, 'dd_random')
%     fprintf('rx pilot avg pwr:%6.3f\n', sqrt(mean(abs(rx_sym_pilot1_usrfrm).^2, 'all')));
% elseif strcmp(chest_option, 'dd_golay_serial') || ...
%         strcmp(chest_option, 'dd_golay_parallel') || ...
%         strcmp(chest_option, 'dd_golay_diag')
%     fprintf('rx pilot1 avg pwr:%6.3f\n', sqrt(mean(abs(rx_sym_pilot1_usrfrm).^2, 'all')));
%     fprintf('rx pilot2 avg pwr:%6.3f\n', sqrt(mean(abs(rx_sym_pilot2_usrfrm).^2, 'all')));
% else
%     error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
% end
% figure, mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, abs(rx_sym_rbs_dd))   % dd real channel
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('RX Pilot Plan')
% if strcmp(chest_option, 'dd_impulse') || ...
%         strcmp(chest_option, 'dd_zc') || ...
%         strcmp(chest_option, 'dd_random')
%     figure, mesh(1:num.num_doppler_pilot_usr, 1:num.num_delay_pilot_usr, abs(rx_sym_pilot1_usrfrm))   % dd real channel
%     xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('RX Pilots')
% elseif strcmp(chest_option, 'dd_golay_serial')
%     figure, mesh(1:num.num_doppler_pilot_usr, 1:num.num_delay_pilot_usr, abs([rx_sym_pilot1_usrfrm; rx_sym_pilot2_usrfrm]))   % dd real channel
%     xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('RX Pilots')
% elseif strcmp(chest_option, 'dd_golay_parallel') || ...
%         strcmp(chest_option, 'dd_golay_diag')
%     figure, mesh(1:num.num_doppler_pilot_usr, 1:num.num_delay_pilot_usr, abs([rx_sym_pilot1_usrfrm, rx_sym_pilot2_usrfrm]))   % dd real channel
%     xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('RX Pilots')
% else
%     error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
% end

% % dump variables
% assignin('base', 'rx_sym_rbs_dd', rx_sym_rbs_dd)
% assignin('base', 'rx_sym_map_base', rx_sym_map_base)
% assignin('base', 'rx_sym_data_usrfrm', rx_sym_data_usrfrm)
% assignin('base', 'rx_sym_pilot1_usrfrm', rx_sym_pilot1_usrfrm)
% assignin('base', 'rx_sym_data1_usrfrm', rx_sym_data1_usrfrm)
% assignin('base', 'rx_sym_data2_usrfrm', rx_sym_data2_usrfrm)
% if strcmp(chest_option, 'dd_golay_serial') || ...
%         strcmp(chest_option, 'dd_golay_parallel')
%     assignin('base', 'rx_sym_pilot2_usrfrm', rx_sym_pilot2_usrfrm)
% elseif strcmp(chest_option, 'dd_golay_serial') || ...
%         strcmp(chest_option, 'dd_golay_parallel') || ...
%         strcmp(chest_option, 'dd_golay_diag')
%     assignin('base', 'rx_sym_pilot2_usrfrm', rx_sym_pilot2_usrfrm)
%     assignin('base', 'rx_sym_data3_usrfrm', rx_sym_data3_usrfrm)
% end
% pause

end
