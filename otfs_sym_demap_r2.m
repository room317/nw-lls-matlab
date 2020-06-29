% demap otfs data and pilot qam symbols
% map plan 1 (d: data, p: pilot, -: pilot, +: guard)
%   (11)guard      :  + + + + d d d d d d d d + + + +
%   (10)guard      :  + + + + d d d d d d d d + + + +
%    (9)data       :  d d d d d d d d d d d d d d d d
%    (8)data       :  d d d d d d d d d d d d d d d d
%    (7)data       :  d d d d d d d d d d d d d d d d
%    (6)data       :  d d d d d d d d d d d d d d d d
%    (5)data&pilot :  - - - - d d d d d d d d + + + +
%    (4)data&pilot :  - - - - d d d d d d d d + + + +
%    (3)data&pilot :  - - - - d d d d d d d d + + + +
%    (2)data&pilot :  p - - - d d d d d d d d + + + +
%    (1)guard      :  + + + + d d d d d d d d + + + +
%    (0)guard      :  + + + + d d d d d d d d + + + +
% 
% map plan 2 (d: data, p: pilot, -: pilot, +: guard)
%   (11)guard      :  + + + + d d d d d d d d + + + +
%   (10)guard      :  + + + + d d d d d d d d + + + +
%    (9)guard      :  + + + + d d d d d d d d + + + +
%    (8)guard      :  + + + + d d d d d d d d + + + +
%    (7)data       :  d d d d d d d d d d d d d d d d
%    (6)data       :  d d d d d d d d d d d d d d d d
%    (5)data       :  d d d d d d d d d d d d d d d d
%    (4)data       :  d d d d d d d d d d d d d d d d
%    (3)data&pilot :  - - - - d d d d d d d d + + + +
%    (2)data&pilot :  - - - - d d d d d d d d + + + +
%    (1)data&pilot :  - - - - d d d d d d d d + + + +
%    (0)data&pilot :  p - - - d d d d d d d d + + + +
% 
% map plan 3 (d: data, p: pilot, -: pilot, +: guard)
%   (11)data       :  d d d d d d d d d d d d d d d d
%   (10)data       :  d d d d d d d d d d d d d d d d
%    (9)data&pilot :  d d d d + + + + - - - - d d d d
%    (8)data&pilot :  d d d d + + + + - - - - d d d d
%    (7)data&pilot :  d d d d + + + + - - - - d d d d
%    (6)data&pilot :  d d d d + + + + p - - - d d d d
%    (5)data&guard :  d d d d + + + + + + + + d d d d
%    (4)data&guard :  d d d d + + + + + + + + d d d d
%    (3)data&guard :  d d d d + + + + + + + + d d d d
%    (2)data&guard :  d d d d + + + + + + + + d d d d
%    (1)data       :  d d d d d d d d d d d d d d d d
%    (0)data       :  d d d d d d d d d d d d d d d d

function [rx_sym_data_subfrm, rx_sym_pilot_subfrm] = otfs_sym_demap_r2(rx_sym_rbs, num, map_plan)

% relocate symbols according to map plans
if map_plan == 1
    map_shift = [ceil(num.num_delay_guard_usr/2), num.num_doppler_guard_usr];
elseif map_plan == 2
    map_shift = [num.num_delay_guard_usr, num.num_doppler_guard_usr];
elseif map_plan == 3
    map_shift = [-ceil(num.num_delay_data_usr/2), -ceil(num.num_doppler_data_usr/2)];
else
    error('Map_plan must be one of these: {1, 2, 3}.')
end
rx_sym_map_basic = circshift(rx_sym_rbs, map_shift);

% demap pilot from user resource block (guard included)
rx_sym_pilot_subfrm = rx_sym_map_basic(1:num.num_delay_guard_usr+num.num_delay_pilot_usr, ...
    1:num.num_doppler_guard_usr+num.num_doppler_pilot_usr);

% demap data from user resource block
rx_sym_data1_subfrm = rx_sym_map_basic(num.num_delay_pilot_usr+num.num_delay_guard_usr+1:end, :);
rx_sym_data2_subfrm = rx_sym_map_basic(1:num.num_delay_pilot_usr+num.num_delay_guard_usr, ...
    num.num_doppler_pilot_usr+num.num_doppler_guard_usr+1:end);

% serialize data symbols
rx_sym_data_subfrm = [rx_sym_data1_subfrm(:); rx_sym_data2_subfrm(:)];

% dump variables
assignin('base', 'rx_sym_rbs', rx_sym_rbs);
assignin('base', 'rx_sym_map_basic', rx_sym_map_basic);
assignin('base', 'rx_sym_data1_subfrm', rx_sym_data1_subfrm);
assignin('base', 'rx_sym_data2_subfrm', rx_sym_data2_subfrm);
assignin('base', 'rx_sym_data_subfrm', rx_sym_data_subfrm);

% % plot results
% figure, mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, abs(rx_sym_map_basic))   % dd real channel
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Pilot Plan')
% figure, plot(abs(rx_sym_data_subfrm))
% xlabel('symbols'), ylabel('magnitude'), title('data symbols')
% pause

end
