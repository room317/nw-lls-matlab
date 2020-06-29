% map otfs data and pilot qam symbols
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

function tx_sym_rbs_dd = otfs_sym_map_r2(tx_sym_data_subfrm, num, map_plan)

% generate pilot symbols
tx_sym_pilot_subfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_usr);
tx_sym_pilot_subfrm(1, 1) = sqrt((num.num_delay_pilot_usr+num.num_delay_guard_usr)*(num.num_doppler_pilot_usr+num.num_doppler_guard_usr));  % boosts pilots

% pad guard symbols to pilot symbols
tx_sym_pilot_guard_subfrm = zeros(num.num_delay_pilot_usr+num.num_delay_guard_usr, num.num_doppler_pilot_usr+num.num_doppler_guard_usr);
tx_sym_pilot_guard_subfrm(num.num_delay_guard_usr+1:end, num.num_doppler_guard_usr+1:end) = tx_sym_pilot_subfrm;

% reshape data symbols
tx_sym_data1_subfrm = reshape(tx_sym_data_subfrm(1:num.num_delay_data_usr*num.num_doppler_usr), num.num_delay_data_usr, []);
tx_sym_data2_subfrm = reshape(tx_sym_data_subfrm(num.num_delay_data_usr*num.num_doppler_usr+1:end), num.num_delay_pilot_usr+num.num_delay_guard_usr, []);

% map resource block
tx_sym_map_basic = [ ...
    tx_sym_pilot_guard_subfrm, tx_sym_data2_subfrm;
    tx_sym_data1_subfrm];

% relocate symbols according to map plans
if map_plan == 1
    map_shift = [-ceil(num.num_delay_guard_usr/2), -num.num_doppler_guard_usr];
elseif map_plan == 2
    map_shift = [-num.num_delay_guard_usr, -num.num_doppler_guard_usr];
elseif map_plan == 3
    map_shift = [ceil(num.num_delay_data_usr/2), ceil(num.num_doppler_data_usr/2)];
else
    error('Map_plan must be one of these: {1, 2, 3}.')
end
tx_sym_rbs_dd = circshift(tx_sym_map_basic, map_shift);

% % calculate single-tone position
% tone_pos = [num.num_delay_guard_usr, num.num_doppler_guard_usr]+map_shift;

% % dump variables
% assignin('base', 'tx_sym_pilot_subfrm', tx_sym_pilot_subfrm);
% assignin('base', 'tx_sym_pilot_guard_subfrm', tx_sym_pilot_guard_subfrm);
% assignin('base', 'tx_sym_data1_subfrm', tx_sym_data1_subfrm);
% assignin('base', 'tx_sym_data2_subfrm', tx_sym_data2_subfrm);
% assignin('base', 'tx_sym_map_basic', tx_sym_map_basic);
% assignin('base', 'tx_sym_rbs', tx_sym_rbs);
% 
% % plot results
% figure, mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, abs(tx_sym_rbs))   % dd real channel
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Pilot Plan')
% fprintf('tone position: (%d, %d)\n', tone_pos)
% pause

end
