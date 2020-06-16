% map otfs data and pilot qam symbols
% map plan 1
%   (11)guard      :  - - d d d d - -
%   (10)guard      :  - - d d d d - -
%    (9)data       :  d d d d d d d d
%    (8)data       :  d d d d d d d d
%    (7)data       :  d d d d d d d d
%    (6)data       :  d d d d d d d d
%    (5)data&pilot :  - - d d d d - -
%    (4)data&pilot :  - - d d d d - -
%    (3)data&pilot :  - - d d d d - -
%    (2)data&pilot :  p - d d d d - -
%    (1)guard      :  - - d d d d - -
%    (0)guard      :  - - d d d d - -
% 
% map plan 2
%   (11)guard      :  - - d d d d - -
%   (10)guard      :  - - d d d d - -
%    (9)guard      :  - - d d d d - -
%    (8)guard      :  - - d d d d - -
%    (7)data       :  d d d d d d d d
%    (6)data       :  d d d d d d d d
%    (5)data       :  d d d d d d d d
%    (4)data       :  d d d d d d d d
%    (3)data&pilot :  - - d d d d - -
%    (2)data&pilot :  - - d d d d - -
%    (1)data&pilot :  - - d d d d - -
%    (0)data&pilot :  p - d d d d - -

function tx_sym_dd_ndft = otfs_sym_map_r1(data_sym_vec, num, map_plan)

% check parameter
if ~ismember(map_plan, [1 2])
    error('Map_plan must be either 1 or 2.')
end

% generate pilot symbols
pilot_sym = zeros(num.num_pilot_delay, num.num_pilot_doppler);
pilot_sym(1, 1) = 1;

% pad guard symbols to pilot symbols
guard_sym = zeros(num.num_guard_delay, num.num_pilot_doppler);
pilot_guard_sym = [guard_sym; pilot_sym]*sqrt((num.num_pilot_delay+num.num_guard_delay)*num.num_pilot_doppler);

% fftshift pilot symbols
pilot_guard_sym_upper = pilot_guard_sym(:, 1:ceil(num.num_pilot_doppler/2));
pilot_guard_sym_lower = pilot_guard_sym(:, ceil(num.num_pilot_doppler/2)+1:end);

% reshape data symbols
data_sym_vec_a = data_sym_vec(1:(num.num_pilot_delay+num.num_guard_delay)*(num.num_ofdmsym_per_subframe-num.num_pilot_doppler));
data_sym_vec_b = data_sym_vec((num.num_pilot_delay+num.num_guard_delay)*(num.num_ofdmsym_per_subframe-num.num_pilot_doppler)+1:end);
data_sym_a = reshape(data_sym_vec_a, num.num_pilot_delay+num.num_guard_delay, []);
data_sym_b = reshape(data_sym_vec_b, num.ndft-(num.num_pilot_delay+num.num_guard_delay), []);

% map resource block
tx_sym_dd_ndft_guard = [ ...
    pilot_guard_sym_upper, data_sym_a, pilot_guard_sym_lower;
    data_sym_b];

% send head guard to tail
if map_plan == 1
    tx_sym_dd_ndft = [ ...
        tx_sym_dd_ndft_guard(ceil(num.num_guard_delay/2)+1:end, :);
        tx_sym_dd_ndft_guard(1:ceil(num.num_guard_delay/2), :)];
else
    tx_sym_dd_ndft = [ ...
        tx_sym_dd_ndft_guard(num.num_guard_delay+1:end, :);
        tx_sym_dd_ndft_guard(1:num.num_guard_delay, :)];
end

% % dump variables
% assignin('base', 'pilot_sym', pilot_sym);
% assignin('base', 'guard_sym', guard_sym);
% assignin('base', 'pilot_guard_sym', pilot_guard_sym);
% assignin('base', 'pilot_guard_sym_up', pilot_guard_sym_up);
% assignin('base', 'pilot_guard_sym_down', pilot_guard_sym_down);
% assignin('base', 'data_sym_a', data_sym_a);
% assignin('base', 'data_sym_b', data_sym_b);
% assignin('base', 'tx_sym_dd_ndft_guard', tx_sym_dd_ndft_guard);
% assignin('base', 'tx_sym_dd_ndft', tx_sym_dd_ndft);
% 
% % plot results
% figure, mesh(1:14, 1:600, abs(tx_sym_dd_ndft))   % dd real channel
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Pilot Plan')
% pause

end
