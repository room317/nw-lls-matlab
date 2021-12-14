% gen_ch_clip generates clipped ideal(real) channel matrix
%     ch = gen_ch_clip(ch, num)
%     ch: (M x N) channel matrix (convolutional channel matrix)
%     num: numerical parameters

function ch_conv_clip = gen_ch_clip(ch_conv_real, num)

% set clipping area
% delay_pilot_ratio = 0.24;        % num_delay_pilot/num_delay_data (40% of available delay grids, even number)
% doppler_pilot_ratio = 1.0;      % num_doppler_pilot/num_doppler_data (100% of available doppler grids)
% delay_guard_a_ratio = 0.44;        % num_delay_guard/num_delay_pilot (20% of pilot delay grids)
% delay_guard_b_ratio = 0.06;        % num_delay_guard/num_delay_pilot (20% of pilot delay grids)
% doppler_guard_ratio = 0.0;      % num_doppler_guard/num_doppler_pilot (20% of pilot doppler grids)

delay_pilot_ratio = 0.4;        % num_delay_pilot/num_delay_data (40% of available delay grids, even number)
doppler_pilot_ratio = 1.0;      % num_doppler_pilot/num_doppler_data (100% of available doppler grids)
delay_guard_a_ratio = 0.45;        % num_delay_guard/num_delay_pilot (20% of pilot delay grids)
delay_guard_b_ratio = 0.05;        % num_delay_guard/num_delay_pilot (20% of pilot delay grids)
doppler_guard_ratio = 0.0;      % num_doppler_guard/num_doppler_pilot (20% of pilot doppler grids)

% user parameters for otfs
num_doppler_usr = num.num_doppler_usr;
num_delay_usr = num.num_delay_usr;
num_doppler_pilot_half_usr = round(num_doppler_usr*(doppler_pilot_ratio/2));    % number of half doppler grids with pilots (for complementary sequence pilot)
num_delay_pilot_half_usr = round(num_delay_usr*(delay_pilot_ratio/2));          % number of half delay grids with pilots (for complementary sequence pilot)
num_doppler_pilot_usr = num_doppler_pilot_half_usr*2;           % number of doppler grids with pilots
num_delay_pilot_usr = num_delay_pilot_half_usr*2;               % number of delay grids with pilots
num_doppler_guard_usr = double(num_doppler_usr~=num_doppler_pilot_usr)*round(num_doppler_pilot_usr*doppler_guard_ratio);    % number of doppler grids for guard (set 0 for full-span pilot)
num_delay_guard_a_usr = double(num_delay_usr~=num_delay_pilot_usr)*round(num_delay_pilot_usr*delay_guard_a_ratio);              % number of delay grids for guard (set even number)
num_delay_guard_b_usr = double(num_delay_usr~=num_delay_pilot_usr)*round(num_delay_pilot_usr*delay_guard_b_ratio);              % number of delay grids for guard (set even number)

% set index
idx_list_delay = ...
    floor(num_delay_usr/2)-num_delay_pilot_half_usr+num_delay_guard_a_usr+1: ...
    floor(num_delay_usr/2)-num_delay_pilot_half_usr+num_delay_pilot_usr-num_delay_guard_b_usr;

idx_list_doppler = ...
    floor(num_doppler_usr/2)-num_doppler_pilot_half_usr+num_doppler_guard_usr+1: ...
    floor(num_doppler_usr/2)-num_doppler_pilot_half_usr+num_doppler_pilot_usr-num_doppler_guard_usr;

% extract pilot area (clipping)
ch_conv_real_a = fftshift(fftshift(ch_conv_real, 1), 2);
ch_conv_clip_a = zeros(size(ch_conv_real_a));
ch_conv_clip_a(idx_list_delay, idx_list_doppler) = ch_conv_real_a(idx_list_delay, idx_list_doppler);
ch_conv_clip = fftshift(fftshift(ch_conv_clip_a, 1), 2);

end
