% demap otfs data and pilot qam symbols

function [rx_sym_dd_ndft_data, rx_sym_dd_ndft_pilot] = otfs_sym_demap_r1(rx_sym_dd_ndft, num, map_plan)

% check parameter
if ~ismember(map_plan, [1 2])
    error('Map_plan must be either 1 or 2.')
end

% send tail guard back to head
if map_plan == 1
    rx_sym_dd_ndft_guard = [ ...
        rx_sym_dd_ndft(end-ceil(num.num_delay_guard/2)+1:end, :);
        rx_sym_dd_ndft(1:end-ceil(num.num_delay_guard/2), :)];
else
    rx_sym_dd_ndft_guard = [ ...
        rx_sym_dd_ndft(end-num.num_delay_guard+1:end, :);
        rx_sym_dd_ndft(1:end-num.num_delay_guard, :)];
end

% demap data symbols
data_sym_a = rx_sym_dd_ndft_guard(1:num.num_delay_pilot+num.num_delay_guard, ceil(num.num_doppler_pilot/2)+1:end-floor(num.num_doppler_pilot/2));
data_sym_b = rx_sym_dd_ndft_guard(num.num_delay_pilot+num.num_delay_guard+1:end, :);

% demap pilot symbols
pilot_sym_upper = rx_sym_dd_ndft_guard(num.num_delay_guard+1:num.num_delay_guard+num.num_delay_pilot, 1:ceil(num.num_doppler_pilot/2));
pilot_sym_lower = rx_sym_dd_ndft_guard(num.num_delay_guard+1:num.num_delay_guard+num.num_delay_pilot, end-floor(num.num_doppler_pilot/2)+1:end);
% pilot_sym_upper = rx_sym_dd_ndft_guard(1:num.num_delay_guard+num.num_delay_pilot, 1:ceil(num.num_doppler_pilot/2));
% pilot_sym_lower = rx_sym_dd_ndft_guard(1:num.num_delay_guard+num.num_delay_pilot, end-floor(num.num_doppler_pilot/2)+1:end);
rx_sym_dd_ndft_pilot = [pilot_sym_upper, pilot_sym_lower];

% serialize data symbols
rx_sym_dd_ndft_data = [data_sym_a(:); data_sym_b(:)];

% % dump variables
% assignin('base', 'rx_sym_dd_ndft', rx_sym_dd_ndft);
% assignin('base', 'rx_sym_dd_ndft_guard', rx_sym_dd_ndft_guard);
% assignin('base', 'data_sym_a', data_sym_a);
% assignin('base', 'data_sym_b', data_sym_b);
% assignin('base', 'rx_sym_dd_ndft_data', rx_sym_dd_ndft_data);
% assignin('base', 'rx_sym_dd_ndft_pilot', rx_sym_dd_ndft_pilot);
% 
% % plot results
% figure, mesh(1:num.num_ofdmsym_per_subframe-num.num_doppler_pilot, 1:num.num_delay_pilot+num.num_delay_guard, abs(data_sym_a))
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Pilot Plan')
% figure, mesh(1:num.num_ofdmsym_per_subframe, 1:num.ndft-(num.num_delay_pilot+num.num_delay_guard), abs(data_sym_b))
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Pilot Plan')
% figure, mesh(1:num.num_doppler_pilot, 1:num.num_delay_pilot, abs(rx_sym_dd_ndft_pilot))
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Pilot Plan')
% pause

end
