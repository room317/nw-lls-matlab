% map otfs data and pilot qam symbols
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

function [tx_sym_rbs_dd, tx_sym_pilot1_usrfrm, tx_sym_pilot2_usrfrm] = otfs_sym_map_r3(tx_sym_data_usrfrm, num, chest_option, test_option)

% set normalized data and pilot symbol power
if strcmp(chest_option, 'dd_impulse')       % use impulse pilot
    % set power for papr reduction
    if test_option.otfs_pilot_pwr_set
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
elseif strcmp(chest_option, 'dd_zc')        % use pilot sequence (zc)
    % set power for sequence spreading
    pwr_data = 1;
    pwr_pilot = num.num_delay_pilot_usr*num.num_doppler_pilot_usr/test_option.zc_seq_len;
elseif strcmp(chest_option, 'dd_random')    % use pilot sequence (random)
    % set power for sequence spreading
    pwr_data = 1;
    pwr_pilot = num.num_delay_pilot_usr*num.num_doppler_pilot_usr/test_option.rand_seq_len;
elseif strcmp(chest_option, 'dd_golay_serial') || ...
        strcmp(chest_option, 'dd_golay_parallel') || ...
        strcmp(chest_option, 'dd_golay_diag')
    % set power for sequence spreading
    pwr_data = 1;
    pwr_pilot = num.num_delay_pilot_usr*num.num_doppler_pilot_usr/test_option.golay_seq_len;
else
    error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
end

% set index
idx_delay_usr = floor(num.num_delay_usr/2)+1;
idx_doppler_usr = floor(num.num_doppler_usr/2)+1;
idx_delay_pilot_usr = num.num_delay_pilot_half_usr+1;       % pilot: even
idx_doppler_pilot_usr = num.num_doppler_pilot_half_usr+1;   % pilot: even

% set pilot resource block
if strcmp(chest_option, 'dd_impulse')       % use impulse pilot
    % generate pilot symbols
    tx_sym_pilot1_usrfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_usr);
    tx_sym_pilot1_usrfrm(idx_delay_pilot_usr, idx_doppler_pilot_usr) = sqrt(pwr_pilot);              % boosts pilots
    tx_sym_pilot2_usrfrm = [];
elseif strcmp(chest_option, 'dd_zc')        % use sequence pilot
    % check numbers
    if num.num_delay_pilot_usr < test_option.zc_seq_len+num.num_delay_guard_a_usr+num.num_delay_guard_b_usr
        fprintf('num_delay_usr = %d\n', num.num_delay_usr)
        fprintf('num_delay_pilot_usr = %d\n', num.num_delay_pilot_usr)
        fprintf('num_delay_guard_a_usr = %d\n', num.num_delay_guard_a_usr)
        fprintf('num_delay_guard_b_usr = %d\n', num.num_delay_guard_b_usr)
        fprintf('test_option.zc_seq_len = %d\n', test_option.zc_seq_len)
        error('Check numbers!')
    end
    
    % map spread sequence
    tx_sym_pilot1_usrfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_usr);
    tx_sym_pilot1_usrfrm(idx_delay_pilot_usr-ceil(test_option.zc_seq_len/2)+1: ...
        idx_delay_pilot_usr+floor(test_option.zc_seq_len/2), ...
        idx_doppler_pilot_usr) = ...
        sqrt(pwr_pilot)*test_option.zc_seq;
    tx_sym_pilot2_usrfrm = [];
elseif strcmp(chest_option, 'dd_random')   % use random pilot
    % check numbers
    if num.num_delay_pilot_usr < test_option.rand_seq_len+num.num_delay_guard_a_usr+num.num_delay_guard_b_usr
        fprintf('num_delay_usr = %d\n', num.num_delay_usr)
        fprintf('num_delay_pilot_usr = %d\n', num.num_delay_pilot_usr)
        fprintf('num_delay_guard_a_usr = %d\n', num.num_delay_guard_a_usr)
        fprintf('num_delay_guard_b_usr = %d\n', num.num_delay_guard_b_usr)
        fprintf('test_option.zc_seq_len = %d\n', test_option.rand_seq_len)
        error('Check numbers!')
    end
    
    % map spread sequence
    tx_sym_pilot1_usrfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_usr);
    tx_sym_pilot1_usrfrm(idx_delay_pilot_usr-ceil(test_option.rand_seq_len/2)+1: ...
        idx_delay_pilot_usr+floor(test_option.rand_seq_len/2), ...
        idx_doppler_pilot_usr) = ...
        sqrt(pwr_pilot)*test_option.rand_seq;
    tx_sym_pilot2_usrfrm = [];
elseif strcmp(chest_option, 'dd_golay_serial')      % use complementary sequence pilot
    % check numbers
    if num.num_delay_pilot_usr < 2*test_option.golay_seq_len+2*num.num_delay_guard_a_usr+2*num.num_delay_guard_b_usr
        fprintf('num_delay_usr = %d\n', num.num_delay_usr)
        fprintf('num_delay_pilot_usr = %d\n', num.num_delay_pilot_usr)
        fprintf('num_delay_guard_a_usr = %d\n', num.num_delay_guard_a_usr)
        fprintf('num_delay_guard_b_usr = %d\n', num.num_delay_guard_b_usr)
        fprintf('test_option.golay_seq_len = %d\n', test_option.golay_seq_len)
        error('Check numbers!')
    end
    
    % map spread sequence (pilot: even number)
    tx_sym_pilot1_usrfrm = zeros(num.num_delay_pilot_half_usr, num.num_doppler_pilot_usr);
    tx_sym_pilot2_usrfrm = zeros(num.num_delay_pilot_half_usr, num.num_doppler_pilot_usr);
    if num.num_doppler_usr == num.num_doppler_pilot_usr             % full doppler spreading
        tx_sym_pilot1_usrfrm((floor(num.num_delay_pilot_half_usr/2)+1)-ceil(test_option.golay_seq_len/2)+1: ...
            (floor(num.num_delay_pilot_half_usr/2)+1)+floor(test_option.golay_seq_len/2), ...
            floor(num.num_doppler_usr/4)+1) = ...
            sqrt(pwr_pilot/2)*test_option.golay_seq_a;
        tx_sym_pilot2_usrfrm((floor(num.num_delay_pilot_half_usr/2)+1)-ceil(test_option.golay_seq_len/2)+1: ...
            (floor(num.num_delay_pilot_half_usr/2)+1)+floor(test_option.golay_seq_len/2), ...
            floor(num.num_doppler_usr/4)+idx_doppler_usr) = ...
            sqrt(pwr_pilot/2)*test_option.golay_seq_b;
    else                                                            % partial doppler spreading
        tx_sym_pilot1_usrfrm(floor(num.num_delay_pilot_half_usr/2)-ceil(test_option.golay_seq_len/2)+1: ...
            floor(num.num_delay_pilot_half_usr/2)+floor(test_option.golay_seq_len/2), ...
            idx_doppler_pilot_usr) = ...
            sqrt(pwr_pilot/2)*test_option.golay_seq_a;
        tx_sym_pilot2_usrfrm(floor(num.num_delay_pilot_half_usr/2)-ceil(test_option.golay_seq_len/2)+1: ...
            floor(num.num_delay_pilot_half_usr/2)+floor(test_option.golay_seq_len/2), ...
            idx_doppler_pilot_usr) = ...
            sqrt(pwr_pilot/2)*test_option.golay_seq_b;
    end
elseif strcmp(chest_option, 'dd_golay_parallel')    % use complementary sequence pilot
    % check numbers
    if num.num_doppler_pilot_usr < 2+3*num.num_doppler_guard_usr
        fprintf('num_doppler_usr = %d\n', num.num_doppler_usr)
        fprintf('num_doppler_pilot_usr = %d\n', num.num_doppler_pilot_usr)
        fprintf('num_doppler_guard_usr = %d\n', num.num_doppler_guard_usr)
        error('Check numbers!')
    end
    if num.num_delay_pilot_usr < test_option.golay_seq_len+num.num_delay_guard_a_usr+num.num_delay_guard_b_usr
        fprintf('num_delay_usr = %d\n', num.num_delay_usr)
        fprintf('num_delay_pilot_usr = %d\n', num.num_delay_pilot_usr)
        fprintf('num_delay_guard_a_usr = %d\n', num.num_delay_guard_a_usr)
        fprintf('num_delay_guard_b_usr = %d\n', num.num_delay_guard_b_usr)
        fprintf('test_option.zc_seq_len = %d\n', test_option.golay_seq_len)
        error('Check numbers!')
    end
    
    % map spread sequence (pilot: even number)
    tx_sym_pilot1_usrfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_half_usr);
    tx_sym_pilot2_usrfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_half_usr);
    tx_sym_pilot1_usrfrm(idx_delay_pilot_usr-ceil(test_option.golay_seq_len/2)+1: ...
        idx_delay_pilot_usr+floor(test_option.golay_seq_len/2), ...
        floor(num.num_doppler_pilot_half_usr/2)+1) = ...
        sqrt(pwr_pilot/2)*test_option.golay_seq_a;
    tx_sym_pilot2_usrfrm(idx_delay_pilot_usr-ceil(test_option.golay_seq_len/2)+1: ...
        idx_delay_pilot_usr+floor(test_option.golay_seq_len/2), ...
        floor(num.num_doppler_pilot_half_usr/2)+1) = ...
        sqrt(pwr_pilot/2)*test_option.golay_seq_b;
elseif strcmp(chest_option, 'dd_golay_diag')        % use complementary sequence pilot
    % check numbers
    if (num.num_doppler_pilot_usr ~= num.num_doppler_usr) || ...
            (num.num_delay_usr < 2*num.num_delay_pilot_usr) || ...
            (num.num_delay_pilot_usr < 2*test_option.golay_seq_len+2*num.num_delay_guard_a_usr+2*num.num_delay_guard_b_usr)
        fprintf('num_delay_usr = %d\n', num.num_delay_usr)
        fprintf('num_doppler_usr = %d\n', num.num_doppler_usr)
        fprintf('num_delay_pilot_usr = %d\n', num.num_delay_pilot_usr)
        fprintf('num_doppler_pilot_usr = %d\n', num.num_doppler_pilot_usr)
        fprintf('num_delay_guard_a_usr = %d\n', num.num_delay_guard_a_usr)
        fprintf('num_delay_guard_b_usr = %d\n', num.num_delay_guard_b_usr)
        fprintf('num_doppler_guard_usr = %d\n', num.num_doppler_guard_usr)
        error('Check numbers!')
    end
    
    % map spread sequence
    tx_sym_pilot1_usrfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_half_usr);
    tx_sym_pilot2_usrfrm = zeros(num.num_delay_pilot_usr, num.num_doppler_pilot_half_usr);
    tx_sym_pilot1_usrfrm(idx_delay_pilot_usr-ceil(test_option.golay_seq_len/2)+1: ...
        idx_delay_pilot_usr+floor(test_option.golay_seq_len/2), ...
        floor(num.num_doppler_pilot_half_usr/2)+1) = ...
        sqrt(pwr_pilot/2)*test_option.golay_seq_a;
    tx_sym_pilot2_usrfrm(idx_delay_pilot_usr-ceil(test_option.golay_seq_len/2)+1: ...
        idx_delay_pilot_usr+floor(test_option.golay_seq_len/2), ...
        floor(num.num_doppler_pilot_half_usr/2)+1) = ...
        sqrt(pwr_pilot/2)*test_option.golay_seq_b;
else
    error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
end

% map data and pilot symbols to resource block
if strcmp(chest_option, 'dd_impulse') || ...
        strcmp(chest_option, 'dd_zc') || strcmp(chest_option, 'dd_random')
    % reshape data symbols
    tx_sym_data1_usrfrm = sqrt(pwr_data)* ...
        reshape(tx_sym_data_usrfrm(1:num.num_delay_pilot_usr*num.num_doppler_data_usr), ...
        num.num_delay_pilot_usr, []);
    tx_sym_data2_usrfrm = sqrt(pwr_data)* ...
        reshape(tx_sym_data_usrfrm(num.num_delay_pilot_usr*num.num_doppler_data_usr+1:end), ...
        num.num_delay_data_usr, []);
    
    % map resource block
    tx_sym_map_base = [ ...
        tx_sym_pilot1_usrfrm, tx_sym_data1_usrfrm;
        tx_sym_data2_usrfrm];
elseif strcmp(chest_option, 'dd_golay_serial')
    % reshape data symbols
    tx_sym_data1_usrfrm = sqrt(pwr_data)* ...
        reshape(tx_sym_data_usrfrm(1:num.num_delay_pilot_usr*num.num_doppler_data_usr), ...
        num.num_delay_pilot_usr, []);
    tx_sym_data2_usrfrm = sqrt(pwr_data)* ...
        reshape(tx_sym_data_usrfrm(num.num_delay_pilot_usr*num.num_doppler_data_usr+1:end), ...
        num.num_delay_data_usr, []);
    
    % map resource block
    tx_sym_map_base = [ ...
        [tx_sym_pilot1_usrfrm; tx_sym_pilot2_usrfrm], tx_sym_data1_usrfrm;
        tx_sym_data2_usrfrm];
elseif strcmp(chest_option, 'dd_golay_parallel')
    % reshape data symbols
    tx_sym_data1_usrfrm = sqrt(pwr_data)* ...
        reshape(tx_sym_data_usrfrm(1:num.num_delay_pilot_usr*num.num_doppler_data_usr), ...
        num.num_delay_pilot_usr, []);
    tx_sym_data2_usrfrm = sqrt(pwr_data)* ...
        reshape(tx_sym_data_usrfrm(num.num_delay_pilot_usr*num.num_doppler_data_usr+1:end), ...
        num.num_delay_data_usr, []);
    
    % map resource block
    tx_sym_map_base = [ ...
        [tx_sym_pilot1_usrfrm, tx_sym_pilot2_usrfrm], tx_sym_data1_usrfrm;
        tx_sym_data2_usrfrm];
elseif strcmp(chest_option, 'dd_golay_diag')        % use golay sequence
    % reshape data symbols
    tx_sym_data1_usrfrm = sqrt(pwr_data)* ...
        reshape(tx_sym_data_usrfrm(1:num.num_delay_pilot_usr*num.num_doppler_pilot_half_usr), ...
        num.num_delay_pilot_usr, []);
    tx_sym_data2_usrfrm = sqrt(pwr_data)* ...
        reshape(tx_sym_data_usrfrm(num.num_delay_pilot_usr*num.num_doppler_pilot_half_usr+1:num.num_delay_pilot_usr*num.num_doppler_pilot_usr), ...
        num.num_delay_pilot_usr, []);
    tx_sym_data3_usrfrm = sqrt(pwr_data)* ...
        reshape(tx_sym_data_usrfrm(num.num_delay_pilot_usr*num.num_doppler_pilot_usr+1:end), ...
        [], num.num_doppler_usr);
    
    % map resource block
    tx_sym_map_base = [ ...
        tx_sym_pilot1_usrfrm, tx_sym_data1_usrfrm;
        tx_sym_data2_usrfrm, tx_sym_pilot2_usrfrm;
        tx_sym_data3_usrfrm];
else
    error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
end

% relocate symbols according to map plans
if strcmp(chest_option, 'dd_impulse') || ...
        strcmp(chest_option, 'dd_zc') || strcmp(chest_option, 'dd_random') || ...
        strcmp(chest_option, 'dd_golay_serial') || strcmp(chest_option, 'dd_golay_parallel')
    map_shift = [idx_delay_usr-idx_delay_pilot_usr, idx_doppler_usr-idx_doppler_pilot_usr];
elseif strcmp(chest_option, 'dd_golay_diag')        % use golay sequence
    map_shift = [idx_delay_usr-num.num_delay_pilot_usr, idx_doppler_usr-idx_doppler_pilot_usr];
else
    error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
end
tx_sym_rbs_dd = circshift(tx_sym_map_base, map_shift);

% % test
% tone_pos1 = [idx_delay_pilot_usr, idx_doppler_pilot_usr];
% tone_pos2 = tone_pos1+map_shift;
% fprintf('total tx avg pwr:%6.3f\n', sqrt(mean(abs(tx_sym_map_base).^2, 'all')))
% if strcmp(chest_option, 'dd_impulse') || ...
%         strcmp(chest_option, 'dd_zc') || ...
%         strcmp(chest_option, 'dd_random')
%     fprintf('tx pilot avg pwr:%6.3f\n', sqrt(mean(abs(tx_sym_pilot1_usrfrm).^2, 'all')));
% elseif strcmp(chest_option, 'dd_golay_serial') || ...
%         strcmp(chest_option, 'dd_golay_parallel') || ...
%         strcmp(chest_option, 'dd_golay_diag')
%     fprintf('tx pilot1 avg pwr:%6.3f\n', sqrt(mean(abs(tx_sym_pilot1_usrfrm).^2, 'all')));
%     fprintf('tx pilot2 avg pwr:%6.3f\n', sqrt(mean(abs(tx_sym_pilot2_usrfrm).^2, 'all')));
% else
%     error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
% end
% figure, mesh(1:num.num_doppler_usr, 1:num.num_delay_usr, abs(tx_sym_rbs_dd))   % dd real channel
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('TX Pilot Plan')
% fprintf('tone position in pilot plain: (%d, %d)\n', tone_pos1)
% fprintf('tone position in resource plain: (%d, %d)\n', tone_pos2)
% 
% % dump variables
% assignin('base', 'tx_sym_data_usrfrm', tx_sym_data_usrfrm);
% assignin('base', 'tx_sym_data1_usrfrm', tx_sym_data1_usrfrm);
% assignin('base', 'tx_sym_data2_usrfrm', tx_sym_data2_usrfrm);
% assignin('base', 'tx_sym_map_basic', tx_sym_map_base);
% assignin('base', 'tx_sym_rbs_dd', tx_sym_rbs_dd);
% if strcmp(chest_option, 'dd_impulse') || ...
%         strcmp(chest_option, 'dd_zc') || ...
%         strcmp(chest_option, 'dd_random')
%     assignin('base', 'tx_sym_pilot1_usrfrm', tx_sym_pilot1_usrfrm);
% elseif strcmp(chest_option, 'dd_golay_serial') || ...
%         strcmp(chest_option, 'dd_golay_parallel')
%     assignin('base', 'tx_sym_pilot1_usrfrm', tx_sym_pilot1_usrfrm);
%     assignin('base', 'tx_sym_pilot2_usrfrm', tx_sym_pilot2_usrfrm);
% elseif strcmp(chest_option, 'dd_golay_diag')
%     assignin('base', 'tx_sym_pilot1_usrfrm', tx_sym_pilot1_usrfrm);
%     assignin('base', 'tx_sym_pilot2_usrfrm', tx_sym_pilot2_usrfrm);
%     assignin('base', 'tx_sym_data3_usrfrm', tx_sym_data3_usrfrm);
% else
%     error('dd-domain ''chest_option'' must be one of these: {dd_impulse, dd_zc, dd_golay_serial, dd_golay_parallel, dd_golay_diag, dd_random}')
% end
% pause

end
