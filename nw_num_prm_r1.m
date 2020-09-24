function nw_num = nw_num_prm_r1(scs_khz, bw_mhz, num_slot, waveform, chest_option, usr_option)

% timing parameter reference
% 1 frame = 10 ms
% 1 subframe = 1 ms
% 1 slot = 14 ofdm symbols
% 1 mini-slot = 7/4/2 ofdm symbols

% terminology
% re: resource element
% rb (lte): resource block, 12 consecutive subcarriers in 1 slot
% rb (nr): resource block, 12 consecutive subcarriers

% set subcarrier spacing parameters
scs_khz_list = [15 30 60];
idx_scs = find(scs_khz_list == scs_khz, 1);
if isempty(idx_scs)
    error('Subcarrier spacing must be one of these: {15, 30, 60}')
end

% set bandwidth parameters
bw_mhz_list = [5 10 15 20 25 30 40 50 60 80 90 100];
idx_bw = find(bw_mhz_list == bw_mhz, 1);
if isempty(idx_bw)
    error('Bandwidth must be one of these: {5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 90, 100}')
end

% set resource block table
% ref.: http://howltestuffworks.blogspot.com/2019/11/5g-nr-resource-blocks.html
% ref.: https://www.sharetechnote.com/html/5G/5G_FR_Bandwidth.html
%     bw (mhz)     |   5 |  10 |  15 |  20 |  25 |  30 |  40 |  50 |  60 |  80 |  90 | 100
%     -------------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+----
%     scs 15 (khz) |  25 |  52 |  79 | 106 | 133 | 160 | 216 | 270 |   0 |   0 |   0 |   0
%     scs 30 (khz) |  11 |  24 |  38 |  51 |  65 |  78 | 106 | 133 | 162 | 217 | 245 | 273
%     scs 60 (khz) |   0 |  11 |  18 |  24 |  31 |  38 |  51 |  65 |  79 | 107 | 121 | 135
rb_tbl = [ ...
    25   52   79  106  133  160  216  270    0    0    0    0
    11   24   38   51   65   78  106  133  162  217  245  273
     0   11   18   24   31   38   51   65   79  107  121  135];

% set numerical parameters
num_rb = rb_tbl(idx_scs, idx_bw);           % number of rbs
num_subc_rb = 12;                           % number of subcarriers per rb
num_subc_bw = num_rb*num_subc_rb;           % number of subcarriers in bandwidth
num_fft = 2^ceil(log2(num_subc_bw));        % fft size
sample_rate = scs_khz*1e3*num_fft;          % sampling rate (hz)
num_cp = round(num_fft*144/2048);           % number of samples in cyclic prefix
num_ofdmsym_slot = 14;                      % number of ofdm symbols per slot
num_ofdmsym = num_ofdmsym_slot*num_slot;    % number of ofdm symbols

% set user parameters (full span, prb_size: 12*14)
if strcmp(usr_option, 'su')                 % full spreading
    num_rb_usr = num_rb;                    % number of resource blocks per user
    num_slot_usr = num_slot;                % number of slots per user
    list_usr = [];                          % user index list
elseif strcmp(usr_option, 'mu')
    num_rb_usr = 3;                         % number of resource blocks per user
    num_slot_usr = 1;                       % number of slots per user
    list_usr = [1 2 3 4];                   % user index list
else
    error('''test_usr'' shall be either ''su'' or ''mu''.\n')
end
num_subc_usr = num_rb_usr*num_subc_rb;              % number of subcarriers per user
num_ofdmsym_usr = num_ofdmsym_slot*num_slot_usr;    % number of ofdm symbols per user (slot-based)
max_usr_rb = floor(num_rb/num_rb_usr);              % max. user rb index
max_usr_slot = floor(num_slot/num_slot_usr);        % max. user slot index
max_usr = max_usr_rb*max_usr_slot;                  % max. number of users
num_usr = length(list_usr);                         % number of users
if sum(double(list_usr > max_usr)) > 0
    error('max. user index shall be %d.\n', max_usr)
end

% set index parameters
if strcmp(waveform, 'ofdm')    % ofdm
    
    % set pilot index parameters
    if strcmp(chest_option, 'tf_lteup')                                 % tf-domain pilots (some symbols for pilots)
%         idx_ofdmsym_pilot_usr = repmat([4 11], 1, num_slot);            % pilot symbol index (should be scalar or vector)
        idx_ofdmsym_pilot_usr = reshape([4; 11]+(0:num_slot_usr-1)*num_ofdmsym_slot, 1, []);    % pilot symbol index (should be scalar or vector)
        idx_subc_pilot_usr = repmat((1:num_subc_usr)', 1, length(idx_ofdmsym_pilot_usr));
    elseif strcmp(chest_option, 'tf_ltedown')                           % tf-domain pilots (some symbols for pilots)
%         idx_ofdmsym_pilot_usr = repmat([1 5 8 12], 1, num_slot);        % pilot symbol index (should be scalar or vector)
        idx_ofdmsym_pilot_usr = reshape([1; 5; 8; 12]+(0:num_slot_usr-1)*num_ofdmsym_slot, 1, []);    % pilot symbol index (should be scalar or vector)
        idx_subc_pilot_usr = repmat([(6:6:num_subc_usr)', (3:6:num_subc_usr)', (6:6:num_subc_usr)', (3:6:num_subc_usr)'], 1, num_slot_usr);
    elseif strcmp(chest_option, 'tf_nr')                                % tf-domain pilots (some symbols for pilots)
%         idx_ofdmsym_pilot_usr = repmat([3 6 9 12], 1, num_slot);        % pilot symbol index (should be scalar or vector)
        idx_ofdmsym_pilot_usr = reshape([3; 6; 9; 12]+(0:num_slot_usr-1)*num_ofdmsym_slot, 1, []);    % pilot symbol index (should be scalar or vector)
        idx_subc_pilot_usr = repmat([(2:2:num_subc_usr)', (2:2:num_subc_usr)', (2:2:num_subc_usr)', (2:2:num_subc_usr)'], 1, num_slot_usr);
    else    % whole resources for data (real, perfect)
        idx_ofdmsym_pilot_usr = [];
        idx_subc_pilot_usr = [];
    end
    num_ofdmsym_pilot_usr = length(idx_ofdmsym_pilot_usr);      % number of ofdm symbols with pilots
    num_subc_pilot_usr = size(idx_subc_pilot_usr, 1);           % number of pilot subcarriers in 1 ofdm symbol with pilots
    
    % set data index parameter
    idx_ofdmsym_data_usr = setxor(1:num_ofdmsym_usr, idx_ofdmsym_pilot_usr);
    idx_subc_data_usr = zeros(num_subc_usr-num_subc_pilot_usr, num_ofdmsym_pilot_usr);
    for i = 1:num_ofdmsym_pilot_usr
        tmp_idx_subc_data_usr = setxor((1:num_subc_usr)', idx_subc_pilot_usr(:, i));
        idx_subc_data_usr(:, i) = tmp_idx_subc_data_usr;
    end
    num_ofdmsym_data_usr = length(idx_ofdmsym_data_usr);        % number of ofdm symbols without pilots (all data subcarrier)
    num_subc_data_usr = num_subc_usr-num_subc_pilot_usr;        % number of data subcarriers in 1 ofdm symbol with pilots
    
    % calculate number of data qam symbols per user
    num_qamsym_usr = (num_subc_usr*num_ofdmsym_usr)-(num_subc_pilot_usr*num_ofdmsym_pilot_usr);
    
    % set unused parameters
    num_doppler_usr = [];
    num_delay_usr = [];
    num_doppler_pilot_usr = [];
    num_delay_pilot_usr = [];
    num_doppler_guard_usr = [];
    num_delay_guard_usr = [];
    num_doppler_data_usr = [];
    num_delay_data_usr = [];
else
    
    % set pilot resource ratio (otfs pilot setup)
    %   - num. otfs pilot(including guard) qam symbols = num. ofdm pilot qam symbols
    %     ex) 600*2 for ofdm (14.29%) and 86*14 for otfs (14.33%)
    delay_pilot_ratio = 0.4;        % num_delay_pilot/num_delay_data (40% of available delay grids, even number)
    doppler_pilot_ratio = 1.0;      % num_doppler_pilot/num_doppler_data (100% of available doppler grids)
    delay_guard_ratio = 0.2;        % num_delay_guard/num_delay_pilot (20% of pilot delay grids)
    doppler_guard_ratio = 0.2;      % num_doppler_guard/num_doppler_pilot (20% of pilot doppler grids)
    
    % user parameters for otfs
    num_doppler_usr = num_ofdmsym_usr;              % corresponds to ofdm symbols (fixed)
    num_delay_usr = num_subc_usr;                   % corresponds to subcarriers (fixed)
    num_doppler_pilot_usr = double(strcmp(chest_option, 'dd_tone'))*round(num_doppler_usr*(doppler_pilot_ratio/2))*2;           % number of doppler grids with pilots
    num_delay_pilot_usr = double(strcmp(chest_option, 'dd_tone'))*round(num_delay_usr*(delay_pilot_ratio/2))*2;                 % number of delay grids with pilots
    num_doppler_guard_usr = double(num_doppler_usr~=num_doppler_pilot_usr)*round(num_doppler_pilot_usr*doppler_guard_ratio);    % number of doppler grids for guard (set 0 for full-span pilot)
    num_delay_guard_usr = double(num_delay_usr~=num_delay_pilot_usr)*round(num_delay_pilot_usr*delay_guard_ratio);              % number of delay grids for guard (set even number)
    num_doppler_data_usr = num_doppler_usr-num_doppler_pilot_usr;   % number of doppler grids without pilot and guard (set even number)
    num_delay_data_usr = num_delay_usr-num_delay_pilot_usr;         % number of delay grids without pilot and guard (set even number)
    
    % calculate number of data qam symbols per user
    num_qamsym_usr = (num_delay_usr*num_doppler_usr)-(num_delay_pilot_usr*num_doppler_pilot_usr);
    
    % set unused parameters
    idx_ofdmsym_pilot_usr = [];
    idx_subc_pilot_usr = [];
    idx_ofdmsym_data_usr = [];
    idx_subc_data_usr = [];
    num_ofdmsym_pilot_usr = [];
    num_subc_pilot_usr = [];
    num_ofdmsym_data_usr = [];
    num_subc_data_usr = [];
end

% set timing parameters
%   - len_frame = len_rb_sym_user * (nfft + num_cp);
%   - t_usrfrm = (length of 1 symbol + CP) * number of symbol * length of symbol
t_usrfrm = (num_fft+num_cp)*num_ofdmsym/sample_rate;   % transmission time interval (sec)

% output
% nw_num.nfft = nfft;
% nw_num.ndft = num_subc_usr;
% nw_num.len_rb_sym_user = num_qamsym_usr;
% nw_num.num_subcarrier = num_subc_subfrm;
% nw_num.sample_rate = sample_rate;
% nw_num.num_ofdmsym_per_subframe = num_ofdmsym_usr;
% nw_num.num_data_ofdmsym_per_subframe = num_ofdmsym_data_usr;
% nw_num.num_pilot_ofdmsym_per_subframe = num_ofdmsym_pilot_usr;
% nw_num.idx_data_ofdmsym = idx_ofdmsym_data_usr;
% nw_num.idx_pilot_ofdmsym = idx_ofdmsym_pilot_usr;
% nw_num.idx_pilot_subc = idx_subc_pilot_usr;
% nw_num.num_cp = num_cp;
% nw_num.t_subframe = t_subfrm;
% nw_num.num_pilot_doppler = num_doppler_pilot;
% nw_num.num_pilot_delay = num_delay_pilot;
% nw_num.num_guard_delay = num_delay_guard;

nw_num.num_fft = num_fft;
nw_num.sample_rate = sample_rate;
nw_num.num_cp = num_cp;
nw_num.num_subc_rb = num_subc_rb;
nw_num.num_subc_bw = num_subc_bw;                   % updated
nw_num.num_ofdmsym_slot = num_ofdmsym_slot;
nw_num.num_ofdmsym = num_ofdmsym;
nw_num.num_qamsym_usr = num_qamsym_usr;
nw_num.t_usrfrm = t_usrfrm;

nw_num.num_rb_usr = num_rb_usr;
nw_num.num_slot_usr = num_slot_usr;
nw_num.list_usr = list_usr;
nw_num.max_usr_rb = max_usr_rb;
nw_num.max_usr_slot = max_usr_slot;
nw_num.num_usr = num_usr;

nw_num.num_subc_usr = num_subc_usr;
nw_num.num_ofdmsym_usr = num_ofdmsym_usr;
nw_num.idx_subc_pilot_usr = idx_subc_pilot_usr;
nw_num.idx_ofdmsym_pilot_usr = idx_ofdmsym_pilot_usr;
nw_num.idx_subc_data_usr = idx_subc_data_usr;
nw_num.idx_ofdmsym_data_usr = idx_ofdmsym_data_usr;
nw_num.num_subc_pilot_usr = num_subc_pilot_usr;
nw_num.num_ofdmsym_pilot_usr = num_ofdmsym_pilot_usr;
nw_num.num_subc_data_usr = num_subc_data_usr;
nw_num.num_ofdmsym_data_usr = num_ofdmsym_data_usr;

nw_num.num_doppler_usr = num_doppler_usr;
nw_num.num_delay_usr = num_delay_usr;
nw_num.num_doppler_pilot_usr = num_doppler_pilot_usr;
nw_num.num_delay_pilot_usr = num_delay_pilot_usr;
nw_num.num_doppler_guard_usr = num_doppler_guard_usr;
nw_num.num_delay_guard_usr = num_delay_guard_usr;
nw_num.num_doppler_data_usr = num_doppler_data_usr;
nw_num.num_delay_data_usr = num_delay_data_usr;

end
