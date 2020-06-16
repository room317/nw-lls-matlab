function nw_num = nw_num_prm(idx_bandwidth, chest_option)

% idx_bandwidth = 4;

% numerical setup
% resource block table (bandwidth_mhz  num_rb)
% bandwidth_mhz |   1.4 |   3 |   5 |   10 |   15 |  20
% --------------+-------+-----+-----+------+------+------
% num_rb        |   6   |  15 |  25 |   50 |   75 |  100
% --------------+-------+-----+-----+------+------+------
% nfft          | 128   | 256 | 512 | 1024 | 1536 | 2048

% set basic parameters
rb_tbl = [6  128;  15  256;  25  512;  50 1024;  75 1536; 100 2048];
num_rb = rb_tbl(idx_bandwidth, 1);
nfft = rb_tbl(idx_bandwidth, 2);
num_subc_rb = 12;                           % number of subcarriers per rb
num_ofdmsym_slot = 7;                         % number of ofdm symbols per slot
sample_rate = 15e3*nfft;                    % sampling rate (hz)
num_cp = nfft * 144 / 2048;                 % number of samples in cyclic prefix

% set subframe parameters
num_subc_bw = num_rb*num_subc_rb;           % number of subcarriers in bandwidth
num_ofdmsym_subfrm = 2*num_ofdmsym_slot;    % number of ofdm symbols per subframe (pilot + data)

% set user parameters (full span)
num_prb_usr = [num_rb 2];                           % number(2d) of physical resource blocks per user (origianlly 3*2 rbs)
num_subc_usr = num_prb_usr(1)*num_subc_rb;          % number of subcarriers per user
num_ofdmsym_usr = num_prb_usr(2)*num_ofdmsym_slot;    % number of ofdm symbols per user (pilot + data)

% set pilot index parameters
if num_ofdmsym_usr ~= 14
    error('There should be 14 ofdm symbols in a user block.')
elseif strcmp(chest_option, 'tf_lteup')         % tf-domain pilots (some symbols for pilots)
    idx_ofdmsym_pilot_usr = [4 11];             % pilot symbol index (should be scalar or vector, all elements should be between 1 and 14)
    idx_subc_pilot_usr = repmat((1:num_subc_usr)', 1, length(idx_ofdmsym_pilot_usr));
elseif strcmp(chest_option, 'tf_ltedown')       % tf-domain pilots (some symbols for pilots)
    idx_ofdmsym_pilot_usr = [1 5 8 12];         % pilot symbol index (should be scalar or vector, all elements should be between 1 and 14)
    idx_subc_pilot_usr = [(6:6:num_subc_usr)', (3:6:num_subc_usr)', (6:6:num_subc_usr)', (3:6:num_subc_usr)'];
elseif strcmp(chest_option, 'tf_nr')            % tf-domain pilots (some symbols for pilots)
    idx_ofdmsym_pilot_usr = [3 6 9 12];         % pilot symbol index (should be scalar or vector, all elements should be between 1 and 14)
    idx_subc_pilot_usr = [(2:2:num_subc_usr)', (2:2:num_subc_usr)', (2:2:num_subc_usr)', (2:2:num_subc_usr)'];
else                                            % whole resources for data (real, perfect, dd_singletone)
    idx_ofdmsym_pilot_usr = [];                 % pilot symbol index (should be scalar or vector, all elements should be between 1 and 14)
    idx_subc_pilot_usr = [];
end
num_ofdmsym_pilot_usr = length(idx_ofdmsym_pilot_usr);      % number of ofdm symbols with pilots
num_subc_pilot_usr = size(idx_subc_pilot_usr, 1);           % number of pilot subcarriers in 1 ofdm symbol with pilots

% set data index parameter
idx_ofdmsym_data_usr = 1:num_ofdmsym_usr;
idx_ofdmsym_data_usr(idx_ofdmsym_pilot_usr) = [];
idx_subc_data_usr = zeros(num_subc_usr-num_subc_pilot_usr, num_ofdmsym_pilot_usr);
for i = 1:num_ofdmsym_pilot_usr
    tmp_idx_subc_data_usr = (1:num_subc_usr)';
    tmp_idx_subc_data_usr(idx_subc_pilot_usr(:, i)) = [];
    idx_subc_data_usr(:, i) = tmp_idx_subc_data_usr;
end
num_ofdmsym_data_usr = length(idx_ofdmsym_data_usr);        % number of ofdm symbols without pilots (all data subcarrier)
num_subc_data_usr = num_subc_usr-num_subc_pilot_usr;        % number of data subcarriers in 1 ofdm symbol with pilots

% set number of data symbols per user
if strcmp(chest_option, 'dd_singletone')     % tf-domain pilots (some symbols for pilots)
    % otfs pilot setup
    %   - num. otfs pilot(+guard) qam symbols = num. ofdm pilot qam symbols
    %   ex) 600*2 and 86*14
    num_doppler_pilot = 14;             % length of pilot rb to otfs symbols (set even number)
    num_delay_pilot = 43;               % length of pilot rb to subcarriers (set even number)
    num_delay_guard = num_delay_pilot;  % length of guard rb to subcarriers (set even number)
    num_qamsym_usr = (num_subc_usr*num_ofdmsym_usr)-(num_doppler_pilot*(num_delay_pilot+num_delay_pilot));    % number of data symbols per user
else
    num_doppler_pilot = [];
    num_delay_pilot = [];
    num_delay_guard = [];
    num_qamsym_usr = (num_subc_usr*num_ofdmsym_usr)-(num_subc_pilot_usr*num_ofdmsym_pilot_usr);   % number of data symbols per user
end

% set timing parameters
%   - len_frame = len_rb_sym_user * (nfft + num_cp);
%   - tti = (length of 1 symbol + CP) * number of symbol * length of symbol
tti = (nfft + num_cp) * num_ofdmsym_subfrm / sample_rate;               % transmission time interval (sec)
t_subfrm = tti;                                                         % subframe length (sec)

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

nw_num.nfft = nfft;
nw_num.sample_rate = sample_rate;
nw_num.num_cp = num_cp;
nw_num.num_subc_bw = num_subc_bw;                   % updated
nw_num.num_ofdmsym_subfrm = num_ofdmsym_subfrm;
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
nw_num.num_doppler_pilot = num_doppler_pilot;
nw_num.num_delay_pilot = num_delay_pilot;
nw_num.num_delay_guard = num_delay_guard;
nw_num.num_qamsym_usr = num_qamsym_usr;
nw_num.t_subfrm = t_subfrm;

end
