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
num_subcarrier = num_rb * 12;
sample_rate = 15e3 * nfft;          % sampling rate (hz)
num_ofdmsym_per_subframe = 14;      % number of ofdm symbols per subframe (pilot + data)

% set data and pilot symbol indices
if strcmp(chest_option, 'tf_pilot') % tf-domain pilots (some symbols for pilots)
    idx_data_ofdmsym = 1 : num_ofdmsym_per_subframe;
    idx_pilot_ofdmsym = [4 11];     % pilot symbol index (should be scalar or vector, all elements should be between 1 and 14)
    idx_data_ofdmsym(idx_pilot_ofdmsym) = [];
else                                % whole resources for data
    idx_data_ofdmsym = 1 : num_ofdmsym_per_subframe;
    idx_pilot_ofdmsym = [];         % pilot symbol index (should be scalar or vector, all elements should be between 1 and 14)
end

% set numerical parameters
num_data_ofdmsym_per_subframe = length(idx_data_ofdmsym);   % number of ofdm symbols per subframe (data)
num_pilot_ofdmsym_per_subframe = length(idx_pilot_ofdmsym); % number of ofdm symbols per subframe (data)
num_ofdmsym_per_subframe = 14;      % number of ofdm symbols per subframe (pilot + data)
num_cp = nfft * 144 / 2048;         % number of samples in cyclic prefix

% set user parameters
num_prb_user = num_rb;          % number of physical resource blocks per user (origianlly 3 rbs)
ndft = num_prb_user * 12;       % number of dft

% set number of data symbols per user
if strcmp(chest_option, 'dd_pilot')     % tf-domain pilots (some symbols for pilots)
    % otfs pilot setup
    % num. otfs pilot(+guard) qam symbols = num. ofdm pilot qam symbols
    % ex) 600*2 and 86*14
    num_pilot_doppler = 14;     % length of pilot rb to otfs symbols (set even number)
    num_pilot_delay = 43;       % length of pilot rb to subcarriers (set even number)
    num_guard_delay = num_pilot_delay;  % length of guard rb to subcarriers (set even number)
    len_rb_sym_user = (ndft*num_ofdmsym_per_subframe)-(num_pilot_doppler*(num_pilot_delay+num_pilot_delay));    % number of data symbols per user
else
    num_pilot_doppler = [];
    num_pilot_delay = [];
    num_guard_delay = [];
    len_rb_sym_user = ndft*num_data_ofdmsym_per_subframe;   % number of data symbols per user
end

% set timing parameters
% len_frame = len_rb_sym_user * (nfft + num_cp);          % frame length
% tti = (length of 1 symbol + CP) * number of symbol * length of symbol
tti = (nfft + num_cp) * 14 / sample_rate;               % transmission time interval (sec)
t_subframe = tti;                                       % subframe length (sec)

% preamble setup (Zadoff-Chu sequence)
% num_zc = 31;
% u = 17;
% m = 0 : num_zc - 1;
% if (mod(num_zc, 2) == 0)
%     z = exp(-1i * pi * u * m ^ 2 / num_zc);
% else
%     z = exp(-1i * pi * u * m(m + 1) / num_zc);
% end
% midamble = z';
% midamble_fseq = [0 ones(1,300) zeros(1,423) ones(1,300)]';
% midamble_fseq = [0 (randi(2,1,300)-1.5).*2 zeros(1,423) (randi(2,1,300)-1.5).*2]';

% % preamble setup (random sequence)
% preamble_seq1 = randi(2, nfft, 1) * 2 - 3;
% preamble_seq2 = randi(2, nfft, 1) * 2 - 3;

% % delay-doppler domain preamble test
% pilot_seq = randi(2, ndft, 1) * 2 - 3;

% % turbo code channel input for awgn channel
% ch_impulse = zeros(nfft, num_sym_per_subframe + 1);
% ch_impulse(725, 1) = sqrt(num_subcarrier * (num_sym_per_subframe + 1));
% 
% mu = 0.01;
% L1 = 20;
% L2 = 3;
% eq_idx1 = 0 : L1 - 1;
% eq_idx2 = 0 : L2 - 1;

% output
nw_num.nfft = nfft;
nw_num.ndft = ndft;
nw_num.len_rb_sym_user = len_rb_sym_user;
nw_num.num_subcarrier = num_subcarrier;
nw_num.sample_rate = sample_rate;
nw_num.num_ofdmsym_per_subframe = num_ofdmsym_per_subframe;
nw_num.num_data_ofdmsym_per_subframe = num_data_ofdmsym_per_subframe;
nw_num.num_pilot_ofdmsym_per_subframe = num_pilot_ofdmsym_per_subframe;
nw_num.idx_data_ofdmsym = idx_data_ofdmsym;
nw_num.idx_pilot_ofdmsym = idx_pilot_ofdmsym;
nw_num.num_cp = num_cp;
nw_num.t_subframe = t_subframe;
% nw_num.len_pilot_otfssym = num_pilot_doppler;
% nw_num.len_pilot_subc = num_pilot_delay;
nw_num.num_pilot_doppler = num_pilot_doppler;
nw_num.num_pilot_delay = num_pilot_delay;
nw_num.num_guard_delay = num_guard_delay;

end
