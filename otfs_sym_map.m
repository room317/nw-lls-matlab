% map otfs data and pilot qam symbols

function [tx_sym_dd_ndft, idx_pilot_sym] = otfs_sym_map(data_sym, num, test_dd_pilot)

% generate pilot symbols
pilot_sym = zeros(num.len_pilot_subc, num.len_pilot_otfssym);
pilot_sym((num.len_pilot_subc/2)+1, (num.len_pilot_otfssym/2)+1) = sqrt(num.len_pilot_subc*num.len_pilot_otfssym);
% pilot_sym((num.len_pilot_subc/2)-4, (num.len_pilot_otfssym/2)+1) = sqrt(num.len_pilot_subc*num.len_pilot_otfssym);

% calculate pilot symbol index
% idx_pilot_sym = [ceil((num.ndft-num.len_pilot_subc)/2) (num.num_ofdmsym_per_subframe-num.len_pilot_otfssym)/2];
idx_pilot_sym = [floor((num.ndft-num.len_pilot_subc)/2) (num.num_ofdmsym_per_subframe-num.len_pilot_otfssym)/2];

% test mode (allocate whole resource block to pilot)
if test_dd_pilot
    pilot_sym((num.len_pilot_subc/2)+1, (num.len_pilot_otfssym/2)+1) = sqrt(num.ndft*num.num_ofdmsym_per_subframe);    % remove this!!! for test only!!!
    data_sym = zeros(size(data_sym));
end

% reshape data symbols
data_sym_left = reshape(data_sym(1:num.ndft*idx_pilot_sym(2)), num.ndft, idx_pilot_sym(2));
data_sym(1:num.ndft*idx_pilot_sym(2)) = [];
data_sym_mid = reshape(data_sym(1:(num.ndft-num.len_pilot_subc)*num.len_pilot_otfssym), (num.ndft-num.len_pilot_subc), num.len_pilot_otfssym);
data_sym(1:(num.ndft-num.len_pilot_subc)*num.len_pilot_otfssym) = [];
data_sym_right = reshape(data_sym, num.ndft, idx_pilot_sym(2));

% map resource block
tx_sym_dd_ndft_left = data_sym_left;
tx_sym_dd_ndft_mid = [data_sym_mid(1:idx_pilot_sym(1),:); pilot_sym; data_sym_mid(idx_pilot_sym(1)+1:end,:)];
tx_sym_dd_ndft_right = data_sym_right;
tx_sym_dd_ndft = [tx_sym_dd_ndft_left, tx_sym_dd_ndft_mid, tx_sym_dd_ndft_right];

% % dump variables
% assignin('base', 'pilot_sym', pilot_sym);
% assignin('base', 'data_sym_left', data_sym_left);
% assignin('base', 'data_sym_mid', data_sym_mid);
% assignin('base', 'data_sym_right', data_sym_right);
% assignin('base', 'tx_sym_dd_ndft_left', tx_sym_dd_ndft_left);
% assignin('base', 'tx_sym_dd_ndft_mid', tx_sym_dd_ndft_mid);
% assignin('base', 'tx_sym_dd_ndft_right', tx_sym_dd_ndft_right);
% assignin('base', 'tx_sym_dd_ndft', tx_sym_dd_ndft);

% % plot results
% figure, mesh(1:14, 1:600, abs(fftshift(fftshift(tx_sym_dd_ndft, 1), 2)))   % dd real channel
% xlabel('Doppler'), ylabel('Delay'), zlabel('Amplitude'), title('Pilot Plan')
% pause

end
