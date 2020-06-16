function rx_sym_dd_ndft_eq = otfs_ch_eq_dd(rx_sym_dd_ndft, ch_est_dd_ndft, num)

% generate block circular channel matrix
len_ch = [num.ndft num.num_ofdmsym_per_subframe];
len_ch_pad_head = len_ch/2;
len_ch_pad_tail = len_ch/2;
len_ch_pad_zero = len_ch;
ch_est_dd_ndft_shift = fftshift(fftshift(ch_est_dd_ndft, 1), 2);
base_ch = repmat(ch_est_dd_ndft_shift(end:-1:1, end:-1:1), 2, 2);
sub_ch = zeros(len_ch(1)+len_ch_pad_zero(1), len_ch(1)+len_ch_pad_zero(1), len_ch(2)+len_ch_pad_zero(2));

for j = 1 : len_ch(2)+len_ch_pad_zero(2)
    for k = 1 : len_ch(1)+len_ch_pad_zero(1)
        sub_ch(k, :, j) = circshift(base_ch(:, j), k-len_ch(1));
    end
end
sub_ch_crop = sub_ch(len_ch_pad_head(1)+1 : end-len_ch_pad_tail(1), 1:len_ch(1), :);
sub_ch_reshape = reshape(sub_ch_crop, len_ch(1), len_ch(1)*(len_ch(2)+len_ch_pad_zero(2)));
new_ch_full = zeros(len_ch(1)*(len_ch(2)+len_ch_pad_zero(2)), len_ch(1)*(len_ch(2)+len_ch_pad_zero(2)));
for k = 1 : len_ch(2)+len_ch_pad_zero(2)
    new_ch_full(len_ch(1)*(k-1)+1 : len_ch(1)*k, :) = circshift(sub_ch_reshape, (k-len_ch(2))*len_ch(1), 2);
end
new_ch = new_ch_full(len_ch_pad_head(2)*len_ch(1)+1 : end-len_ch_pad_tail(2)*len_ch(1), 1:len_ch(2)*len_ch(1)) / sqrt(num.ndft*num.num_ofdmsym_per_subframe);

% equalize channel in dd domain
rx_sym_dd_ndft_eq_vec = new_ch \ rx_sym_dd_ndft(:);
rx_sym_dd_ndft_eq = reshape(rx_sym_dd_ndft_eq_vec, num.ndft, num.num_ofdmsym_per_subframe);

% % dump variables
% assignin('base', 'rx_sym_dd_ndft', rx_sym_dd_ndft);
% assignin('base', 'ch_est_dd_ndft', ch_est_dd_ndft);
% assignin('base', 'base_ch', base_ch);
% assignin('base', 'new_ch_full', new_ch_full);
% assignin('base', 'new_ch', new_ch);
% assignin('base', 'rx_sym_dd_ndft_eq', rx_sym_dd_ndft_eq);

end
