% gen_conv_ch generates doubly-convolutional matrix from block circulant matrix
%     ch = gen_conv_ch(ch_eff)
%     ch: (M x N) channel matrix
%     ch_eff: (MN x MN) channel matrix
%     tip) reshape(permute(reshape(reshape(permute(a(:, :, [1 3 2 4]), [1 3 2]), 8, 3), 4, 2, 3), [1 3 2]), 4, 6)

function ch = gen_conv_ch(ch_eff, n)

if (size(ch_eff, 1) ~= size(ch_eff, 2)) || (mod(size(ch_eff, 1), n) ~= 0)
    error('Check ''ch_eff'' matrix size.')
end

m = size(ch_eff, 1)/n;

ch = zeros(n, m, size(ch_eff, 1));
% for i = 1:m
%     for j = 1:n
%         k = (i-1)*n+j;
%         ch(:, :, k) = circshift(reshape(ch_eff(:, k), n, m), [-(j-1) -(i-1)])*sqrt(size(ch_eff, 1));
%     end
% end

% for i = 1:m
%     for j = 1:n
%         k = (i-1)*n+j;
%         ch(:, :, k) = circshift(reshape(fliplr(ch_eff(k, :)), n, m), [j i])*sqrt(size(ch_eff, 1));
%     end
% end

for i = 1:m
    for j = 1:n
        k = (i-1)*n+j;
        l = (j-1)*m+i;
        ch(:, :, l) = circshift(reshape(ch_eff(:, k), n, m), [-(j-1) -(i-1)])*sqrt(size(ch_eff, 1));
    end
end

end
