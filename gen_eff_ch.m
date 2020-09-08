% gen_eff_ch generates block circulant matrix, or effective
% doubly-convolutional matrix.
%     eff_ch = gen_eff_ch(ch)
%     ch: M by N channel matrix
%     eff_ch: MN by MN channel matrix
%     tip) reshape(permute(reshape(reshape(permute(a(:, :, [1 3 2 4]), [1 3 2]), 8, 3), 4, 2, 3), [1 3 2]), 4, 6)

function eff_ch = gen_eff_ch(ch)

num_el = numel(ch);     % size(ch,1)*size(ch,2)
ch0 = zeros(size(ch,1),size(ch,1),size(ch,2));
for i = 1:size(ch,2)
    ch0(:,:,i) = toeplitz(ch(:,i),circshift(flipud(ch(:,i)),1));
end
idx = reshape(toeplitz(1:size(ch,2),circshift(fliplr(1:size(ch,2)),1)),[],1);
ch1 = reshape(permute(ch0(:,:,idx),[1 3 2]),num_el*size(ch,2),size(ch,1));
ch2 = reshape(permute(reshape(ch1,numel(ch),size(ch,2),size(ch,1)),[1 3 2]),num_el,num_el);
eff_ch = ch2/sqrt(num_el);

end
