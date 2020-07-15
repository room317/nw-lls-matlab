% IEEE 802.11ad single carrier linear feedback shift register for soft input

function lfsr_out = lfsr_soft(lfsr_in, lfsr_init)

lfsr_shift_reg = lfsr_init;
lfsr_out = zeros(size(lfsr_in));
for lfsr_idx = 1 : length(lfsr_in)
    lfsr_shift_reg = [xor(lfsr_shift_reg(4), lfsr_shift_reg(7)) lfsr_shift_reg(1 : 6)];
    if lfsr_shift_reg(1) == 1
        lfsr_out(lfsr_idx) = (-1) * lfsr_in(lfsr_idx);
    else
        lfsr_out(lfsr_idx) = lfsr_in(lfsr_idx);
    end
%     figure, plot(real(lfsr_in), '-b.'), hold on, plot(real(lfsr_out), ':r.'), hold off, pause
end
