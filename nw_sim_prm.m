function nw_sim = nw_sim_prm(len_tb_bit, num_sim, num, cc, rm)

nw_sim.len_tb_bit = len_tb_bit;
nw_sim.t_sim = rm.t_tb * num_sim;
nw_sim.block_delay = double(cc.C > 1);
nw_sim.cw_delay = double(rm.num_data_sym_per_cw > num.len_rb_sym_user);

end