% cell-based

function nw_sim = nw_sim_prm_r1(len_tb_bit, num_sim, num, cc, rm)

nw_sim.len_tb_bit = len_tb_bit;
nw_sim.t_sim = rm.N_RB*num.t_usrfrm*num_sim;
nw_sim.block_delay = double(cc.C > 1);
nw_sim.cw_delay = double(rm.num_data_sym_cb > num.num_qamsym_usr);

end