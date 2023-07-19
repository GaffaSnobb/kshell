#!/bin/sh 
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
# ulimit -s unlimited

# ---------- V50_gxpf1a --------------
echo "start running log_V50_gxpf1a_m0p.txt ..."
cat > V50_gxpf1a_0p.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "V50_gxpf1a_p.ptn"
  fn_save_wave = "V50_gxpf1a_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 5.585, -3.826, 
  hw_type = 1
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 10
  n_restart_vec = 15
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe V50_gxpf1a_0p.input > log_V50_gxpf1a_m0p.txt 2>&1 

rm -f tmp_snapshot_V50_gxpf1a_p.ptn_0_* tmp_lv_V50_gxpf1a_p.ptn_0_* V50_gxpf1a_0p.input 


# --------------- transition probabilities --------------

echo "start running log_V50_gxpf1a_tr_m0p_m0p.txt ..."
cat > V50_gxpf1a_0p_0p.input <<EOF
&input
  fn_int   = "gxpf1a.snt"
  fn_ptn_l = "V50_gxpf1a_p.ptn"
  fn_ptn_r = "V50_gxpf1a_p.ptn"
  fn_load_wave_l = "V50_gxpf1a_m0p.wav"
  fn_load_wave_r = "V50_gxpf1a_m0p.wav"
  hw_type = 1
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 5.585, -3.826
&end
EOF
nice ./transit.exe V50_gxpf1a_0p_0p.input > log_V50_gxpf1a_tr_m0p_m0p.txt 2>&1 

rm -f V50_gxpf1a_0p_0p.input


echo "Finished computing V50_gxpf1a"
echo
