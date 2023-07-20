#!/bin/sh 
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
# ulimit -s unlimited

# ---------- Ne22_usda --------------
echo "start running log_Ne22_usda_m0p.txt ..."
cat > Ne22_usda_0p.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "Ne22_usda_p.ptn"
  fn_save_wave = "Ne22_usda_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 5.585, -3.826, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 10
  n_restart_vec = 15
&end
EOF
nice ./kshell.exe Ne22_usda_0p.input > log_Ne22_usda_m0p.txt 2>&1 

rm -f tmp_snapshot_Ne22_usda_p.ptn_0_* tmp_lv_Ne22_usda_p.ptn_0_* Ne22_usda_0p.input 


# --------------- transition probabilities --------------

echo "start running log_Ne22_usda_tr_m0p_m0p.txt ..."
cat > Ne22_usda_0p_0p.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "Ne22_usda_p.ptn"
  fn_ptn_r = "Ne22_usda_p.ptn"
  fn_load_wave_l = "Ne22_usda_m0p.wav"
  fn_load_wave_r = "Ne22_usda_m0p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 5.585, -3.826
&end
EOF
nice ./transit.exe Ne22_usda_0p_0p.input > log_Ne22_usda_tr_m0p_m0p.txt 2>&1 

rm -f Ne22_usda_0p_0p.input


echo "Finished computing Ne22_usda"
echo
