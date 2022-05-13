# gA channel simulation


mpirun -np 32 ./PNP $c_bulk $c_bulk \
 -dof_type P1 -options_file advanced1.options \
 -solute_region 1 \
 -solvent_region 2 \
 -fn_mesh 1MAG_2.mesh \
 -fn_pqr 1MAG.pqr \
 -fn_bcmap channel_old.bc_dirichlet \
 -fn_ions 2ions_HR_SRG.ions \
 -channel_top 7.0 \
 -channel_bottom 14.0 \
 -V_D 2 \
 -Method 2 \
 -MOD 0 \
 -Stabilize 3 \
 -tol 1.0e-6 \
 -u_bulk -1.50 \
 -Alpha 0.22 \
 -current_calculation 2 \
 -section_value -3.5 \
 -dof_section_value_calculation 1 \
 -out_dof_section_file dof_section_gA_PNP-EAFE_550Na650K_537Cl_3.1a0_con0.1M_minus-0.15V.out \
 -out_current_file dof_current_gA_PNP-EAFE_550Na650K_537Cl_3.1a0_con0.1M_minus-0.15V.out
 
