# Molecular sphere simulation

mpirun -np 32 ./PNP $c_bulk $c_bulk \
 -dof_type P1 -options_file ./solver_Axb/advanced1.options \
 -solute_region 1 \
 -solvent_region 2 \
 -fn_mesh ./mesh/qiaoyu_sphere80.mesh \
 -fn_pqr ./mesh/qiaoyu_sphere80.pqr \
 -fn_bcmap ./mesh/srg_sphere_old.bc_dirichlet \
 -fn_ions ./ions_infor/2ions_HR_SRG.ions \
 -MOD 0 \
 -Stabilize 1 \
 -tol 1.0e-6 \
 -u_bulk 0.0 \
 -Alpha 0.1 \
 -dof_circle_value_calculation 1 \
 -out_dof_circle_file dof_sphere_PNP_FEM-SUPG_15q_0.1M.out 
 
 
 #-dof_circle_value_calculation 1 \
 #-out_dof_section_file dof_section_smpnp.out 
 
#-Stabilize 3 ----IAFEM 方法； -Stabilize 1  SUPG


#-MOD 6 -Stabilize 3  IAFEM ---> 为了算smpnp，但可退化成 pnp---所以先验证pnp
#-MOD 7 --Primitive for smpnp by  Qiaoyu's style

#-MOD 0 -Stabilize 3 原pnp code EAFE  ---作为结果对照 Stabilize 0 --Primitive for pnp
#-MOD 1 -Primitive for smpnp

