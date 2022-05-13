# gA channel simulation


mpirun -np 32 ./PNP $c_bulk $c_bulk \
 -dof_type P1 -options_file ./solver_Axb/advanced1.options \
 -solute_region 1 \
 -solvent_region 2 \
 -fn_mesh ./mesh/1MAG_2.mesh \
 -fn_pqr ./mesh/1MAG.pqr \
 -fn_bcmap ./mesh/channel_old.bc_dirichlet \
 -fn_ions ./ions_infor/2ions_HR_SRG.ions \
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

 
 #dof_current_gA_PNP-EAFE_KCl_con01_minus-015V.out



# �?code 暂时只能计算pnp的电流，smpnp还需要加上尺寸项产生的电流，待补�?#-Method 2 就是 标准FEM 方法 而已�?define Primitive   2  
# INT Method = Primitive;




# 2ions_HR.ions ---两种离子信息，通过扩散系数值分�?# 3ions_HR.ions ---三种离子信息，通过扩散系数值分�?
 
 #-Stabilize 3 \  ----EAFE 方法
 # 1-Slotboom         2-Primitive （但 mod=6时，亦是它，即poisson方程线性求解）

#-MOD 6 -Stabilize 3  新书�?EAFE ---> 为了算smpnp，但可退化成 pnp---所以先验证pnp

#-MOD 0 -Stabilize 3 原code EAFE (pnp) ---作为结果对照
# a0 2.5  溶剂水尺�?
#-MOD 1  --> smpnp primitive 标准有限元法
#-V_D 1  --> 定常扩散系数D�?#-V_D 2  --> 变扩散系�?D(X)--(V_D_type = 0, 默认，算gA 通道)