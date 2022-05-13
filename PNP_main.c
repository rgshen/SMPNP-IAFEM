/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This program is based on programs by PhD. xieyan with PHG, which	 */
/* solves the PNP equation consist of two equations written below.	 */
/* I'm trying to make this program more legible and scalable.		 */
/* ------------------------------------2015.03.10-Beijing-Scarlet-Maple	 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Nernst-Planck equation:(Size-Modified in last term)			 */
/* in Omega_s(solvent_region)						 */
/* dpi/dt = div_Ji	(Ji is flux defined below, which is current)	 */
/*     Ji = Di * (grad_pi + beta * qi * pi * grad_ph			 */
/*	  + K_i * pi * Sigma(al^3 * grad_pl) / (1 - Sigma(al^3 * pl)))	 */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Poisson equation:							 */
/* if using CGSE							 */
/* in Omega(solute_region & solvent_region)				 */
/* -div_Epsilon( grad_ph ) = 4 * M_PI * Lambda * Sigma( qi * pi )	 */
/*			   + 4 * M_PI * Sigma(qfi * Delta(x - xi))	 */
/*   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   */
/* if using SI								 */
/* in Omega(solute_region & solvent_region)                              */
/* -div_Epsilon( grad_ph ) * Epsilon_vaccum = Lambda * Sigma( qi * pi )	 */
/*                         		    + Sigma(qfi * Delta(x - xi)) */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Boundary condition:	(there would be more bdrys in special case)	 */
/*	     ph = 0			on patial_Omega			 */
/*	     pi = pi_bulk		on patial_Omega			 */
/* (nonreactive)	 Neumann					 */
/*	(n, Ji) = 0			on bmlc_bdry			 */
/* (complete reactive)	 Neumann & Dirichlet				 */
/*	(n, Ji) = 0			on bmlc_bdry \ reac_bdry	 */
/*	     pi = 0			on reac_bdry			 */
/* (incomplete reactive) Neumann & Robin				 */
/*	(n, Ji) = 0			on bmlc_bdry \ reac_bdry	 */
/*      (n, Ji) = rr_Alpha * pi		on reac_bdry (rr_Alpha = 1.0e8)	 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* We can choose different kind of ways to solve PNP in this program	 */
/* system of units: 0-CGSE, 1-SI (default -SOU 0)			 */
/* method:	0-for test		(default -Method 1)		 */
/* 	 	1-Slotboom 		    2-Primitive 			 */
/*		3-Time_Slotboom	 	4-Time_Primitive		 */
/*		5-Mono_Slotboom		6-Mono_Primitive		 */
/* stabilization scheme: 0-NONE, 1-SUPG, 2-PRFB (default -Stabilize 0)	 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* PNP_func.h 		declares all the functions in program		 */
/* PNP_coefficient.c 	uses SI and CGSE to recalculate coefficients	 */
/* PNP_build_solver.c 	builds differents solvers for u and H, P, p	 */
/* PNP_more.c 		gives basic functions, can be defined by users,  */
/* 			to simplify integrations in building solvers	 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Attention!								 */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* All the coefficients are based on the two equation before.		 */
/* Using different systems of units, the values of functions without	 */
/* derivatives are the same as in the functions in equation, but the 	 */
/* system of units changes from cm(or m) to A, so every derivative of A  */
/* generates a cm2A = 1.0e8 (or m2A = 1.0e10), for example in CGSE	 */
/* df/dr = dx/dr * df/dx = cm2A * df/dx (x in A, r in cm), same in SI.	 */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Mainly, all {u, p, P} are used for new calculation step, while all	 */
/* {TMP_u, TMP_p, TMP_P} are used for last step both steady state and	 */
/* time-dependent situation.						 */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* In order to improve the scalability, all solver should prepare for 	 */
/* calculating any number of ions with any kind of mesh, so information	 */
/* about the ions and mesh should input before program in *.bc and *.ion */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Many command line options are offered in README or PNP_coefficient.c	 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*-------------  rgshen 全部改造于 2021-6-13 --为了计算 SMPNP------------ */

#include"phg.h"
#include"PNP_func.h"
#include"PNP_coefficient.c"
#include"PNP_build_solver.c"
#include"PNP_more.c"

#include"PNP.h" // rgshen



int main(int argc, char **argv){
	INT mem_max = 20000;
	size_t mem_peak;
	//size_t mem, mem_peak; //rgshen
	SOLVER *solver;
	GRID *g;
	DOF *size_p; //size_p is always < 1 - Sigma(a[n]^3, p[n]) >
	DOF *ph, *ph_u_D, *u, *u_D, *tu, *tmp_u, *TMP_u, *tmp_u_D, *TMP_u_D, *G, *H, *grad_G, *grad_H, *indicator, *error, *chp0, *chp1;
	DOF *p[MNION], *P[MNION], *TMP_p[MNION], *TMP_P[MNION];
	DOF *D[MNION];
    DOF *Epsilon;
    DOF *Psi; //rgshen
	DOF *size_c[MNION]; //rgshen-2021-6-11 为了smpnp primitive 方法而设

	INT i;

	/* command line options preprocess */
	initialize_options(mem_max);

/* Initialize the constants, mesh, fixed charges, concentration */
	phgInit(&argc, &argv);
	phgPrintf("\nPoisson-Nernst-Planck Solver (version 0.2.8)\n\n");
	phgPrintf("--------------------------\n");
	Read_parameter(read_parameter);
	switch(Method) {
		case Primitive:		phgPrintf("     Primitive Method    \n"); break;
		case Slotboom:      phgPrintf("     Slotboom  Method    \n"); break; //rgshne2021-4-25
	}
	phgPrintf("--------------------------\n");


	Init_coefficient();

	/* if using surface charge, *.pqr is not necessary */
	if(surface_charge) {
		phgPrintf("\nUse surface charge\n");
		phgPrintf("density = %.4e C/m^2\n", scd);
	}
	else {
		Read_pqr(fn_pqr);
	}
	Read_bcmap(fn_bcmap); //读取边界标记信息
	Read_ions(fn_ions);  //读取离子信息：化合价、边界体积浓度c_bulk、扩散系数 Di 等


	/* adjust concentration data */

	if(input_data) {
		//problem: wrong number of bulks
		if(argc == NION + 1||1) {
			for(i = 0; i < NION; i++) {
				ion[i].c = atof(argv[i + 1]);
			}
		}
		else
			phgError(1, "\nnumber of input values %d is not equal to number of ions %d", argc - 1, NION);
	}

	Is = 0.0;
	phgPrintf("NION = %d\n", NION);

        if(MOD == 0) {
            phgPrintf("\n MOD = %d :-->Solve the standard PNP equations\n\n", MOD);
        }
        else if(MOD == 1) {
            phgPrintf("\n MOD = %d :-->>Solve SMPNP--Size effect considered\n\n", MOD);
        }
		
		//else if(MOD == 2) {
        else if(MOD == 2 ) { //rgshen VDPNP -eafe 法
                phgPrintf("\n MOD = %d :-->Variable dielectric coefficient considered\n\n", MOD);
        }
		
		else if(MOD == 6) {  //written by rgshen 2021-4-21
                phgPrintf("\n MOD = %d :-->Solve SMPNP_with_Slotboomtransform_averageEAFE\n\n", MOD);
		}
		else if(MOD == 7) {  //written by rgshen 2021-7-1, Qiaoyu's style FEM for smpnp
                phgPrintf("\n MOD = %d :-->Solve SMPNP_with_primitive by Qiaoyu\n\n", MOD);
		}

	for(i = 0; i < NION; i++) {
		Is += 0.5 *  ion[i].c * ion[i].Z * ion[i].Z;
		/* system of units change mol/L ----> N/L ----> N/A^3 or N/nm^3 ----> N/cm^3 or N/m^3 */
		    bulk[i] = ion[i].c * NA / Pow(change_mesh_unit, 3) * Pow(change_unit, 3);   // 这是liuxuejiao的方式--已经验证二者结果一致
			  //\bulk[i] = ion[i].c ;   // 直接 锁定浓度单位为 mol/L，用到边界赋值：func_pb
            if(MOD == 0) {
                     //phgPrintf("ion[%d].c = %.2e mol/L\n", i, ion[i].c);
			    phgPrintf("ion[%d].c = %.2e mol/L     ion[%d].a = %.2e A\n", i, ion[i].c, i, ion[i].a);  //rgshen
             }
             else if(MOD == 1 || MOD == 6 || MOD == 7) {
                    phgPrintf("ion[%d].c = %.2e mol/L     ion[%d].a = %.2e A\n", i, ion[i].c, i, ion[i].a);
             } 
	}

	phgPrintf("a0 = %.2e A \n", a0); // rgshen2021-5-28
	phgPrintf("u_bulk = %.2e V\n", u_bulk);
	phgPrintf("Is = %.2e mol/L\n", Is);
	phgPrintf("chan_ratio = %.4e \n", chan_ratio);
	phgPrintf("Alpha = %.4e \n", Alpha);

	g = phgNewGrid(-1);
	phgImportSetBdryMapFunc(bc_map);
	if (!phgImport(g, fn_mesh, FALSE))
		phgError(1, "can't read file \"%s\".\n", fn_mesh);

	/* check the grid */
	if(face_calculation) {
		INT face_num, tmp_num;
		FLOAT face_area, tmp_area;
		face_measure(face_calculation, g, &tmp_num, &tmp_area);
		MPI_Reduce(&tmp_num, &face_num, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&tmp_area, &face_area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		phgPrintf("\nthis boundary has %d faces and area %lf\n", face_num, face_area);
	}

	if(region_calculation) {
		DOF *test_I = phgDofNew(g, DOF_DEFAULT, 1, "I", DofInterpolation);
		char *vtk_name = "PNP_region_test";
		change_file_path(tmp_path, &vtk_name);
		phgDofSetDataByValue(test_I, 1);
		phgPrintf("\nCreating \"%s\".\n", phgExportVTK(g, vtk_name, test_I, NULL));
		phgDofFree(&test_I);
		INT region_num, tmp_num;
		FLOAT region_volume, tmp_volume;
		region_measure(region_calculation, g, &tmp_num, &tmp_volume);
		MPI_Reduce(&tmp_num, &region_num, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&tmp_volume, &region_volume, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		phgPrintf("\nregion %d has %d elements and volume %lf\n", region_calculation, region_num, region_volume);
	}

	if(region_level >= 0) {
		mark_region(g, refine_region, region_level);
		phgRefineMarkedElements(g);
	}
	phgRefineAllElements(g, all_level);

	/* After every refinement, new grid should be rebalanced */
	if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
		phgPrintf("\nRepartition mesh, load imbalance: %lg\n", (double)g->lif);

	for(which_ion = 0; which_ion < NION; which_ion++) {
		K_i[which_ion] = Pow(ion[which_ion].a / a0, 3);
		phgPrintf("K_i[%d] = %e\n", which_ion, K_i[which_ion]); //rgshen2021-5-28
		phgPrintf("\n--------\n"); //rgshen2021-5-28
	}
	for(which_ion = 0; which_ion < NION; which_ion++) {
		phgPrintf("ion[%d].Z = %e\n", which_ion, ion[which_ion].Z);
	}

	/* initialize dofs */
	u = phgDofNew(g, DOF_DEFAULT, 1, "u", DofInterpolation);
	u_D = phgDofNew(g, DOF_DEFAULT, 1, "u_D", DofInterpolation);
	tu = phgDofNew(g, DOF_DEFAULT, 1, "tu", DofInterpolation);
	tmp_u = phgDofNew(g, DOF_DEFAULT, 1, "tmp_u", DofInterpolation);
	TMP_u = phgDofNew(g, DOF_DEFAULT, 1, "TMP_u", DofInterpolation);
	tmp_u_D = phgDofNew(g, DOF_DEFAULT, 1, "tmp_u_D", DofInterpolation);
	TMP_u_D = phgDofNew(g, DOF_DEFAULT, 1, "TMP_u_D", DofInterpolation);

	/*Initialize G, Grad_G, H*/
	G = phgDofNew(g, DOF_DEFAULT, 1, "G", func_G);
	grad_G = phgDofNew(g, DOF_DEFAULT, 3, "grad_G", func_grad_G);
	H = phgDofNew(g, DOF_DEFAULT, 1, "H", DofInterpolation);
	/**************************************************/

	size_p = phgDofNew(g, DOF_DEFAULT, 1, "size_p", DofInterpolation);
	indicator = phgDofNew(g, DOF_P0, 1, "indicator", DofNoAction);
	error = phgDofNew(g, DOF_P0, 1, "error", DofNoAction);
	chp0 = phgDofNew(g, DOF_DEFAULT, 1, "chp0", DofInterpolation);
	chp1 = phgDofNew(g, DOF_DEFAULT, 1, "chp1", DofInterpolation);

	Psi = phgDofNew(g, DOF_DEFAULT, 1, "Psi", DofInterpolation); //rgshen2021-4-27

	if(NION < MNION) {
		p[NION] = NULL;
		P[NION] = NULL;
		TMP_p[NION] = NULL;
		TMP_P[NION] = NULL;

		size_c[NION] = NULL; //rgshen2021-6-11
	}
	else phgError(1, "\nerror: The number of DOFs is more than MNION, set a larger MNION!");

	/* name the dofs which may output */
	for(i = 0; i < NION; i++) {
		char p_name[5] = "c[0]";    // 小写p ->c[]
		char P_name[5] = "C[0]";    // 大写P ->C[]
		p_name[2] += i;
		P_name[2] += i;
		p[i] = phgDofNew(g, DOF_DEFAULT, 1, p_name, DofInterpolation);
		P[i] = phgDofNew(g, DOF_DEFAULT, 1, P_name, DofInterpolation);
		TMP_p[i] = phgDofNew(g, DOF_DEFAULT, 1, "tmp", DofInterpolation);
		TMP_P[i] = phgDofNew(g, DOF_DEFAULT, 1, "tmp", DofInterpolation);
		size_c[i] = phgDofNew(g, DOF_DEFAULT, 1, "size_c", DofInterpolation); //rgshen2021-6-11
		if(V_D) {
			which_ion = i;
			char D_name[5] = "D[0]";
			D_name[2] += i;
			D[i] = phgDofNew(g, DOF_DEFAULT, 1, D_name, func_D);
		}
	}

	/* Initialize dofs as 0 */
	phgDofSetDataByValue(H, 0.0);
	phgDofSetDataByValue(u, 0.0);
	phgDofSetDataByValue(u_D, 0.0);
	phgDofSetDataByValue(tu, 0.0);
	phgDofSetDataByValue(tmp_u, 0.0);
	phgDofSetDataByValue(TMP_u, 0.0);
	phgDofSetDataByValue(tmp_u_D, 0.0);
	phgDofSetDataByValue(TMP_u_D, 0.0);
	phgDofSetDataByValue(chp0, 0.0);
	phgDofSetDataByValue(chp1, 0.0);
	phgDofSetDataByValue(Psi, 0.0); //rgshen2021-4-27

	for(i = 0; i < NION; i++) {
		phgDofSetDataByValue(TMP_p[i], 0.0);
		phgDofSetDataByValue(TMP_P[i], 0.0);
		phgDofSetDataByValue(p[i], 0.0);
		phgDofSetDataByValue(P[i], 0.0);
		phgDofSetDataByValue(size_c[i], 0.0); //rgshen2021-6-11
	}

        Epsilon = phgDofNew(g, DOF_DEFAULT, 1, "Epsilon", DofInterpolation);
        phgDofSetDataByValue(Epsilon, Epsilon_s);


	phgPrintf("\nAll the unit of derivatives and integrates below is in %f Angstrom.\n", mesh_unit);

	//get the solution of H (Assert no need to change solver, so fixed in PCG)
	//surface charge can replace G and H
	if(!surface_charge) {
		phgPrintf("||G||_2 = %.4e, ||grad_G||_2 = %.4e\n", phgDofNormL2(G), phgDofNormL2(grad_G));
		solver = phgSolverCreate(SOLVER_DEFAULT, H, NULL);
		//solver = phgSolverCreate(SOLVER_PCG, H, NULL);
		Build_H_Solver(solver, H);   //Solving Harmonic equantion---rgshen
		phgSolverSolve(solver, TRUE, H, NULL);
		phgPrintf("||H||_2 = %.4le, ", phgDofNormL2(H));
		grad_H = phgDofGradient(H, NULL, NULL, NULL);
		phgPrintf("||grad_H||_2 = %.4e, ", phgDofNormL2(grad_H));
		phgPrintf("nits = %d, residual = %.4le\n", solver->nits, solver->residual);
		phgSolverDestroy(&solver);
	}
	else grad_H = phgDofGradient(H, NULL, NULL, NULL);
	phgPrintf("\nAll coefficients initialization completed\n");

//get_size_p(&size_p, p, func_p_k);  //

	INT N = DofGetDataCount(u);
	FLOAT t0 = phgGetTime(NULL);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* before every algorithm, there is a initialization prepared for special situation */


/* ---------------开始非线性迭代循环-------------------------rgshen */
	while( count < 10 * loop_maxit) {
		int reflag = 0;//for refinement use
		count ++;

		phgPrintf("STEP = %d\n",count); ///rgshen 6-13
      //-------------------------------------------------------
		if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
			phgPrintf("Repartition mesh, load imbalance: %lg\n",(double)g->lif);
		    phgPrintf("%d DOF, %d elements, %d submeshes, load imbalance: %lg\n",
			DofGetDataCountGlobal(u), g->nleaf_global, g->nprocs, (double)g->lif);

		N = DofGetDataCount(u);
	  //------------------------------------------------------

        //Solving Poisson equation --->eletrostatic potential u
     /*-------------------------------------------------------------------*/
		//p[i]->u
		phgDofAXPBY(1.0, u, 0.0, &TMP_u); //TMP_u = u_old
		solver = phgSolverCreate(SOLVER_DEFAULT, u, NULL);
		if(MOD == 0 || MOD == 1 || MOD == 6 || MOD == 7) {  // 不管什么方程、方法poisson方程都线性求解6-13

			Build_Primitive_Solver(solver, u, size_p, p, grad_G, grad_H);
			    //phgPrintf(" MOD  = %d, I am the standard PNP--(smpnp-EAFE) solver for poisson equation", MOD);// rgshen
		}		
		phgSolverSolve(solver, 1, u, NULL);

		phgDofAXPBY(1.0 - Alpha, TMP_u, Alpha, &u); 
		phgDofAXPY(-1.0, u, &TMP_u);
		    phgPrintf("\n||u||_2 = %.4e, ||delta_u||_2 = %.4e, nits = %d, residual = %.4e\n",
			            phgDofNormL2(u), phgDofNormL2(TMP_u), solver->nits, solver->residual);
		phgSolverDestroy(&solver);

     /*-------------------------------------------------------------------*/
        // Solving Nernst-Planck equations --->concentrations p1, p2, ..., pk

		//u -> p[i]
		for(i = 0; i < NION; i++) {
			which_ion = i;
			if(bulk[i] != 0) {
				phgDofAXPBY(1.0, p[i], 0.0, &TMP_p[i]); //TMP_p[i] is the old p[i]

				solver = phgSolverCreate(SOLVER_DEFAULT, p[i], NULL);

				if(MOD == 0) { // pnp 方程标准有限元法 FEM 和 IAFEM(Stabilize == 3时)
					Build_p_Solver(solver, u, size_p, TMP_p, D[i]);  //可调IAFEM方法
					    phgPrintf("I AM FOR PNP--2021-6-11");
				}
                else if(MOD == 1) {//标准 FEM(primitive) for smpnp -通过2021-6-11
			        get_size_c(&size_c[i], p);    /*-->k_i / (1 - \sum a_l^3 p_l) */	  
                    SMPNP_NP_Slover_Primitive(solver, u, size_c, TMP_p, D[i], size_p);   //标准 FEM for smpnp -通过2021-6-11
					    phgPrintf("I AM Primitive FOR SMPNP--2021-6-11");
				}
				else if(MOD == 7) {	//标准 FEM(primitive) for smpnp, Qiaoyu's style
					SMPNP_PRIMITIVE_NP(solver, u, p[i], TMP_p, D[i]);  //这是qiaoyu标准FEM 求解smpnp--> D[i]  只有求 通道问题时才有用，并修改内部内容
                                        phgPrintf("I AM Primitive by Qiaoyu FOR SMPNP--2021-7-1");
			        }
				 //---MOD = 6, SMPNP-IAFEM-----------------------------------
				else if(MOD == 6) { /*SMPNP_IAFE-2021-4-21 */
					get_size_u_p(&Psi, u, p);  // 对于Psi = u （即标准pnp, 正确）；psi= zi*u - ki* ln((1 - Sigma(al^3 * pl)))
					SMPNP_NP_IAFEM(solver, Psi, TMP_p, D[i]);
						phgPrintf("i am  here for smpnp-IAFE 5-28---333\n");
				}
				
                //-----------------------------------------------
				//p[i] is new
				phgSolverSolve(solver, 1, p[i], NULL);

				dof_positive(p[i]); // if p[i]<0,then let p->data[i] = 1e-8; 

				phgDofAXPBY(1.0 - Alpha, TMP_p[i], Alpha, &p[i]);
                phgDofAXPY(-1.0, p[i], &TMP_p[i]);
                    phgPrintf("||c[%d]||_2 = %.4e, ||delta_c[%d]|| = %.4e, nits = %d, residual = %.4e\n",
                                i, phgDofNormL2(p[i]) / NA, i, phgDofNormL2(TMP_p[i]), solver->nits, solver->residual / NA);

					/* phgPrintf("||c[%d]||_2 = %.4e, ||delta_c[%d]|| = %.4e, nits = %d, residual = %.4e\n",
                                   i, phgDofNormL2(p[i]), i, phgDofNormL2(TMP_p[i]), solver->nits, solver->residual); */	//rgshen

				phgSolverDestroy(&solver);

				//get_size_p(&size_p, p, func_p_k); // 如果用的话 需要更新它(经验证，跟我用size_c 差不多)

			}
			else {
				phgPrintf("||c[%d]||_2 = %.4e\n", i, phgDofNormL2(p[i]) / NA);
				//\phgPrintf("||c[%d]||_2 = %.4e\n", i, phgDofNormL2(p[i])); //rgshen
			}
		} // 结束离子循环
            /**************************************************/
		    //----跳出判断
		if(phgDofNormL2(TMP_u) < phgDofNormL2(u) * tol) {
			phgPrintf("\n ------------------------------------\n");
            phgPrintf("\n||u||_2 = %.4e, ||delta_u||_2 = %.4e, STEP = %d \n",
			            phgDofNormL2(u), phgDofNormL2(TMP_u),  count);
			break;
		}
	} //结束非线性迭代



/* informatin printed for final result */

/* calculate energy */
	if(energy_calculation)
		Calculate_energy(u, H);

/* calculate current */
	phgDofNormL2(u);
        get_size_p(&size_p, p, func_p_k);
	if(MOD != 4) {
		get_Epsilon(&Epsilon, p);  //.....
	}
	if(current_calculation == 1) { //分别计算 x，y，z三个方向的电流
		if(Method == Primitive) {
            //for(i = 0; i < 3; i++) {
            for(i = 0; i < 99; i++) {  //for cylinder_nanopore
				//section_value = -25 + (i+1) * 0.5;
			    //section_value = -5 + i * 5;   //for nanopore    //纳米孔
				//section_value = 0.0 + i * 5;   //for 2JK4
			  section_value = -3.5 + i * 3.5; //for gA    //gA 通道 rgshen
				//section_value = 30 + i * 5; //for connexin
        	    VD_Calculate_current_u_D_born(solvent_region, current_direction, u, u_D, p, section_value, size_p, Epsilon);
			}
		}
	}
	else if(current_calculation == 2) { //只计算 z 方向(z_axis)的电流：section_value=0.0 由开关传入
		VD_Calculate_current_u_D_born(solvent_region, current_direction, u, u_D, p, section_value, size_p, Epsilon);   //暂时只能计算pnp方程电流哈
	}

        for(i = 0; i < N; i++) {
            //change N/cm^3 or N/m^3 -> N/A^3 -> mol/L
            for(which_ion = 0; which_ion < NION; which_ion++) {
                p[which_ion]->data[i] *= Pow(change_mesh_unit, 3) / Pow(change_unit, 3) / NA;
            }
        }     

	if(chemical_potential_calculation) {
		Calculate_chemical_potential(&chp0, u, p, Epsilon, 0);
		Calculate_chemical_potential(&chp1, u, p, Epsilon, 1);
	}
	if(potential_energy_calculation) {
		Calculate_potential_energy(&chp0, u, Epsilon, 0);
		Calculate_potential_energy(&chp1, u, Epsilon, 1);
	}

/* calculate dof surface value*/
	if (dof_circle_value_calculation) {
		/*FLOAT ra =0.0;
		FLOAT rb =2.0; //for the cylinder's radius = 2 A
		FLOAT tr =0.1;
		FLOAT radius, dof_value_u, dof_value_u_D, dof_value_p1, dof_value_p2, dof_value_p3, dof_value_eps, value_result_phi;
        	INT J = (INT)((rb-ra)/tr);
		int j;*/
//---------------沈瑞刚改------------------------
		FLOAT ra = 10; //球心位置，小球半径为 1 A
		FLOAT rb = 40.5; //40.05; //69.05; //球心到盒子外边界的距离
		FLOAT tr = 0.5; //0.05;	//取值步长
		FLOAT radius, dof_value_u, dof_value_u_D, dof_value_p1, dof_value_p2, dof_value_p3, dof_value_eps, value_result_phi;
        	//INT J = (INT)((rb-ra)/tr) + 1; //为了取值到盒子边界
        	INT J = (INT)((rb-ra)/tr); //为了取值到盒子边界
		int j;
		radius = ra;
//--------------------------------------------------------------
//----------沈瑞刚添加--为了存储分子小球数据--------------------
		if (g->rank == 0 && out_dof_circle_file != NULL) {
			FILE *f_in;
            f_in = fopen(out_dof_circle_file, "a");
            if(f_in == NULL) {
                phgError(1, "Can not open output file %s.\n", out_dof_circle_file);
            }
            fprintf(f_in, "z		potential(V)	epsilon		p0(mol/L)	p1(mol/L)\n");
            fclose(f_in);
		}
//---------------------------------------------------------------
		    for(j = 0; j < J; j++) {  //for sphere rgshen-20201118
			//section_value = 0.0;
            //for(i = 0; i < 60; i++) {  //for k_channel
            //for(i = 0; i < 43; i++) {  //for gA
            //for(i = 0; i < 101; i++) {  //for cylinder_nanopore
			   //section_value = 0.0 + i * 0.5; //for the cylinder's length = 20A (0, 20)  //rgshen
			    //section_value = -25 + i * 0.5; //for the cylinder's length = 50A (-25, 25)
			    //section_value = -14 + i * 0.5; //for the gA's length = 21A (-14, 7)
			    //section_value = -15 + i * 0.5; //for the k_channel's length = 30A (-15, 15)
			dof_value_u = Calculate_dof_circle_value(solvent_region, u, radius, section_value);
				//dof_value_u_D = Calculate_dof_circle_value(solvent_region, u_D, radius, section_value);
			value_result_phi = dof_value_u / (Beta * ec);  // 转换为电势输出	rgshen

				dof_value_p1 = Calculate_dof_circle_value(solvent_region, p[0], radius, section_value);
				dof_value_p2 = Calculate_dof_circle_value(solvent_region, p[1], radius, section_value);
				if (NION == 3)
					dof_value_p3 = Calculate_dof_circle_value(solvent_region, p[2], radius, section_value);
				    dof_value_eps = Calculate_dof_circle_value(solvent_region, Epsilon, radius, section_value);
	                if (g->rank == 0 && out_dof_circle_file != NULL) {
        	            FILE *f_in;
                	    f_in = fopen(out_dof_circle_file, "a");
                        if(f_in == NULL) {
							phgError(1, "Can not open output file %s.\n", out_dof_circle_file);
                        }
					    if (NION == 3) {
                            //\fprintf(f_in, "%f  %f  %lf  %lf  %lf  %lf  %lf\n", radius, section_value, dof_value_u, dof_value_eps, dof_value_p1, dof_value_p2, dof_value_p3);
							fprintf(f_in, "%f  %lf    %lf   %lf    %lf    %lf\n", radius, value_result_phi, dof_value_eps, dof_value_p1, dof_value_p2, dof_value_p3);//rgshen--sphere
					    }
					    else {
                            fprintf(f_in, "%f   %lf  %lf  %lf  %lf\n", radius, dof_value_u, dof_value_eps, dof_value_p1, dof_value_p2);
							//fprintf(f_in, "%f	%lf	%lf	%lf	%lf\n", radius, value_result_phi, dof_value_eps, dof_value_p1, dof_value_p2); //rgshen--sphere
					    }
						fclose(f_in);
                    }
			//}
			radius += tr;
		}
	}

/* calculate section value*/
    if(dof_section_value_calculation) {
		FLOAT value_result_u, value_result_phi, value_result_u_D, value_result_p1, value_result_p2, value_result_p3;
		FLOAT value_result_eps, value_chp0, value_chp1;
        if (g->rank == 0 && out_dof_section_file != NULL) {
			FILE *f_in;
            f_in = fopen(out_dof_section_file, "a");
            if(f_in == NULL) {
                phgError(1, "Can not open output file %s.\n", out_dof_section_file);
            }
            fprintf(f_in, "z		potential(V)	epsilon		p0(mol/L)	p1(mol/L)\n");
            fclose(f_in);
		}
        //\for(i = 0; i < 181; i++) {  //for gA (-45,45)
		for(i = 0; i < 81; i++) {  //for gA (-20,20)
        //for(i = 0; i < 201; i++) {  //for cylinder_nanopore
            //for(i=0;i <102; i++){
        //\for(i = 0; i < 199; i++) {  //for 1bl8 channel
        //\for(i = 0; i < 81; i++) {  //for 1bl8 channel
		    //section_value = -20 + i*0.5; //for test
           //section_value = -50 + i*0.5; //for cylinder_nanopore
		   //section_value = 0.0 + i*0.5; //for cylinder_nanopore  //rgshen
            //section_value = -49.5 + i*0.5; //for 1bl8 channel
          //\section_value = -45 + i * 0.5;  //for gA // 因子0.5决定数值间距，也决定取多少个数据，取值：(-45,45)
		  section_value = -20 + i * 0.5;  //for gA // 因子0.5决定数值间距，也决定取多少个数据，取值：(-20,20)
            value_result_u = Calculate_dof_value(solvent_region, current_direction, u, section_value);
			value_result_phi = value_result_u / (Beta * ec);
			//value_result_u_D = Calculate_dof_value(solvent_region, current_direction, u_D, section_value);
            value_result_eps = Calculate_dof_value(solvent_region, current_direction, Epsilon, section_value);
            if(chemical_potential_calculation || potential_energy_calculation) {
                value_chp0 = Calculate_dof_value(solvent_region, current_direction, chp0, section_value);
                value_chp1 = Calculate_dof_value(solvent_region, current_direction, chp1, section_value);
            }
            value_result_p1 = Calculate_dof_positive_value(solvent_region, current_direction, p[0], section_value);
            value_result_p2 = Calculate_dof_positive_value(solvent_region, current_direction, p[1], section_value);
            if (NION == 3)
                value_result_p3 = Calculate_dof_positive_value(solvent_region, current_direction, p[2], section_value);
            if (g->rank == 0 && out_dof_section_file != NULL) {
                FILE *f_in;
                f_in = fopen(out_dof_section_file, "a");
                if(f_in == NULL) {
                    phgError(1, "Can not open output file %s.\n", out_dof_section_file);
                }
                if(NION == 3) {
                    if(chemical_potential_calculation || potential_energy_calculation) {
                        fprintf(f_in, "%f	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", section_value, value_result_u, value_result_eps, value_result_p1, value_result_p2, value_result_p3, value_chp0, value_chp1);
                    }
                    else {
                        //fprintf(f_in, "%f	%lf	%lf	%lf	%lf	%lf\n", section_value, value_result_u, value_result_eps, value_result_p1, value_result_p2, value_result_p3);
						fprintf(f_in, "%f	%lf	%lf	%lf	%lf	%lf\n", section_value, value_result_phi, value_result_eps, value_result_p1, value_result_p2, value_result_p3);
                    }
                }
                else { //通常只计算两个离子情形                             //静电势phi
                    fprintf(f_in, "%f	%lf	%lf	%lf	%lf\n", section_value, value_result_phi, value_result_eps, value_result_p1, value_result_p2);
                }
                fclose(f_in);
            }
        }
    }


/* test calculation */
	if(0)
		for(which_ion = 0; which_ion < NION; which_ion++)
			phgPrintf("!!!!%le", P_B_Dof_Int(solvent_region, u));

/* calculate reactive rates */

/* change system of units for output */
	ph = phgDofNew(g, DOF_DEFAULT, 1, "ph", DofInterpolation);
	ph_u_D = phgDofNew(g, DOF_DEFAULT, 1, "ph_u_D", DofInterpolation);
	for(i = 0; i < N; i++) {
		ph->data[i] = u->data[i] / (Beta * ec); /* kBT ~= 0.6 */
		ph_u_D->data[i] = u_D->data[i] / (Beta * ec); /* kBT ~= 0.6 */
		//change N/cm^3 or N/m^3 -> N/A^3 -> mol/L
		for(which_ion = 0; which_ion < NION; which_ion++) {
			P[which_ion]->data[i] *= Pow(change_mesh_unit, 3) / Pow(change_unit, 3) / NA;
		}   
	}

/* about boundary */
	DOF *epsl, *bdry;
	//epsl = phgDofNew(g, DOF_P0, 1, "epsl", DofNoAction);
	bdry = phgDofNew(g, DOF_DEFAULT, 1, "bdry", DofInterpolation);
	epsl = phgDofNew(g, DOF_DEFAULT, 1, "epsl", DofInterpolation);
        //epsl = phgDofCopy(Epsilon, NULL, NULL, NULL);
	phgDofSetDataByValue(epsl, 0.0);
	phgDofSetDataByValue(bdry, 0.0);
	get_Epsilon_all(&epsl, Epsilon, p);
	See_bdry(&bdry);
/* vtk output */
	if(off_vtk) {
		count = -1;
		dof_save(dof_file, &count, ph, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], NULL);
	}
	else if(Method == Primitive) {
		//char *vtk_name = "PNP_Primitive_result";
		char vtk_name[100];
		//sprintf(vtk_name, "PNP-M%d_nE%d_a456_cc-2-30-%4.2f.vtk", MOD, nonlinear_Epsilon, u_bulk);
		//sprintf(vtk_name, "PNP-M%d-%4.2f.vtk", MOD, u_bulk);

		sprintf(vtk_name, "smpnp-gA-M%d-%4.2f.vtk", MOD, u_bulk); //rgshen

		//sprintf(vtk_name, "VDPNP_type1-sphere-M%d-Vol%4.2f-Con%4.2f.vtk", MOD, u_bulk, Is);  // 固定浓度改变电压，用这个
		//sprintf(vtk_name, "VDPNP-sphere-M%d-Con%4.2f-Vol%4.2f.vtk", MOD, Is, u_bulk); // 固定电压改变浓度，用这个
		change_file_path(tmp_path, &vtk_name);
            if(V_D) {
                phgPrintf("Creating \"%s\".\n\n", phgExportVTK(g, vtk_name, ph, u, D[0], G, H,
                            p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], NULL));
            }
            else {
                phgPrintf("Creating \"%s\".\n\n", phgExportVTK(g, vtk_name, ph, u, G, H,
                            p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], NULL));
            }
                //phgPrintf("Creating \"%s\".\n\n", phgExportVTK(g, vtk_name, ph, u, epsl, bdry, G, H,
                //phgPrintf("Creating \"%s\".\n\n", phgExportVTK(g, vtk_name, ph, u, epsl, bdry, D[0], G, H,

/*
            sprintf(vtk_name, "PNP-1bl8-u-%4.2f.vtk", u_bulk);
            change_file_path(tmp_path, &vtk_name);
            phgPrintf("Creating \"%s\".\n\n", phgExportVTK(g, vtk_name, u, NULL));

            sprintf(vtk_name, "PNP-1bl8-c1-%4.2f.vtk", u_bulk);
            change_file_path(tmp_path, &vtk_name);
            phgPrintf("Creating \"%s\".\n\n", phgExportVTK(g, vtk_name, p[0], NULL));

            sprintf(vtk_name, "PNP-1bl8-c2-%4.2f.vtk", u_bulk);
            change_file_path(tmp_path, &vtk_name);
            phgPrintf("Creating \"%s\".\n\n", phgExportVTK(g, vtk_name, p[1], NULL));

            sprintf(vtk_name, "PNP-1bl8-c3-%4.2f.vtk", u_bulk);
            change_file_path(tmp_path, &vtk_name);
            phgPrintf("Creating \"%s\".\n\n", phgExportVTK(g, vtk_name, p[2], NULL));
*/
	}

	for(i = 0; i < N; i++) {
		for(which_ion = 0; which_ion < NION; which_ion++) {
			P[which_ion]->data[i] /= Pow(change_mesh_unit, 3) / Pow(change_unit, 3) / NA;
			p[which_ion]->data[i] /= Pow(change_mesh_unit, 3) / Pow(change_unit, 3) / NA;
		}
	}  

	phgMemoryUsage(g, &mem_peak);
	FLOAT t1 = phgGetTime(NULL);
	phgPrintf("\n Walltime = %.4e s    mem = %0.2lf MB\n\n", (t1 - t0), (double)mem_peak / (1024.0 * 1024.0));

	phgDofFree(&bdry);
	phgDofFree(&epsl);
	phgDofFree(&ph);
	phgDofFree(&ph_u_D);
	phgDofFree(&u);
	phgDofFree(&tu);
	phgDofFree(&tmp_u);
	phgDofFree(&TMP_u);
	phgDofFree(&u_D);
	phgDofFree(&tmp_u_D);
	phgDofFree(&TMP_u_D);
	phgDofFree(&G);
	phgDofFree(&H);
	phgDofFree(&grad_G);
	phgDofFree(&grad_H);
	phgDofFree(&chp0);
	phgDofFree(&chp1);
	phgDofFree(&Psi); //rgshen2021-4-27
	for(which_ion = 0; which_ion < NION; which_ion++) {
		phgDofFree(&TMP_p[which_ion]);
		phgDofFree(&TMP_P[which_ion]);
		phgDofFree(&p[which_ion]);
		phgDofFree(&P[which_ion]);
		phgDofFree(&size_c[which_ion]); //rgshen2021-6-11
		if(V_D) {
			phgDofFree(&D[which_ion]);
		}
	}
        phgDofFree(&Epsilon);
	phgDofFree(&size_p);
	phgDofFree(&indicator);
	phgDofFree(&error);
	phgFreeGrid(&g);
	phgFinalize();

	return 0;
}



