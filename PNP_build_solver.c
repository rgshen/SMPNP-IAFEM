#include"phg.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#include"PNP.h"  //rgshen

#include"PNP_func.h"
#include"PNP_quad.h"




static void Build_H_Solver(SOLVER *solver, DOF *H) {
	int i, j;
	DOF *HH;
	GRID *g = H->g;
	SIMPLEX *e;
	HH = phgDofCopy(H, NULL, NULL, NULL);

	int N = HH->type->nbas;
	FLOAT **A, *buffer, *rhs;
	INT *I;
	A = (FLOAT **)malloc(N * sizeof(FLOAT *));
	for(i = 0; i < N; i++) {
		A[i] = (FLOAT *)malloc(N * sizeof(FLOAT));
	}
	rhs = (FLOAT *)malloc(N * sizeof(FLOAT));
	buffer = (FLOAT *)malloc(N * sizeof(FLOAT));
	I = (INT *)malloc(N * sizeof(INT));

	ForAllElements(g, e) {
		bzero(rhs, N * sizeof(*rhs));

		for(i = 0; i < N; i++) {
			I[i] = phgSolverMapE2L(solver, 0, e, i);
		}

//		if(e->region_mark == solute_region) {
                if(e->region_mark == solute_region || e->region_mark == membrane_region) {
			for(i = 0; i < N; i++) {
				for(j = 0; j <= i; j++)
					A[j][i] = A[i][j] = phgQuadGradBasDotGradBas(e, H, i, H, j, QUAD_DEFAULT);
			}
		}
		else if(e->region_mark == solvent_region) {
			for(i = 0; i < N; i++) {
				bzero(A[i], N * sizeof(**A));
			}
			for(i = 0; i < N; i++) {
				if(phgDofGetElementBoundaryType(H, e, i) & bmlc_bdry)
					A[i][i] = 0.0;
				else
					A[i][i] = 1.0;
			}
		}
		phgDofSetDirichletBoundaryMask(HH, bmlc_bdry | c_bdry);
		//phgDofSetDirichletBoundaryMask(H, bmlc_bdry | c_bdry);
		for(i = 0; i < N; i++) {
			BOOLEAN is_bdry = FALSE;
//			if(e->region_mark == solute_region) {
            if(e->region_mark == solute_region || e->region_mark == membrane_region) {
				is_bdry = phgDofDirichletBC(HH, e, i, func_G, buffer, rhs + i, DOF_PROJ_NONE);
				if(is_bdry)
					*(rhs + i) *= -1.0;
			}
			if(is_bdry)
				phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer);
			else
				phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
	phgDofFree(&HH);
	for(i = 0; i < N; i++) {
		free(A[i]);
	}
	free(rhs);
	free(buffer);
	free(I);
	free(A);
}





/*----------------------------------------------------------------------*/
//-----求解 pnp 方程  新归纳在一起的---2021-6-14  -----重新放了位置
/*
	matrix		u	p1	p2	p3		rhs
	v		.	.	.	.		u	.
	ph1		.	.	.	.		p1	.
	ph2		.	.	.	.		p2	.
	ph3		.	.	.	.		p3	.
*/

//Primitive ------------------- Poisson solver ---------------
//This solver always needs small coefficient of relaxation at the beginning of iteration
static void Build_Primitive_Solver(SOLVER* solver, DOF *u, DOF *size_p, DOF **p, DOF *grad_G, DOF *grad_H) {
	GRID *g = u->g;
	SIMPLEX *e;
	int n;
	int i, j;
	int N = u->type->nbas;
	FLOAT **A, *buffer, *rhs;
	INT *I;
	DOF *E;
	if(surface_charge) {
		if(V_scd) {
			E = phgDofNew(g, DOF_DEFAULT, 1, "scd", func_scd);
		}
		else {
			E = phgDofNew(g, DOF_CONSTANT, 1, "unit", DofNoAction);
			phgDofSetDataByValue(E, 1.0);
		}
	}

	A = (FLOAT **)malloc(N * sizeof(FLOAT *));
	for(i = 0; i < N; i++) {
		A[i] = (FLOAT *)malloc(N * sizeof(FLOAT));
	}
	rhs = (FLOAT *)malloc(N * sizeof(FLOAT));
	buffer = (FLOAT *)malloc(N * sizeof(FLOAT));
	I = (INT *)malloc(N * sizeof(INT));

	FLOAT Epsilon = 0.0;
	FLOAT Kappa = 0.0;
	FLOAT value = 0.0;

	ForAllElements(g, e) {
		for(i = 0; i < N; i++) {
			memset(A[i], 0, N * sizeof(**A));
		}
		memset(buffer, 0, N * sizeof(*buffer));
		memset(rhs, 0, N * sizeof(*rhs));

		get_Kappa(e, &Epsilon, &Kappa);

		for(i = 0; i < N; i++) {
			I[i] = phgSolverMapE2L(solver, 0, e, i);
		}
		//matrix
		for(i = 0; i < N; i++) {
			for(j = 0; j <= i; j++) {
				A[i][j] = A[j][i] = Epsilon * phgQuadGradBasDotGradBas(e, u, j, u, i, 3);
			}
		}
		//rhs
		for(i = 0; i < N; i++) {
			for(j = 0; j < NBCMAP; j++) {
				phgDofSetDirichletBoundaryMask(u, BDRY_USER0 * pow(2, j));
                if(voltage_test) {
                    if(bc_tab[j] == 1) {
                    	if(phgDofDirichletBC(u, e, i, func_ub1, buffer, rhs + i, DOF_PROJ_NONE)) {
                    		phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
                    		break;
                    	}
                    }
                    else if(bc_tab[j] == 0 || bc_tab[j] == 2 || bc_tab[j] == 5) {
                    	if(phgDofDirichletBC(u, e, i, func_ub, buffer, rhs + i, DOF_PROJ_NONE)) {
                    		phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
                    		break;
                    	}
                    }
                }
                else{
                    if(bc_tab[j] == 1) {
                            if(phgDofDirichletBC(u, e, i, NULL, buffer, rhs + i, DOF_PROJ_NONE)) {
                                    phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
                                    break;
                            }
                    }
                    else if(bc_tab[j] == 0 || bc_tab[j] == 2 || bc_tab[j] == 5) {
                            if(phgDofDirichletBC(u, e, i, func_ub, buffer, rhs + i, DOF_PROJ_NONE)) {
                                    phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
                                    break;
                            }
                    }
                }
			}
			if(j == NBCMAP) {
				phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]);
				for(n = 0; n < NION; n++) {
					phgQuadDofTimesBas(e, p[n], u, i, 3, &value);
					    rhs[i] += Kappa * ion[n].Z * value / Pow(change_unit, 2);  // 这是LIUXUEJIAO无量纲方式
					       //\rhs[i] +=  Kappa * NA / Pow(10, -3) * Pow(10, -20) * ion[n].Z * value ;  // rgshen-->直接始终用 mol/L， 用这条
				}
//				if(e->region_mark == solute_region && !surface_charge) {
                if((e->region_mark == solute_region || e->region_mark == membrane_region) && !surface_charge) {
					for(j = 0; j < NFace; j++) {
						if(e->bound_type[j] & bmlc_bdry) {
							rhs[i] -= grow_Lambda * Epsilon_m * phgQuadFaceDofDotBas(e, j,
								grad_G, DOF_PROJ_DOT, u, i, QUAD_DEFAULT);
							rhs[i] -= grow_Lambda * Epsilon_m * phgQuadFaceDofDotBas(e, j,
								grad_H, DOF_PROJ_DOT, u, i, QUAD_DEFAULT);
						}
					}
				}
				if(e->region_mark == solvent_region && surface_charge) {
					for(j = 0; j < NFace; j++) {
						if(e->bound_type[j] & scd_bdry) {
							//rhs[i] += scd * Kappa * change_unit * phgQuadFaceDofDotBas(e, j,
							if(SOU){
								if(V_scd) {
									rhs[i] += Beta * ec / Epsilon_vacuum / change_unit
										* phgQuadFaceDofDotBas(e, j, E, DOF_PROJ_NONE, u, i, QUAD_DEFAULT);
								}
								else {
									rhs[i] += scd * Beta * ec / Epsilon_vacuum / change_unit
										* phgQuadFaceDofDotBas(e, j, E, DOF_PROJ_NONE, u, i, QUAD_DEFAULT);
								}
							}
							else{
								rhs[i] += scd * Beta * ec / change_unit
									* phgQuadFaceDofDotBas(e, j, E, DOF_PROJ_NONE, u, i, QUAD_DEFAULT);
							}
						}
					}
				}
			}
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}

	if(surface_charge) phgDofFree(&E);
	for(i = 0; i < N; i++) {
		free(A[i]);
	}
	free(rhs);
	free(buffer);
	free(I);
	free(A);
	phgDofSetDirichletBoundaryMask(u, DIRICHLET);
}


//----需要小松弛参数
//ok  ----------这里面包括 IAFE 方法，written in rgshen2021-4-19
          //Build_p_Solver(        solver,      u,      size_p,   TMP_p,      D[i]);  对应调用位置
static void Build_p_Solver(SOLVER *solver, DOF *u, DOF *size_p, DOF **p, DOF *D) {
	GRID *g = u->g;
	SIMPLEX *e;
        DOF *grad_u = phgDofGradient(u, NULL, NULL, NULL);

        DOF *TMP, *tmp;
        tmp = phgDofNew(g, DOF_DEFAULT, 1, "tmp", DofInterpolation);
        phgDofSetDataByValue(tmp, 0.0);

	INT i, j;

        if(MOD == 0) {
                TMP = phgDofGradient(tmp, NULL, NULL, NULL);//TMP = 0
        }

        phgDofFree(&tmp);

	int N = p[which_ion]->type->nbas;
	FLOAT **A, *buffer, *rhs, value;
	INT *I;

    FLOAT psi0, psi1, B, BB; //rgshen
	A = (FLOAT **)malloc(N * sizeof(FLOAT *));
	for(i = 0; i < N; i++) {
		A[i] = (FLOAT *)malloc(N * sizeof(FLOAT));
	}
	rhs = (FLOAT *)malloc(N * sizeof(FLOAT));
	buffer = (FLOAT *)malloc(N * sizeof(FLOAT));
	I = (INT *)malloc(N * sizeof(INT));

	ForAllElements(g, e) {
		for(i = 0; i < N; i++) {
			memset(A[i], 0, N * sizeof(**A));
		}
		memset(buffer, 0, N * sizeof(*buffer));
		memset(rhs, 0, N * sizeof(*rhs));
		for(i = 0; i < N; i++) {
			I[i] = phgSolverMapE2L(solver, 0, e, i);
		}

//		if(e->region_mark == solute_region) {
        if(e->region_mark == solute_region || e->region_mark == membrane_region) {
			for(i = 0; i < N; i++) {
				if(phgDofGetElementBoundaryType(p[which_ion], e, i) & bmlc_bdry) {
					A[i][i] = 0;
				}
				else {
					A[i][i] = 1;
				}
			}
		}
		else if(e->region_mark == solvent_region) {
            if(Stabilize == EAFE) { // EAFE = 3, from zhangqr, written in rgshen2021-4-19
				for(i = 0; i < N; i++) {
                    psi0 = *(DofVertexData(u, e->verts[i]));
                    A[i][i] = 0.;
                    for(j = 0; j < N; j++) {
                        if(j != i) {
                            psi1 = *(DofVertexData(u, e->verts[j]));
                         //----------------------------------------------始终打开，若只算俩离子情形，就走这里--------------
                            if (NION == 2) {
							     B = which_ion == 0 ? Ber(psi0 - psi1) : Ber(psi1 - psi0);
                                 BB = which_ion == 0 ? Ber(psi1 - psi0) : Ber(psi0 - psi1);     //这是原来只能算俩离子(正、负)情形
							}
						 //-------------------------------------------------------------------------------------------------
					     //-------为了能算三离子改的（添加的）- 若只算两个方程，就只走上面两条即可---2021-6-3-------------------
                            if (NION == 3) {
					        	if(which_ion == 0 || which_ion ==1) {
					        		B = Ber(psi0 - psi1);
                                    BB = Ber(psi1 - psi0);   //这两条执行正离子方程，如 Na， K实则是相同的对pnp
					        	}
                                else{
					        	    B = Ber(psi1 - psi0);
                                    BB = Ber(psi0 - psi1);   //这条执行负离子方程
					        	}
                            }
					     //----------------------------------------------
                            if(V_D==2) { //若是变扩散系数，走这支，gA 通道是变系数问题
                                A[i][j] = B * phgQuadGradBasAGradBas(e, p[which_ion], i, D, p[which_ion], j, QUAD_DEFAULT);
                                A[i][i] -= BB * phgQuadGradBasAGradBas(e, p[which_ion], i, D, p[which_ion], j, QUAD_DEFAULT);
                             }
                             else{
		    			        A[i][j] = B * ion[which_ion].D * phgQuadGradBasDotGradBas(e, p[which_ion], i, p[which_ion], j, QUAD_DEFAULT);
		    			        A[i][i] -= BB * ion[which_ion].D * phgQuadGradBasDotGradBas(e, p[which_ion], i, p[which_ion], j, QUAD_DEFAULT);
                             }
                        }
                    }
                }
            }
		    else { //rgshen --下面是标准有限元法 FEM
			    for(i = 0; i < N; i++) {
			    	for(j = 0; j < N; j++) {
			    		/*对称部分:( grad p, grad v)*/
			    		A[i][j] = ion[which_ion].D * phgQuadGradBasDotGradBas(e, p[which_ion], i, p[which_ion], j, QUAD_DEFAULT);
			    		/*非线性部分 : (p*grad_u, grad v)*/
			    		A[i][j] += ion[which_ion].D * ion[which_ion].Z * PNP_Quad_1_D_Bas_G_D_G_B(e, func_unit, p[which_ion],
			    				        p[which_ion], j, grad_u, p[which_ion], i, QUAD_DEFAULT);
			    		//-----------------------------------------------------------------------------------------------------
			    		if(Stabilize == SUPG) {  // 若用SUPG = 1 稳定化，再加上该项 shenruig
			    			Region_Quad_SUPG_term_Primitive(2, e, grad_u, p[which_ion], j, i, 4, &value);
			    			A[i][j] += value;
			    		}
			    	}
			    }
		    } //rgshen
	    }

		for(i = 0; i < N; i++) {
			for(j = 0; j < NBCMAP; j++) {
				phgDofSetDirichletBoundaryMask(p[which_ion], BDRY_USER0 * pow(2, j));
				if(bc_tab[j] == 2) {
					if(phgDofDirichletBC(p[which_ion], e, i, NULL, buffer, rhs + i, DOF_PROJ_NONE)) {
						phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
						break;
					}
				}
				else if(bc_tab[j] == 0 || bc_tab[j] == 1) {
					if(phgDofDirichletBC(p[which_ion], e, i, func_pb, buffer, rhs + i, DOF_PROJ_NONE)) {
						phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
						break;
					}
				}
			}
			if(j == NBCMAP) phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
	phgDofFree(&TMP);
	phgDofFree(&grad_u);
	for(i = 0; i < N; i++) {
		free(A[i]);
	}
	free(rhs);
	free(buffer);
	free(I);
	free(A);
	phgDofSetDirichletBoundaryMask(p[which_ion], DIRICHLET);
}

//--------上面是原始书写 继承的



//000000000000000000000000
void get_size_u_p(DOF **tmp, DOF *u, DOF **p) {
    INT i, N;
    int j;
    N = DofGetDataCount(*tmp);            //a0 = 11.5; /* size of H2O */ ---已经在内部初始定义
	                                      //K_i[which_ion] = Pow(ion[which_ion].a, 3) / Pow(a0, 3)
    for(i = 0; i < N; i++) {
		(*tmp)->data[i] = 1.0;
        for(j = 0; j < NION; j++) {
			(*tmp)->data[i] -= Pow(ion[j].a, 3) * p[j]->data[i] / Pow(change_unit * mesh_unit, 3);
			//(*tmp)->data[i] -= asize_ratio * Pow(ion[j].a, 3) * p[j]->data[i];  //  直接 mol/L 传入时用这个 
        }
		if( (*tmp)->data[i] < -1e-8 ) (*tmp)->data[i] = - ( (*tmp)->data[i] );
		(*tmp)->data[i] = u->data[i] - K_i[which_ion] * log( (*tmp)->data[i] ); //Psi = u - ki* ln(1 - Sigma(al^3 * pl)
		//phgPrintf("-----%f", (*tmp)->data[i]);
	}                               //只要"u->data[i]"就是标准pnp(或许尺寸取0，也退化成pnp) //因为下面的EAFE 方法已经考虑了指数前面的正负号，所以u前面就是正号 ？
}                                 


/* --通过，-2021-4-20 SMPNP -平均技术法-------这里面包括 EAFE 方法，form张倩茹  written in rgshen */
	      //SMPNP_NP_IAFEM(        solver,     Psi,    TMP_p,   D[i]);
static void SMPNP_NP_IAFEM(SOLVER *solver, DOF *u,  DOF **p, DOF *D) {
	GRID *g = u->g;
	SIMPLEX *e;

	INT i, j;

	int N = p[which_ion]->type->nbas;
	FLOAT **A, *buffer, *rhs, value;
	INT *I;

    FLOAT psi0, psi1, B, BB; //rgshen

	A = (FLOAT **)malloc(N * sizeof(FLOAT *));
	for(i = 0; i < N; i++) {
		A[i] = (FLOAT *)malloc(N * sizeof(FLOAT));
	}
	rhs = (FLOAT *)malloc(N * sizeof(FLOAT));
	buffer = (FLOAT *)malloc(N * sizeof(FLOAT));
	I = (INT *)malloc(N * sizeof(INT));

	ForAllElements(g, e) {
		for(i = 0; i < N; i++) {
			memset(A[i], 0, N * sizeof(**A));
		}
		memset(buffer, 0, N * sizeof(*buffer));
		memset(rhs, 0, N * sizeof(*rhs));
		for(i = 0; i < N; i++) {
			I[i] = phgSolverMapE2L(solver, 0, e, i);
		}

//		if(e->region_mark == solute_region) {
        if(e->region_mark == solute_region || e->region_mark == membrane_region) {
			for(i = 0; i < N; i++) {
				if(phgDofGetElementBoundaryType(p[which_ion], e, i) & bmlc_bdry) {
					A[i][i] = 0;
				}
				else {
					A[i][i] = 1;
				}
			}
		}
		else if(e->region_mark == solvent_region) {
            if(Stabilize == EAFE) { // EAFE , from zhangqiaru, written in rgshen2021-4-19
                for(i = 0; i < N; i++) {
                    psi0 = *(DofVertexData(u, e->verts[i])); // 这里的u = zi*u - ki* ln((1 - Sigma(al^3 * pl)))
                    A[i][i] = 0.;
                    for(j = 0; j < N; j++) {
                        if(j != i) {
                            psi1 = *(DofVertexData(u, e->verts[j]));
						//----------------------------------------------始终打开，若只算俩离子情形，就走这里--------------
                            if (NION == 2) {
							    B = which_ion == 0 ? Ber(psi0 - psi1) : Ber(psi1 - psi0);
                                BB = which_ion == 0 ? Ber(psi1 - psi0) : Ber(psi0 - psi1);     //这是原来只能算俩离子(正、负)情形
							}
						//-------------------------------------------------------------------------------------------------
					    //-------为了能算三离子改的（添加的）- 若只算两个方程，就只走上面两条即可---2021-6-3-------------------
                            if (NION == 3) {
					        	if(which_ion == 0 || which_ion ==1) {
					        		B = Ber(psi0 - psi1);
                                    BB = Ber(psi1 - psi0);   //这两条执行正离子方程，如 Na， K实则是相同的对pnp
					        	}
                                else{
					        	    B = Ber(psi1 - psi0);
                                    BB = Ber(psi0 - psi1);   //这条执行负离子方程
					        	}
                            }
					    //----------------------------------------------
                            if(V_D == 2) { //若是变扩散系数，走这支，gA 通道是变系数问题
                                //phgPrintf("i am here -6-19--11 \n");
                                A[i][j] = B * phgQuadGradBasAGradBas(e, p[which_ion], i, D, p[which_ion], j, QUAD_DEFAULT);
                                A[i][i] -= BB * phgQuadGradBasAGradBas(e, p[which_ion], i, D, p[which_ion], j, QUAD_DEFAULT);
                            }
                            else{
		    			        A[i][j] = B * ion[which_ion].D * phgQuadGradBasDotGradBas(e, p[which_ion], i, p[which_ion], j, QUAD_DEFAULT);
		    			        A[i][i] -= BB * ion[which_ion].D * phgQuadGradBasDotGradBas(e, p[which_ion], i, p[which_ion], j, QUAD_DEFAULT);
                            }
                        }
                    }
                }
            }
	    }

		for(i = 0; i < N; i++) {
			for(j = 0; j < NBCMAP; j++) {
				phgDofSetDirichletBoundaryMask(p[which_ion], BDRY_USER0 * pow(2, j));
				if(bc_tab[j] == 2) {
					if(phgDofDirichletBC(p[which_ion], e, i, NULL, buffer, rhs + i, DOF_PROJ_NONE)) {
						phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
						break;
					}
				}
				else if(bc_tab[j] == 0 || bc_tab[j] == 1) {
					if(phgDofDirichletBC(p[which_ion], e, i, func_pb, buffer, rhs + i, DOF_PROJ_NONE)) {
						phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
						break;
					}
				}
			}
			if(j == NBCMAP) phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}

	for(i = 0; i < N; i++) {
		free(A[i]);
	}
	free(rhs);
	free(buffer);
	free(I);
	free(A);
	phgDofSetDirichletBoundaryMask(p[which_ion], DIRICHLET);
}






/* ---SMPNP primitive method 2021-6-11----通过-----------------------*/
static void get_size_c(DOF **tmp, DOF **p) {
    INT i, N;
    int j;
    N = DofGetDataCount(*tmp);            //a0 = 3.1; /* size of H2O */ ---已经在内部初始定义
	                                      //K_i[which_ion] = Pow(ion[which_ion].a, 3) / Pow(a0, 3)
    for(i = 0; i < N; i++) {
		(*tmp)->data[i] = 1.0;
        for(j = 0; j < NION; j++) {
			(*tmp)->data[i] -= Pow(ion[j].a, 3) * p[j]->data[i] / Pow(change_unit * mesh_unit, 3);  //liuxuejiao 无量纲方式
			//(*tmp)->data[i] -= asize_ratio * Pow(ion[j].a, 3) * p[j]->data[i];  //--> (1 - Sigma(al^3 * pl))
        }
		if( (*tmp)->data[i] < -1e-8 ) (*tmp)->data[i] = - ( (*tmp)->data[i] );
		(*tmp)->data[i] = K_i[which_ion] / (*tmp)->data[i];  //--> k_i /(1 - Sigma(al^3 * pl))
	}
}


/* ---SMPNP primitive method 2021-6-11----通过-----目前最多只能实现三个离子情形！！！6-16 0:20-------*/
static void SMPNP_NP_Slover_Primitive(SOLVER *solver, DOF *u, DOF **size_c, DOF **p, DOF *D, DOF *size_p) {
	GRID *g = u->g;
	SIMPLEX *e;
        DOF *grad_u = phgDofGradient(u, NULL, NULL, NULL);
//--------------------------------------------------------------------------------------------------
	    DOF *grad_p1 = phgDofGradient(p[0], NULL, NULL, NULL);
		DOF *grad_p2 = phgDofGradient(p[1], NULL, NULL, NULL);


	phgDofAXPBY(Pow(ion[0].a, 3), grad_p1, Pow(ion[1].a, 3), &grad_p2);  // \sum a_l^3 grad c_l,  l=1,2
	phgDofAXPBY(0.0, grad_p1, 1.0 / Pow(change_unit * mesh_unit, 3), &grad_p2);
	//\phgDofAXPBY(0.0, grad_p1, asize_ratio, &grad_p2);  //进一步乘以无量纲比率系数

	DOF *grad_p1p2 = phgDofCopy(grad_p2, NULL, NULL, NULL);

	phgDofFree(&grad_p1);
	phgDofFree(&grad_p2);
//--------------------------------------------------------------------------------------------------
	DOF *TMP;
	if(NION == 3) { // 若是三个离子情形，则再加上第三个离子的梯度
	    DOF *grad_p3 = phgDofGradient(p[2], NULL, NULL, NULL);
		phgDofAXPBY(1.0, grad_p1p2, 1.0 / Pow(change_unit * mesh_unit, 3) * Pow(ion[2].a, 3), &grad_p3);
		//\phgDofAXPBY(1.0, grad_p1p2, asize_ratio * Pow(ion[2].a, 3), &grad_p3);  // \sum a_l^3 grad c_l, l=1,3
		TMP = phgDofCopy(grad_p3, NULL, NULL, NULL);
	    phgDofFree(&grad_p3);
	}
//---------------------------------------------------------------------------------------------------
//------新写法-循环离子求和型（实现任意多个离子情形）2021-6-16, 21:24-----经验证，跟我上述使用 size_c 传入差不多，就用上一种了
    //\\DOF *TMP1, *tmp1;
    //\\tmp1 = phgDofNew(g, DOF_DEFAULT, 1, "tmp", DofInterpolation);
    //\\phgDofSetDataByValue(tmp1, 0.0);
    //\\INT n;
	//\\for(n = 0; n < NION; n++) {
    //\\    phgDofAXPY(K_i[which_ion] * asize_ratio * Pow(ion[n].a, 3), p[n], &tmp1);
    //\\}
    //\\TMP1 = phgDofGradient(tmp1, NULL, NULL, NULL);//TMP = Sigma(al^3 * grad_pl)
	//\\
    //\\dof_pow(&tmp1, size_p, -1);
    //\\dof_coef_vec(&TMP1, tmp1, TMP1);//TMP = Sigma(al^3 * grad_pl) / (1 - Sigma(al^3 * pl))
    //\\phgDofSetDataByValue(tmp1, 0.0);
//-------------------------------------


	INT i, j;

	int N = p[which_ion]->type->nbas;
	FLOAT **A, *buffer, *rhs, value;
	INT *I;

	A = (FLOAT **)malloc(N * sizeof(FLOAT *));
	for(i = 0; i < N; i++) {
		A[i] = (FLOAT *)malloc(N * sizeof(FLOAT));
	}
	rhs = (FLOAT *)malloc(N * sizeof(FLOAT));
	buffer = (FLOAT *)malloc(N * sizeof(FLOAT));
	I = (INT *)malloc(N * sizeof(INT));

	ForAllElements(g, e) {
		for(i = 0; i < N; i++) {
			memset(A[i], 0, N * sizeof(**A));
		}
		memset(buffer, 0, N * sizeof(*buffer));
		memset(rhs, 0, N * sizeof(*rhs));
		for(i = 0; i < N; i++) {
			I[i] = phgSolverMapE2L(solver, 0, e, i);
		}

//		if(e->region_mark == solute_region) {
        if(e->region_mark == solute_region || e->region_mark == membrane_region) {
			for(i = 0; i < N; i++) {
				if(phgDofGetElementBoundaryType(p[which_ion], e, i) & bmlc_bdry) {
					A[i][i] = 0;
				}
				else {
					A[i][i] = 1;
				}
			}
		}
		else if(e->region_mark == solvent_region) {
		    for(i = 0; i < N; i++) {
		    	for(j = 0; j < N; j++) {
		    		/*对称部分:( grad p, grad v)*/
		    		A[i][j] = ion[which_ion].D * phgQuadGradBasDotGradBas(e, p[which_ion], i, p[which_ion], j, QUAD_DEFAULT);

		    		/*非线性部分 : Z_i(p*grad_u, grad v)*/
		    		A[i][j] += ion[which_ion].D * ion[which_ion].Z * PNP_Quad_1_D_Bas_G_D_G_B(e, func_unit, p[which_ion],
		    				   p[which_ion], j, grad_u, p[which_ion], i, QUAD_DEFAULT);

		    		/* smpnp 尺寸效应部分：(p* k_i / (1 - \sum a_l^3 p_l) (\sum a_l^3 \grad p_l), grad v)*/
                    if(NION == 2) {
		    	        A[i][j] += ion[which_ion].D * PNP_Quad_2_D_Bas_G_D_G_B(e, func_unit, p[which_ion], size_c[which_ion],
		    	                   p[which_ion], j, grad_p1p2, p[which_ion], i, QUAD_DEFAULT);

						/* A[i][j] += ion[which_ion].D * PNP_Quad_1_D_Bas_G_D_G_B(e, func_unit, p[which_ion],
		    	                   p[which_ion], j, TMP1, p[which_ion], i, QUAD_DEFAULT); */ //经过验证，跟我上述使用差不多
		    		}
					else if(NION == 3) {
                        A[i][j] += ion[which_ion].D * PNP_Quad_2_D_Bas_G_D_G_B(e, func_unit, p[which_ion], size_c[which_ion],
		    	                   p[which_ion], j, TMP, p[which_ion], i, QUAD_DEFAULT);
					}
		    	}
		    }
		}


		for(i = 0; i < N; i++) {
			for(j = 0; j < NBCMAP; j++) {
				phgDofSetDirichletBoundaryMask(p[which_ion], BDRY_USER0 * pow(2, j));
				if(bc_tab[j] == 2) {
					if(phgDofDirichletBC(p[which_ion], e, i, NULL, buffer, rhs + i, DOF_PROJ_NONE)) {
						phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
						break;
					}
				}
				else if(bc_tab[j] == 0 || bc_tab[j] == 1) {
					if(phgDofDirichletBC(p[which_ion], e, i, func_pb, buffer, rhs + i, DOF_PROJ_NONE)) {
						phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
						break;
					}
				}
			}
			if(j == NBCMAP) phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}

	phgDofFree(&grad_u);
	phgDofFree(&grad_p1p2);
	if(NION == 3) {
	    phgDofFree(&TMP);
	}
	for(i = 0; i < N; i++) {
		free(A[i]);
	}
	free(rhs);
	free(buffer);
	free(I);
	free(A);
	phgDofSetDirichletBoundaryMask(p[which_ion], DIRICHLET);
}



//------下面是来源乔瑜  标准 FEM 求解 smpnp 方程的2021-7-1, 但是边界信息改成 本套code 匹配的--通过
//static void myquad(QUAD *quad, SIMPLEX *e, int i, int j, DOF *p, FLOAT *DATA2, FLOAT *DATA5, FLOAT **DATAc, FLOAT **DATAGc, FLOAT h, FLOAT *value){
static void smpnp_myquad(QUAD *quad, SIMPLEX *e, int i, int j, DOF *p, FLOAT *DATA2, FLOAT **DATAc, FLOAT **DATAGc, FLOAT h, FLOAT *value){	
	
	FLOAT sum1, sum2, sum3, res = 0.0;
	int I, J, k;
	FLOAT *DATA1, *DATA3, *DATA4;
	FLOAT tm, tmp;
	
	DATA1 = phgQuadGetBasisValues(e, p, j, quad);
	DATA3 = phgQuadGetBasisGradient(e, p, i, quad);
	DATA4 = phgQuadGetBasisGradient(e, p, j, quad);

	I = quad->npoints;
	
	for( J = 0; J < I; J++ ){
		sum1 = DATA4[3*J+0]*DATA3[3*J+0] + DATA4[3*J+1]*DATA3[3*J+1] + DATA4[3*J+2]*DATA3[3*J+2];
		sum2 = ion[which_ion].Z * DATA1[J] * (DATA2[3*J+0]*DATA3[3*J+0] + DATA2[3*J+1]*DATA3[3*J+1] + DATA2[3*J+2]*DATA3[3*J+2]);
		tm = 0.0;
		for(k = 0; k < NION; k++){
			tm += pow(ion[k].a, 3) * DATAc[k][J];
		}
		tmp = 1 - tm * NA *1e-27;
		if(tmp < -1e-8) tmp = -tmp;
		tm = 0.0;
		for(k = 0; k < NION; k++){
			tm += pow(ion[k].a, 3) * (DATAGc[k][3*J+0]*DATA3[3*J+0]+DATAGc[k][3*J+1]*DATA3[3*J+1]+DATAGc[k][3*J+2]*DATA3[3*J+2]);
		}
		sum3 = K_i[which_ion] * DATA1[J] * tm * NA * 1e-27/tmp;

		//res += DATA5[J] * (sum1 + sum2 + sum3) * quad->weights[J];  //变扩散系数情形，如gA 通道时
		res += ion[which_ion].D * (sum1 + sum2 + sum3) * quad->weights[J];
	}
	*value = res * h;
}


//build linear system for the NP equation
static void SMPNP_PRIMITIVE_NP(SOLVER *solver, DOF *u, DOF *p, DOF **c, DOF *d){
	GRID *g = u->g;
	SIMPLEX *e;
	QUAD *quad;
	quad = phgQuadGetQuad3D(3);
	DOF *Grad_u, *pp, **Grad_c;
	Grad_u  = phgDofGradient(u,  NULL, NULL, NULL);
	pp = phgDofCopy(p, NULL, NULL, NULL);
	if(!(Grad_c = phgAlloc(NION * sizeof(*Grad_c))))
		phgError(1, "Error in memory allocation for gradient of concentration c.\n");
	int i, j;
	for(i = 0; i < NION; i++)
		Grad_c[i] = phgDofGradient(c[i], NULL, NULL, NULL);
	FLOAT value;

	ForAllElements(g,e){
		int N = DofGetNBas(p, e);
		FLOAT A[N][N], rhs[N];
		FLOAT buffer[N];
		int I[N];

		memset(A, 0, N*N*sizeof(**A));
		memset(rhs, 0, N*sizeof(*rhs));
       
		for( i = 0; i < N; i++ )
			I[i] = phgSolverMapE2L(solver, 0, e, i);
          
     // if(e->region_mark == solute_region) {
        if(e->region_mark == solute_region || e->region_mark == membrane_region) {
			for( i = 0; i < N; i++ ){
				if(phgDofGetElementBoundaryType(p, e, i) & bmlc_bdry)
 					A[i][i] = 0;
				else A[i][i] = 1;
			}
		}
		else {	
			FLOAT *DATA2, *DATA5, h, **DATAc, **DATAGc;
			h = phgGeomGetVolume(g, e);
			DATA2 = phgQuadGetDofValues(e, Grad_u, quad);
			
			//\DATA5 = phgQuadGetDofValues(e, d, quad); //本套代码 要算 gA 通道时才打开他，并修改想要位置
		
			if(!(DATAc = phgAlloc(NION * sizeof(*DATAc))))
				phgError(1, "Error in memory allocation for c.\n");
			if(!(DATAGc = phgAlloc(NION * sizeof(*DATAGc))))
				phgError(1, "Error in memory allocation for Grad_c.\n");
			for(i = 0; i < NION; i++){
				DATAc[i] = phgQuadGetDofValues(e, c[i], quad);
				DATAGc[i] = phgQuadGetDofValues(e, Grad_c[i], quad);
			}
			for(i = 0; i < N; i++){
				for(j = 0; j < N; j++){
					//myquad(quad, e, i, j, p, DATA2, DATA5, DATAc, DATAGc, h, &value); // 扩散系数 D 是变系数情形时候 才用它
					smpnp_myquad(quad, e, i, j, p, DATA2, DATAc, DATAGc, h, &value);  // 该函数里就是求数值积分
					A[i][j] = value;
				}    
			}     
			phgFree(DATAc);
			DATAc = NULL;  
			phgFree(DATAGc);
			DATAc = NULL;                   
		}
		          
        for(i = 0; i < N; i++) {
			for(j = 0; j < NBCMAP; j++) {
				phgDofSetDirichletBoundaryMask(p, BDRY_USER0 * pow(2, j));
				if(bc_tab[j] == 2) {
					if(phgDofDirichletBC(p, e, i, NULL, buffer, rhs + i, DOF_PROJ_NONE)) {
						phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
						break;
					}
				}
				else if(bc_tab[j] == 0 || bc_tab[j] == 1) {
					if(phgDofDirichletBC(p, e, i, func_pb, buffer, rhs + i, DOF_PROJ_NONE)) {
						phgSolverAddMatrixEntries(solver, 1, I + i , N, I, buffer);
						break;
					}
				}
			}
			if(j == NBCMAP) phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);		
	}
	
	phgDofFree(&Grad_u);   
	phgDofFree(&pp);
	for(i = 0; i < NION; i++){
		phgDofFree(Grad_c + i);
	}
}



