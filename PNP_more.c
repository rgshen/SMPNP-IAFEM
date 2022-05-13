#include"phg.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdarg.h>
#include<math.h>
#include"PNP_func.h"

void func_id(FLOAT *value) {
}

void func_exp_1(FLOAT *value) {
	*value = Exp(-ion[which_ion].Z * *value);
}

void func_exp_2(FLOAT *value) {
	*value = ion[which_ion].Z * Exp(-ion[which_ion].Z * *value);
}

void func_exp_3(FLOAT *value) {
	*value = Pow(ion[which_ion].Z, 2) * Exp(-ion[which_ion].Z * *value);
}

void func_exp_4(FLOAT *value) {
	*value = Exp(*value);
}

void func_unit(FLOAT *value) {
	*value = 1;
}
/* change measurement unit from N/cm^3 or N/m^3 -> N/A^3 -> mol/L for the formula */
//static void func_p_k(FLOAT *p, FLOAT *value) {
void func_p_k(FLOAT *p, FLOAT *value) {
        int i;
        *value = 0;
        for(i = 0; i < NION; i++) {
        	if (MOD == 2 && nonlinear_Epsilon == 2) {
                	*value += p[i];
		}
		else{
                	*value += Pow(ion[i].a, 3) * p[i];
		}
	}
        if (MOD == 2 && nonlinear_Epsilon == 2) {
		*value = 1 + A / NION * 1.0 / (Pow(10, 3) * NA) * *value;
	}
        else if (MOD == 2 && nonlinear_Epsilon == 1) {
                *value = 1 + 1.0 / Pow(change_unit * mesh_unit, 3) * *value;
        }
        else {
                *value = 1 - 1.0 / Pow(change_unit * mesh_unit, 3) * *value;
                if(*value < 1e-8) {
                        *value = 1e-8;
       		}
        }
}

/* set boundary of value of DOF *p with zero and DOF *P: 0 <= p <= P */
static void dof_positive(DOF *p) {
	INT i, N;
	N = DofGetDataCount(p);
	for(i = 0; i < N; i++) {
		if(p->data[i] < 0)
			//p->data[i] = 0;
			p->data[i] = 1e-8;
	}
}

static void dof_restriction(DOF *p, DOF *P) {
	INT i, N;
	N = DofGetDataCount(p);
	if(N != DofGetDataCount(P))
		phgError(1, "different dof used in dof_positive!!!\n");
	for(i = 0; i < N; i++) {
		if(p->data[i] > P->data[i]) p->data[i] = P->data[i];
	}
}

/* choose one direction of vector dof */
static void dof_direction(DOF **tmp, DOF *u, int d) {
	INT i, N;
	N = DofGetDataCount(*tmp);
	for(i = 0; i < N; i++) {
		(*tmp)->data[i] = u->data[3*i+d];
	}
}

static void dof_const(DOF **tmp, DOF *u, FLOAT k) {
	INT i, N;
	N = DofGetDataCount(*tmp);
	for(i = 0; i < N; i++) {
		(*tmp)->data[i] = u->data[i] * k;
	}
}

/*the DataCount has something wrong*/
static void dof_coef_vec(DOF **tmp, DOF *u, DOF *v) {
	int k;
	INT i, N;
	N = DofGetDataCount(u);
	for(i = 0; i < N; i++) {
		for(k = 0; k < 3; k++) {
			(*tmp)->data[3*i+k] = v->data[3*i+k] * u->data[i];
		}
	}
}

static void dof_pow(DOF **tmp, DOF *dof, FLOAT a) {
	INT i, N;
 	N = DofGetDataCount(*tmp);
	for(i = 0; i < N; i++) {
		if(dof->data[i] <= 0 && ((a - (int)a)!=0 || a <= 0)) {
			printf("\n error:%le can't be powered", dof->data[i]);
			exit(1);
		}
		(*tmp)->data[i] = Pow(dof->data[i], a);
	}
}

static void dof_plus(DOF **tmp, DOF *dof, ...) {
	int ndof;
	INT i, N;
	DOF **dofs;
	FLOAT res;
	va_list ap;
	dofs = phgAlloc(256 * sizeof(*dofs));
	N = DofGetDataCount(*tmp);
	DOF *redof = dof;
	for(i = 0; i < N; i++) {
		va_start(ap, dof);
		res = 0.0;
		for(ndof = 0; ndof < 256; ndof++) {
			if(dof == NULL)break;
			dofs[ndof] = dof;
			res += dofs[ndof]->data[i];
			dof = va_arg(ap, DOF *);
		}
		(*tmp)->data[i] = res;
		va_end(ap);
		dof = redof;
	}
	phgFree(dofs);
}

static void dof_multiple(DOF **tmp, DOF *dof, ...) {
	int ndof;
	INT i, N;
	DOF **dofs;
	FLOAT res;
	va_list ap;
	dofs = phgAlloc(256 * sizeof(*dofs));
	N = DofGetDataCount(*tmp);
	DOF *redof = dof;
	for(i = 0; i < N; i++) {
		res = 1.0;
		va_start(ap, dof);
		for(ndof = 0; ndof < 256; ndof++) {
			if(dof == NULL)break;
			dofs[ndof] = dof;
			res *= dofs[ndof]->data[i];
			dof = va_arg(ap, DOF *);
		}
		(*tmp)->data[i] = res;
		va_end(ap);
		dof = redof;
	}
	phgFree(dofs);
}

//static void Region_Quad_D_B_G_D_G_B(int region, SIMPLEX *e, DOF *dof1, DOF *u, int i, DOF *dof2, DOF *v, int j, int order, FLOAT *value){

static void Region_Quad_Dof_G_D_Bas_G_B(int region, SIMPLEX *e, DOF *dof, DOF *dof_vec, DOF *u, int i, DOF *v, int j, int order, FLOAT *value ) {
	INT I, J, K;
        QUAD *quad = phgQuadGetQuad3D(order);
	GRID *g = u->g;
	const FLOAT *data1, *data2, *data3, *data4;
	FLOAT res, h = phgGeomGetVolume(g, e);
	I = quad->npoints;

	data1 = phgQuadGetDofValues(e, dof, quad);
	data2 = phgQuadGetDofValues(e, dof_vec, quad);
	data3 = phgQuadGetBasisValues(e, u, i, quad);
	data4 = phgQuadGetBasisGradient(e, v, j, quad);

	*value = 0.0;

	if(e->region_mark == region || region < 0) {
		for(J = 0; J < I; J++) {
			res = 0.0;
			for(K = 0; K < 3; K++) {
				res += data2[3*J+K] * data4[3*J+K];
			}
			*value += res * data1[J] * data3[J] * h * (quad->weights[J]);
		}
	}
}


static void Region_Quad_Bas_Dof_G_B(int region, SIMPLEX *e, DOF *dof, DOF *u, int i, DOF *v, int j, int order, FLOAT *value) {
	INT I, J, K;
        QUAD *quad = phgQuadGetQuad3D(order);
	GRID *g = u->g;
	const FLOAT *data1, *data2, *data3;
	FLOAT res, h = phgGeomGetVolume(g, e);
	I = quad->npoints;

	data1 = phgQuadGetDofValues(e, dof, quad);
	data2 = phgQuadGetBasisValues(e, u, i, quad);
	data3 = phgQuadGetBasisGradient(e, v, j, quad);

	*value = 0.0;

	if(e->region_mark==region||region < 0) {
		for(J = 0; J < I; J++) {
			res = 0.0;
			for(K = 0; K < 3; K++) {
				res += data1[3*J+K] * data3[3*J+K];
			}
			*value += res * data2[J] * h * (quad->weights[J]);
		}
	}
}

static void Region_Quad_Dof_G_B_G_B(int region, SIMPLEX *e, DOF *dof, DOF *u, int i, DOF *v, int j, int order, FLOAT *value) {
	INT I, J, K;
        QUAD *quad = phgQuadGetQuad3D(order);
	GRID *g = u->g;
	const FLOAT *data1, *data2, *data3;
	FLOAT res, h = phgGeomGetVolume(g, e);
	I = quad->npoints;

	data1 = phgQuadGetDofValues(e, dof, quad);
	data2 = phgQuadGetBasisGradient(e, u, i, quad);
	data3 = phgQuadGetBasisGradient(e, v, j, quad);

	*value = 0.0;

	if(e->region_mark==region||region < 0) {
		for(J = 0; J < I; J++) {
			res = 0.0;
			for(K = 0; K < 3; K++) {
				res += data2[3*J+K] * data3[3*J+K];
			}
			*value += res * data1[J] * h * (quad->weights[J]);
		}
	}
}

static void Region_Quad_Dof_Bas_Bas(int region, SIMPLEX *e, DOF *dof, DOF *u, int i, DOF *v, int j, int order, FLOAT *value) {
	INT I, J;
        QUAD *quad = phgQuadGetQuad3D(order);
	GRID *g = u->g;
	const FLOAT *data1, *data2, *data3;
	FLOAT h = phgGeomGetVolume(g, e);
	I = quad->npoints;

	data1 = phgQuadGetDofValues(e, dof, quad);
	data2 = phgQuadGetBasisValues(e, u, i, quad);
	data3 = phgQuadGetBasisValues(e, v, j, quad);

	*value = 0.0;

	if(e->region_mark==region||region < 0) {
		for(J = 0; J < I; J++) {
			*value += data1[J] * data2[J] * data3[J] * h * (quad->weights[J]);
		}
	}
}

static void Region_Quad_SUPG_term_Primitive(int region, SIMPLEX *e, DOF *grad_u, DOF *p, int i, int j, int order, FLOAT *value) {
	INT I, J, K;
	QUAD *quad = phgQuadGetQuad3D(order);
	GRID *g = p->g;
	const FLOAT *DATA1, *DATA2, *DATA3, *DATA4;
	FLOAT h = phgGeomGetVolume(g, e);
	FLOAT d = phgGeomGetDiameter(g, e);
	DATA1 = phgQuadGetBasisValues(e, p, i, quad);
        DATA2 = phgQuadGetDofValues(e, grad_u, quad);
        DATA3 = phgQuadGetBasisGradient(e, p, j, quad);
        DATA4 = phgQuadGetBasisGradient(e, p, i, quad);
        I = quad->npoints;
	FLOAT Pe, Sigma, res = 0.0;
	FLOAT sum, a_L2 = 0.0;
	if(e->region_mark == region) {
		for(J = 0; J < I; J++) {
			for(K = 0; K < 3; K++) {
				a_L2 += (DATA2[3*J+K]*DATA2[3*J+K]);
			}
			a_L2 *= quad->weights[J];
		}
		if(a_L2 == 0)
			*value = 0;
		else {
			a_L2 = Sqrt(a_L2 * h) * ion[which_ion].D * ion[which_ion].Z;
			Pe = a_L2 * d / 6.0 / ion[which_ion].D;
				if(Pe > 1.0)
				Pe = 1.0;
			Sigma = d * Pe / a_L2 / 2.0;
		        for(J = 0; J < I; J++) {
				sum = DATA2[3*J+0]*DATA3[3*J+0] + DATA2[3*J+1]*DATA3[3*J+1] + DATA2[3*J+2]*DATA3[3*J+2];
				sum *= (DATA2[3*J+0]*DATA4[3*J+0] + DATA2[3*J+1]*DATA4[3*J+1] + DATA2[3*J+2]*DATA4[3*J+2]);
				res += (-Pow(ion[which_ion].D * ion[which_ion].Z,2) * Sigma * sum)* quad->weights[J];
			}
			*value = res * h;
		}
	}
	else *value = 0.0;
}

static void Region_Quad_SUPG_term_Slotboom(int region, SIMPLEX *e, DOF *grad_u, DOF *u, DOF *size_p, DOF *P, int i, int j, int order, FLOAT *value) {
	INT I, J, K;
	QUAD *quad = phgQuadGetQuad3D(order);
	GRID *g = P->g;
	const FLOAT *DATA0, *DATA1, *DATA2, *DATA3, *DATA4;
	FLOAT h = phgGeomGetVolume(g, e);
	FLOAT d = phgGeomGetDiameter(g, e);
	DATA0 = phgQuadGetDofValues(e, size_p, quad);
	DATA1 = phgQuadGetDofValues(e, u, quad);
        DATA2 = phgQuadGetDofValues(e, grad_u, quad);
        DATA3 = phgQuadGetBasisGradient(e, P, j, quad);
        DATA4 = phgQuadGetBasisGradient(e, P, i, quad);
        I = quad->npoints;
	FLOAT Pe, Sigma, res = 0.0;
	FLOAT sum, a_L2 = 0.0;
	if(e->region_mark == region) {
		for(J = 0; J < I; J++) {
			for(K = 0; K < 3; K++) {
				a_L2 += (DATA2[3*J+K]*DATA2[3*J+K]);
			}
			a_L2 *= quad->weights[J];
		}
		if(a_L2 != 0) {
			a_L2 = Sqrt(a_L2 * h) * ion[which_ion].D * ion[which_ion].Z;
			Pe = a_L2 * d / 6.0 / ion[which_ion].D;
			if(Pe > 1.0)
				Pe = 1.0;
			Sigma = d * Pe / a_L2 / 2.0;
		        for(J = 0; J < I; J++) {
				sum = DATA2[3*J+0]*DATA3[3*J+0] + DATA2[3*J+1]*DATA3[3*J+1] + DATA2[3*J+2]*DATA3[3*J+2];
				sum *= (DATA2[3*J+0]*DATA4[3*J+0] + DATA2[3*J+1]*DATA4[3*J+1] + DATA2[3*J+2]*DATA4[3*J+2]);
				sum *= Pow(DATA0[J], K_i[which_ion]) * Exp(-1.0 * ion[which_ion].Z * DATA1[J]);
				res += (-Pow(ion[which_ion].D * ion[which_ion].Z, 2) * Sigma * sum)* quad->weights[J];
			}
			*value = res * h;
		}
	}
	else *value = 0.0;
}

static void get_exp_u(DOF **tmp, DOF *u, PNP_FUNC_1 func) {
	INT i, N;
	FLOAT value;
	N = DofGetDataCount(*tmp);
	for(i = 0; i < N; i++) {
		value = u->data[i];
		func(&value);
		(*tmp)->data[i] = value;
	}
}

//static void get_size_p(DOF **tmp, DOF **p, PNP_FUNC_2 func) {
void get_size_p(DOF **tmp, DOF **p, PNP_FUNC_2 func) {
	INT i, N;
	int j;
	N = DofGetDataCount(*tmp);
	FLOAT *data;
	data = (FLOAT *)malloc(MNION * sizeof(FLOAT));
	memset(data, 0, MNION * sizeof(*data));

	for(i = 0; i < N; i++) {
		for(j = 0; j < NION; j++) {
			data[j] = p[j]->data[i];
		}
		func(data, ((*tmp)->data+i));
//		printf("%lf",(*tmp)->data[i]);
                //if(MOD == 0 && nonlinear_Epsilon == 0){
                if(MOD == 0 || (MOD ==2 && nonlinear_Epsilon == 0)){
                        (*tmp)->data[i] = 1;
		}
	}
	free(data);
}

static void p_Transform(FLOAT Alpha, DOF *p, DOF *P, DOF *u, DOF *size_p) {
	DOF *tmp = phgDofCopy(p, NULL, NULL, NULL);
	DOF *tmp_exp_u = phgDofCopy(u, NULL, NULL, NULL);
	DOF *size_p_k = phgDofCopy(p, NULL, NULL, NULL);;
	get_exp_u(&tmp_exp_u, u, func_exp_1);
	dof_pow(&size_p_k, size_p, K_i[which_ion]);
	dof_multiple(&tmp, tmp_exp_u, size_p_k, P, NULL);
	phgDofAXPBY(1.0 - Alpha, tmp, Alpha, &p);

	phgDofFree(&tmp);
	phgDofFree(&tmp_exp_u);
	phgDofFree(&size_p_k);
}

static void P_Transform(FLOAT Alpha, DOF *p, DOF *P, DOF *u, DOF *size_p) {
	DOF *tmp = phgDofCopy(P, NULL, NULL, NULL);
	DOF *tmp_exp_u = phgDofCopy(u, NULL, NULL, NULL);
	DOF *size_p_k = phgDofCopy(p, NULL, NULL, NULL);;
	get_exp_u(&tmp_exp_u, u, func_exp_1);
	dof_pow(&tmp_exp_u, tmp_exp_u, -1);
	dof_pow(&size_p_k, size_p, -K_i[which_ion]);
	dof_multiple(&tmp, tmp_exp_u, size_p_k, p, NULL);
	phgDofAXPBY(1.0 - Alpha, tmp, Alpha, &P);

	phgDofFree(&tmp);
	phgDofFree(&tmp_exp_u);
	phgDofFree(&size_p_k);
}

static void t_Transform(FLOAT Alpha, DOF *p, DOF *T) {
	DOF *tmp = phgDofCopy(p, NULL, NULL, NULL);
	get_exp_u(&tmp, T, func_exp_4);
	dof_const(&tmp, tmp, bulk[which_ion]);
	phgDofAXPBY(1.0 - Alpha, tmp, Alpha, &p);

	phgDofFree(&tmp);
}

static void mark_region(GRID *g, int region, int level) {
	SIMPLEX *e;
	int i;
	ForAllElements(g, e) {
		if(e->region_mark == region) {
			for(i = 0; i < NFace; i++) {
				if(e->bound_type[i] & bmlc_bdry) {
					e->mark = level;
					break;
				}
			}
		}
	}
}

void get_Epsilon_all(DOF **Full_Epsilon, DOF *Epsilon, DOF **p) {
        INT i, N;
        int j;
	FLOAT tmp;
        phgDofAXPBY(1.0, Epsilon, 0.0, Full_Epsilon);
        N = DofGetDataCount(*Full_Epsilon);

        for(i = 0; i < N; i++) {
		tmp = 0.0;  // for nonlinear form
                for(j = 0; j < NION; j++) {
			tmp += Fabs(p[j]->data[i]);
		}
		if(!tmp) {
			(*Full_Epsilon)->data[i] = Epsilon_m;
		}
	}
}

static void get_Epsilon_all_test(DOF **Full_Epsilon, DOF *Epsilon) {
        GRID *g = (Epsilon)->g;
        SIMPLEX *e;
        int i = 0;
	int bdry_flag = 0;
	INT count;
	phgDofAXPBY(1.0, Epsilon, 0.0, Full_Epsilon);
	ForAllElements(g, e) {
//		if(e->region_mark == solute_region) {
                if(e->region_mark == solute_region || e->region_mark == membrane_region) {
			bdry_flag = 0;
			for(i = 0; i < NFace; i++) {
				if(e->bound_type[i] == bmlc_bdry/* write solute boundary here: Omega_m */)
					bdry_flag = 1;
			}
			for(i = 0; i < NVert; i++) {
				count = e->verts[i];
				if(bdry_flag == 0) {
					(*Full_Epsilon)->data[count] = Epsilon_m;
				}
				else if(e->bound_type[i] == bmlc_bdry/* same solute boundary */) {
					(*Full_Epsilon)->data[count] = Epsilon_m;
				}
			}
		}
	}
}


//00000000000000002021-11-14 rgshen
void get_Epsilon_rgshen(DOF **tmp, DOF **p) {
        INT i, N;
        int j;
        N = DofGetDataCount(*tmp);

        for(i = 0; i < N; i++) {
                
                if (MOD == 2 && nonlinear_Epsilon == 1) {
                    (*tmp)->data[i] = 1.0;  // for nonlinear form
                    for(j = 0; j < NION; j++) {
                            (*tmp)->data[i] += asize_ratio * Pow(ion[j].a, 3) * p[j]->data[i];  // for nonlinear
                    }
                    (*tmp)->data[i] = (Epsilon_s-Epsilon_m) / (*tmp)->data[i] + Epsilon_m;  // for nonlinear
                }
                
        }
		
		phgPrintf("asize_ratio = %.8f \n", asize_ratio);
}


//0000000000000000000000000




// for MOD = 2 or MOD = 3 i.e. variable dielectric model: epsilon is the function of concentration p in the solvent region
//static void get_Epsilon(DOF **tmp, DOF **p) {
void get_Epsilon(DOF **tmp, DOF **p) {
        INT i, N;
        int j;
        N = DofGetDataCount(*tmp);

        for(i = 0; i < N; i++) {
                if (MOD == 2 && nonlinear_Epsilon == 2) {
                        (*tmp)->data[i] = 1.0;  // for nonlinear form
                        for(j = 0; j < NION; j++) {
                                (*tmp)->data[i] += A * p[j]->data[i] / (Pow(10, 3) * NA) / NION;  // for nonlinear
                        }
                        (*tmp)->data[i] = Epsilon_s / (*tmp)->data[i];  // for nonlinear
                }
                if (MOD == 2 && nonlinear_Epsilon == 1) {
                        (*tmp)->data[i] = 1.0;  // for nonlinear form
                        for(j = 0; j < NION; j++) {
                                (*tmp)->data[i] += Pow(ion[j].a, 3) * p[j]->data[i] / Pow(change_unit * mesh_unit, 3);  // for nonlinear
                        }
                        (*tmp)->data[i] = (Epsilon_s-Epsilon_m) / (*tmp)->data[i] + Epsilon_m;  // for nonlinear
                }
                //else{
                if (MOD == 2 && nonlinear_Epsilon == 0) {
                        (*tmp)->data[i] = Epsilon_s;  // for linear form
                        for(j = 0; j < NION; j++) {
                                (*tmp)->data[i] -= (Epsilon_s - Epsilon_m) * Pow(ion[j].a, 3) * p[j]->data[i] / Pow(change_unit * mesh_unit, 3);  // for linear
                        }
                }
        }
}
/* for MOD =4 i.e. epsilon is the function of space z in the solvent region */
//static void get_Epsilon_r(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
void get_Epsilon_r(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
	FLOAT a_1 = Epsilon_c / Epsilon_s;
	FLOAT z_0 = (channel_top + channel_top - buffer_length_eps) / 2.0;
	if(Epsilon_type == 1) {
		if(z <= -channel_bottom || z >= channel_top ) {
			*value = Epsilon_s;
        	}
        	else if(z >= -(channel_bottom - buffer_length_eps) && z <= (channel_top - buffer_length_eps)) {
			*value = Epsilon_c;
        	}
		else if(z >= - channel_bottom && z <= -(channel_bottom - buffer_length_eps)) {
			*value = (Epsilon_s - Epsilon_c)/(2*atan(s_0)) * atan(2 * s_0 / (-buffer_length_eps) * z - s_0 * (-(channel_bottom - buffer_length_eps) - channel_bottom) / (-buffer_length_eps)) + (Epsilon_s + Epsilon_c) / 2.0;
        	}
		else if(z >= (channel_top - buffer_length_eps) && z <= channel_top) {
                	*value = (Epsilon_s - Epsilon_c)/(2*atan(s_0)) * atan(2 * s_0 / buffer_length_eps * z - s_0 * (channel_top - buffer_length_eps + channel_top) / buffer_length_eps) + (Epsilon_s + Epsilon_c) / 2.0;
        	}
	}
	else if(Epsilon_type == 2) {
		*value = Epsilon_s * (a_1 + (1.0 - a_1) / (1.0 + exp(-(abs(z) - z_0) / delta_z_E)));
	}
}

void get_G_u_DOT_G_v(DOF **tmp, DOF *u, DOF *v) {
        INT i, N;
        int j;
        N = DofGetDataCount(*tmp);
        DOF *grad_u = phgDofGradient(u, NULL, NULL, NULL);
        DOF *grad_v = phgDofGradient(v, NULL, NULL, NULL);
        //N = DofGetDataCount(grad_u);
	//N = N/3;

        for(i = 0; i < N; i++) {
                (*tmp)->data[i] = 0;
                for(j = 0; j < 3; j++) {
                        (*tmp)->data[i] += grad_u->data[3*i+j] * grad_v->data[3*i+j];
                }
		(*tmp)->data[i] *= Pow(change_unit * mesh_unit, 2);
        }
        phgDofFree(&grad_u);
        phgDofFree(&grad_v);
}

//static void get_G_D_DOT_G_D(DOF **tmp, DOF *dof) {
void get_G_D_DOT_G_D(DOF **tmp, DOF *dof) {
        INT i, N;
        int j;
        N = DofGetDataCount(*tmp);
        DOF *grad_dof = phgDofGradient(dof, NULL, NULL, NULL);

        for(i = 0; i < N; i++) {
                (*tmp)->data[i] = 0;
                for(j = 0; j < 3; j++) {
                        (*tmp)->data[i] += Pow(grad_dof->data[3*i+j], 2) * Pow(change_unit * mesh_unit, 2);
                }
        }
        phgDofFree(&grad_dof);
}

static void See_bdry(DOF **bdry) {
        GRID *g = (*bdry)->g;
        SIMPLEX *e;
        INT i, N;

        ForAllElements(g, e) {
                N = DofGetNBas(*bdry, e);
                for(i = 0; i < N; i++ ) {
//			printf("%d\n",phgDofGetElementBoundaryType(*bdry, e, i));
			if(phgDofGetElementBoundaryType(*bdry, e, i) & BDRY_USER8) {
				int all = phgDofMapE2D(*bdry, e, i);
				(*bdry)->data[all] = 8;
			}
			if(phgDofGetElementBoundaryType(*bdry, e, i) & BDRY_USER7) {
				int all = phgDofMapE2D(*bdry, e, i);
				(*bdry)->data[all] = 7;
			}
			if(phgDofGetElementBoundaryType(*bdry, e, i) & BDRY_USER6) {
				int all = phgDofMapE2D(*bdry, e, i);
				(*bdry)->data[all] = 6;
			}
			if(phgDofGetElementBoundaryType(*bdry, e, i) & BDRY_USER5) {
				int all = phgDofMapE2D(*bdry, e, i);
				(*bdry)->data[all] = 5;
			}
			if(phgDofGetElementBoundaryType(*bdry, e, i) & BDRY_USER4) {
				int all = phgDofMapE2D(*bdry, e, i);
				(*bdry)->data[all] = 4;
			}
			if(phgDofGetElementBoundaryType(*bdry, e, i) & BDRY_USER3) {
				int all = phgDofMapE2D(*bdry, e, i);
				(*bdry)->data[all] = 3;
			}
			if(phgDofGetElementBoundaryType(*bdry, e, i) & BDRY_USER2) {
				int all = phgDofMapE2D(*bdry, e, i);
				(*bdry)->data[all] = 2;
			}
			if(phgDofGetElementBoundaryType(*bdry, e, i) & BDRY_USER1){
				int all = phgDofMapE2D(*bdry,e,i);
				(*bdry)->data[all] = 1;
			}
			if(phgDofGetElementBoundaryType(*bdry, e, i) & BDRY_USER0){
				int all = phgDofMapE2D(*bdry,e,i);
				(*bdry)->data[all] = -1;
			}
//			*((*bdry)->data+all)=phgDofGetElementBoundaryType(*bdry, e, i);
		}
	}

}

//approximation by Sharp and Honig, 1990: 0.5 * Sigma(Rho_f * Phi)
static void Calculate_energy(DOF *u, DOF *H) {
	int i, j;
	COORD position[NATOMS];
	FLOAT u_value[NATOMS], H_value[NATOMS];
	FLOAT res = 0.0, sum = 0.0;
	for(i = 0; i < NATOMS; i++){
		position[i][0] = atom[i].x;
		position[i][1] = atom[i].y;
		position[i][2] = atom[i].z;
	}
	phgInterGridDofEval(u, NATOMS, position, u_value, 0);
	phgInterGridDofEval(H, NATOMS, position, H_value, 0);
	for(i = 0; i < NATOMS; i++) {
		sum += (u_value[i] + H_value[i]) * atom[i].Z;
		phgPrintf("\nu_value[%d] = %le  H_value[%d] = %le  atom[%d] = %le", i, u_value[i], i, H_value[i], i, atom[i].Z);
		//energy between fixed ions
		for(j = i + 1; j < NATOMS; j++) {
			res += atom[i].Z * atom[j].Z / r2(atom[i].x, atom[i].y, atom[i].z, atom[j].x, atom[j].y, atom[j].z);
		}
	}
	if(!SOU) {
		res *= Pow(ec, 2) * change_unit / Epsilon_m;
	}
	else {
		res *= Pow(ec, 2) * change_unit / (4 * M_PI * Epsilon_m * Epsilon_vacuum);
	}
	sum *= 0.5 * NA / Beta / kcal;
	res *= 0.5 * NA / kcal;
	phgPrintf("\n\nSolvation energy = %le	kcal/mol\n", sum);
	phgPrintf("Fixed atoms energy = %le	kcal/mol\n", res);
}

static void Calculate_chemical_potential(DOF **tmp, DOF *u, DOF **p, DOF *Epsilon, int j) {
	INT i, N;
	N = DofGetDataCount(u);
	for(i=0; i < N; i++) {
		if(p[j]->data[i] < 0.0000) {
			p[j]->data[i] = 1e-6;
                }
		//if(ion[j].c==0)phgError(1, "ionc");
		//if(ion[j].a==0)phgError(1, "iona");
		if(p[j]->data[i]) {
			if(MOD == 4) {
				(*tmp)->data[i] = ion[j].Z * u->data[i] + Log(p[j]->data[i] / ion[j].c) + Mu * Kappa_2 * change_unit * Pow(ion[j].Z,2) / (2 * ion[j].a) * (1/(Epsilon->data[i]) - 1);
			}
			else if(MOD == 0) {
				(*tmp)->data[i] = ion[j].Z * u->data[i] + Log(p[j]->data[i] / ion[j].c);
			}
		}
		else (*tmp)->data[i] = 0;
	}
}

static void Calculate_potential_energy(DOF **tmp, DOF *u, DOF *Epsilon, int j) {
        INT i, N;
        N = DofGetDataCount(u);
        for(i = 0; i < N; i++) {
                if(MOD == 4) {
                	(*tmp)->data[i] = ion[j].Z * u->data[i] + Mu * Kappa_2 * change_unit * Pow(ion[j].Z,2) / (2 * ion[j].a) * (1/(Epsilon->data[i]) - 1);
		}
		else if(MOD == 0) {
			(*tmp)->data[i] = ion[j].Z * u->data[i];
		}
	}
}


FLOAT PNP_quad_channel_section(ELEMENT *e, FLOAT npoint, COORD *points, INT grad_dofs, DOF *dof, ...) {
	//points are coordinate of A, B, C
	//(the three points of the triangle ABC)
	//points[3][ ] for A B C or points[4][ ] for A B C D
	//points[ ][3] for x y z
	//grad_dofs is used to calculate gradient of dof and select z axis


	GRID *g = dof->g;
	int d = 2;

	//get dofs
	DOF **dofs;
	int i, j, k, ndof;
	dofs = phgAlloc(8 * sizeof(*dofs));

	va_list ap;
	va_start(ap, dof);
	for(ndof = 0; ndof < 256; ndof++) {
		if(dof == NULL) break;
		dofs[ndof] = dof;
		dof = va_arg(ap, DOF *);
	}
	dofs[ndof] = NULL;
	va_end(ap);

	//change coordinate
	FLOAT point_lambda[3][4];
	FLOAT *tmp_point;
	for(i = 0; i < 3; i++) {
		tmp_point = phgGeomXYZ2Lambda(g, e, points[i][0], points[i][1], points[i][2]);
		for(j = 0; j < Dim + 1; j++)
			point_lambda[i][j] = tmp_point[j];
	}

	FLOAT area;
	FLOAT *X = (FLOAT *)malloc(3 * sizeof(FLOAT));
	FLOAT *Y = (FLOAT *)malloc(3 * sizeof(FLOAT));
	for(i = 0; i < 3; i++) {
		X[i] = points[i][(d+1)%3];
		Y[i] = points[i][(d+2)%3];
	}
	triangle_area(X, Y, &area);

	//face quad
	int order;
	QUAD *quad;
	FLOAT res, sum, lambda[Dim + 1];
	FLOAT dof_value[1], dof_grad[Dim];
	const FLOAT *p, *w;

	order = 0;
	for(i = 0; i < ndof; i++) {
//phgPrintf("\n DofOrder = %d\n", DofTypeOrder(dofs[i], e));
		if(DofTypeOrder(dofs[i], e) > 0)
			order += DofTypeOrder(dofs[i], e);
	}
	order -= grad_dofs;
//phgPrintf("\n QuadOrder = %d\n", order);
	//if(order < 1) order = 1;
	if(test_current == 1){
		//if(MOD = 2){
		order = order + 1;
		phgPrintf("\n QuadOrder = %d\n", order);
		//}
		//else if(MOD == 0){
		//	order = 4;
		//}
	}
//phgPrintf("\n QuadOrder = %d\n", order);
	quad = phgQuadGetQuad2D(order);
	p = quad->points;
	w = quad->weights;

	sum = 0.0;
	for(i = 0; i < quad->npoints; i++) {
		for(j = 0; j < Dim + 1; j++) {
		    lambda[j] = 0.0;
		    for(k = 0; k < 3; k++) {
			lambda[j] += point_lambda[k][j] * p[k + i * 3];
			if(lambda[j] > 1.1) phgError(1, "lambda > 1");
		    }
		}
		res = 1.0;
		for(j = 0; j < ndof; j++) {
		    if(grad_dofs && !j) {
			phgDofEvalGradient(dofs[j], e, lambda, NULL, dof_grad);
			res *= dof_grad[d];
		    }
		    else {
			phgDofEval(dofs[j], e, lambda, dof_value);
			res *= dof_value[0];
		    }
		}
		sum += res * *(w++);
	}

	phgFree(dofs);

	return sum * area;
}

// rgshen已将下述函数调整为标准（规范）缩进对齐---2020-11-03
static void VD_Calculate_current_u_D_born(int region, int d, DOF *u, DOF *u_D, DOF **p, FLOAT value, DOF *size_p, DOF *Epsilon) {
        FLOAT sum_current = 0;
        FLOAT p_current[NION];
        FLOAT current[NION];
        FLOAT sum_current1 = 0;
        FLOAT sum_current2 = 0;
        FLOAT sum_current3 = 0;
        FLOAT sum_current4 = 0;
        FLOAT p_current1[NION];
        FLOAT p_current2[NION];
        FLOAT p_current3[NION];
        FLOAT p_current4[NION];
        FLOAT current1[NION];
        FLOAT current2[NION];
        FLOAT current3[NION];
        FLOAT current4[NION];  //for born correct

        SIMPLEX *e;
        GRID *g = u->g;
        FLOAT area, volume, *X, *Y;
        FLOAT *lambda, e_lambda[4] = {0.25, 0.25, 0.25, 0.25};
        FLOAT sum_area = 0;
        FLOAT total_area = 0;

        COORD triangle_verts[4];
        COORD tetrahetron_verts[5];

        DOF *tmp, *tmp_eps, *E, *G_u_D_G_v, *div_grad_eps;
	    DOF *tmp_tmp, *multiple;
        DOF *TMP_u_D = phgDofCopy(u_D, NULL, NULL, NULL);
        DOF *tmp_u_D = phgDofCopy(u_D, NULL, NULL, NULL);
        tmp = phgDofNew(g, DOF_DEFAULT, 1, "tmp", DofInterpolation);
        tmp_eps = phgDofNew(g, DOF_DEFAULT, 1, "tmp_eps", DofInterpolation);
        if(!strcmp(u->type->name, "P1")) {
			G_u_D_G_v = phgDofNew(g, DOF_DG0, 1, "G_u_D_G_v", DofInterpolation);
	        div_grad_eps = phgDofNew(g, DOF_DG0, 1, "div_grad_eps", DofInterpolation);
               	tmp_tmp = phgDofNew(u->g, DOF_DG1, 1, "tmp_tmp", DofInterpolation);
               	multiple = phgDofNew(u->g, DOF_DG1, 1, "tmp_grad_u_2", DofInterpolation);
		}
        if(!strcmp(u->type->name, "P2")) {
			G_u_D_G_v = phgDofNew(g, DOF_DG1, 1, "G_u_D_G_v", DofInterpolation);
        	div_grad_eps = phgDofNew(g, DOF_DG0, 1, "div_grad_eps", DofInterpolation);
			if(born_correct) {
				tmp_tmp = phgDofNew(u->g, DOF_DG1, 1, "tmp_tmp", DofInterpolation);
                multiple = phgDofNew(u->g, DOF_DG1, 1, "tmp_grad_u_2", DofInterpolation);
			}
			else {
                tmp_tmp = phgDofNew(u->g, DOF_DG2, 1, "tmp_tmp", DofInterpolation);
                multiple = phgDofNew(u->g, DOF_DG2, 1, "tmp_grad_u_2", DofInterpolation);
			}
		}
        if(!strcmp(u->type->name, "P3")) {
			G_u_D_G_v = phgDofNew(g, DOF_DG2, 1, "G_u_D_G_v", DofInterpolation);
        	div_grad_eps = phgDofNew(g, DOF_DG1, 1, "div_grad_eps", DofInterpolation);
		    if(born_correct) {
                tmp_tmp = phgDofNew(u->g, DOF_DG2, 1, "tmp_tmp", DofInterpolation);
                multiple = phgDofNew(u->g, DOF_DG2, 1, "tmp_grad_u_2", DofInterpolation);
			}
			else {
                tmp_tmp = phgDofNew(u->g, DOF_DG3, 1, "tmp_tmp", DofInterpolation);
                multiple = phgDofNew(u->g, DOF_DG3, 1, "tmp_grad_u_2", DofInterpolation);
			}
		}

        phgDofSetDataByValue(tmp, 0.0);
        phgDofSetDataByValue(tmp_eps, 0.0);
        phgDofSetDataByValue(G_u_D_G_v, 0.0);
        phgDofSetDataByValue(div_grad_eps, 0.0);
        phgDofSetDataByValue(tmp_tmp, 0.0);
        phgDofSetDataByValue(multiple, 0.0);
        phgDofAXPBY(1.0, u, -1, &tmp_u_D);    //tmp_u_D = u - u_D
        phgDofAXPBY(1.0, u, -0.5, &TMP_u_D);  //TMP_u_D = u - 0.5*u_D
        if(MOD == 2) {
			get_G_u_DOT_G_v(&G_u_D_G_v, u, tmp_u_D);
            //TMP = phgDofGradient(G_u_D_G_u, NULL, NULL, NULL);//TMP = grad(grad_u dot grad_u)
            if (nonlinear_Epsilon == 1 || nonlinear_Epsilon == 2) {
                dof_pow(&tmp, size_p, -2);
            }
        	phgDofCopy(tmp, &tmp_tmp, NULL, NULL);
        	phgDofCopy(G_u_D_G_v, &multiple, NULL, NULL);
        	dof_multiple(&multiple, tmp_tmp, multiple, NULL);
        }
	    if(MOD == 4) {
	    	if(born_correct) {
                dof_pow(&tmp, Epsilon, -2);
                DOF *grad_eps = phgDofGradient(Epsilon, NULL, NULL, NULL);
                div_grad_eps = phgDofDivergence(grad_eps, NULL, NULL, NULL);

                //get_G_D_DOT_G_D(&G_u_D_G_v, Epsilon);
	    		phgDofCopy(tmp, &tmp_tmp, NULL, NULL);
	    		phgDofCopy(div_grad_eps, &multiple, NULL, NULL);
                dof_multiple(&multiple, tmp_tmp, multiple, NULL);

                E = phgDofNew(g, DOF_CONSTANT, 1, "unit", DofNoAction);
                phgDofSetDataByValue(E, 1.0);
                dof_pow(&tmp_eps, Epsilon, -1);
                phgDofAXPY(-1/Epsilon_m, E, &tmp_eps); //tmp_eps = 1/Epsilon - 1/Epsilon_m
	    		phgDofFree(&grad_eps);
            }
            else {
	            E = phgDofNew(g, DOF_CONSTANT, 1, "unit", DofNoAction);
                phgDofSetDataByValue(E, 1.0);
                dof_pow(&tmp_eps, Epsilon, -1);
                phgDofAXPY(-1/Epsilon_m, E, &tmp_eps); //tmp_eps = 1/Epsilon - 1/Epsilon_m
	    	}
	    }
        int i, j, k, l, vert;
        int type_flag;//number of p_tab
        int p_tab[NVert];//save the tab of verts of a ELEMENT which has coordinate on d direction larger than value
        int n_tab[NVert];//if smaller than value
        // phgPrintf("\n%d", type_flag);
        FLOAT grad_u[Dim];
        FLOAT grad_p[NION][Dim];
        FLOAT grad_eps[Dim];
        FLOAT p_value[NION];
        FLOAT grad_G2_value[Dim];
        FLOAT tmp_value[1];
        FLOAT multiple_value[Dim];

        memset(current, 0, (NION)*sizeof(FLOAT));
        memset(current1, 0, (NION)*sizeof(FLOAT));
        memset(current2, 0, (NION)*sizeof(FLOAT));
        memset(current3, 0, (NION)*sizeof(FLOAT));
        memset(current4, 0, (NION)*sizeof(FLOAT));
        memset(p_current, 0, (NION)*sizeof(FLOAT));
        memset(p_current1, 0, (NION)*sizeof(FLOAT));
        memset(p_current2, 0, (NION)*sizeof(FLOAT));
        memset(p_current3, 0, (NION)*sizeof(FLOAT));
        memset(p_current4, 0, (NION)*sizeof(FLOAT));
        //phgPrintf("\n%d", type_flag);
        ForAllElements(g, e) {
			if(e->region_mark == region) {
                memset(p_tab, 0, NVert*sizeof(int));
                memset(n_tab, 0, NVert*sizeof(int));
                type_flag = 0;
                i = 0;
                //phgPrintf("@@@@@@%\n%d", NVert);
                for(vert = 0; vert < NVert; vert++) {
                    if(g->verts[e->verts[vert]][d] >= value) {
						p_tab[type_flag++] = vert;
                    }
                    else {
                        n_tab[i++] = vert;
                    }
                    //phgPrintf("\n%d", type_flag);
                }
                //phgPrintf("\n%d", type_flag);

                if(type_flag == 1 || type_flag == 3) {
                    //COORD *triangle_verts = (COORD *)malloc(4 * sizeof(COORD));
                    memset(triangle_verts[3], 0, Dim * sizeof(FLOAT));
                    for(i = 0; i < 3; i++) {
                        if(type_flag == 1)
                            linear_interpolation(g->verts[e->verts[p_tab[0]]],
                                                 g->verts[e->verts[n_tab[i]]],
                                                 triangle_verts[i], value, d);
                        if(type_flag == 3)
                            linear_interpolation(g->verts[e->verts[n_tab[0]]],
                                                 g->verts[e->verts[p_tab[i]]],
                                                 triangle_verts[i], value, d);
                        for(j = 0; j < Dim; j ++) {
                            triangle_verts[3][j] += triangle_verts[i][j] / 3.0;
                        phgPrintf("****, %e", triangle_verts[3][j]);
                        }
                    }
                    lambda = phgGeomXYZ2Lambda(g,e,triangle_verts[3][0],triangle_verts[3][1],triangle_verts[3][2]);
                    X = (FLOAT *)malloc(3 * sizeof(FLOAT));
                    Y = (FLOAT *)malloc(3 * sizeof(FLOAT));
                    for(i = 0; i < 3; i++) {
                        X[i] = triangle_verts[i][(d+1)%3];
                        Y[i] = triangle_verts[i][(d+2)%3];
                    }
                    triangle_area(X, Y, &area);
                    //free(triangle_verts);
                    //phgDofEval(flux, e, e_lambda, &value);
                    //current += area * value;
                }
                else if(type_flag == 2) {
                    //COORD *tetrahetron_verts = (COORD *)malloc(5 * sizeof(COORD));
                    memset(tetrahetron_verts[4], 0, Dim * sizeof(FLOAT));
                    for(i = 0; i < 2; i++) {
                        for(j = 0; j < 2; j++) {
                            linear_interpolation(g->verts[e->verts[p_tab[i]]],
                                                 g->verts[e->verts[n_tab[j]]],
                                                 tetrahetron_verts[2*i+j], value, d);
                            for(k = 0; k < Dim; k ++) {
                                tetrahetron_verts[4][k] += tetrahetron_verts[2*i+j][k] / 4.0;
                            }
                        }
                    }
                    lambda = phgGeomXYZ2Lambda(g,e,tetrahetron_verts[4][0],tetrahetron_verts[4][1],tetrahetron_verts[4][2]);
                    X = (FLOAT *)malloc(4 * sizeof(FLOAT));
                    Y = (FLOAT *)malloc(4 * sizeof(FLOAT));
                    for(i = 0; i < 4; i++) {
                        X[i] = tetrahetron_verts[i][(d+1)%3];
                        Y[i] = tetrahetron_verts[i][(d+2)%3];
                    }
                    tetrahetron_area(X, Y, &area);
                    //free(tetrahetron_verts);
                    //phgPrintf("\n%e",area);
                    //phgDofEval(flux, e, e_lambda, &value);
                    //current += area * value;
                }
                if(type_flag > 0 && type_flag < 4) {
                    if(MOD == 2) {
                        if (nonlinear_Epsilon == 0) {
                            sum_area += area;
                            phgDofEvalGradient(TMP_u_D, e, lambda, NULL, grad_u);
                            phgDofEvalGradient(G_u_D_G_v, e, lambda, NULL, grad_G2_value);
                            for(i = 0; i < NION; i++) {
                                phgDofEval(p[i], e, lambda, p_value + i);
                                phgDofEvalGradient(p[i], e, lambda, NULL, grad_p[i]);
                            }
                            for(i = 0; i < NION; i++) {
							    if(current_switch) {
                                    current[i] += area * (grad_p[i][d] + ion[i].Z * p_value[i] * grad_u[d] + p_value[i] * (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2) * grad_G2_value[d]);
                                    current1[i] += area * grad_p[i][d];
                                    current2[i] += area * ion[i].Z * p_value[i] * grad_u[d];
                                    current3[i] += area * p_value[i] * (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2) * grad_G2_value[d];
								    //printf("\nold c1 = %le, c2 = %le, c3 = %le", current1[i], current2[i], current3[i]);
							    }
						        else {
								    if(type_flag != 2) {
									    current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL);
									    current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL);
									    current3[i] += (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2)
									    	           * PNP_quad_channel_section(e, 3, triangle_verts, 1, G_u_D_G_v, p[i], NULL);
								    }
								    else {
								    	for(j = 0; j < 4; j++) {
								    		for(k = 0; k < 3; k++) {
								    			for(l = 0; l < 3; l++) {
								    				triangle_verts[k][l] = tetrahetron_verts[(j+k)%4][l];
								    			}
								    		}
								    		current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL) / 2.0;
								    		current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL) / 2.0;
								    		current3[i] += (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2)
								    			            * PNP_quad_channel_section(e, 3, triangle_verts, 1, G_u_D_G_v, p[i], NULL) / 2.0;
								    	}
								    }
								    current[i] = current1[i] + current2[i] + current3[i];
								    //printf("  new c1 = %le, c2 = %le, c3 = %le\n", current1[i], current2[i], current3[i]);
							    }
                            }
                        }
                        else if (nonlinear_Epsilon == 1) {
                            sum_area += area;
                            //phgDofEval(tmp, e, lambda, tmp_value);
                            phgDofEvalGradient(TMP_u_D, e, lambda, NULL, grad_u);
                            phgDofEvalGradient(multiple, e, lambda, NULL, multiple_value);
                            for(i = 0; i < NION; i++) {
                                phgDofEval(p[i], e, lambda, p_value + i);
                                phgDofEvalGradient(p[i], e, lambda, NULL, grad_p[i]);
                            }
                            for(i = 0; i < NION; i++) {
							    if(current_switch) {
                                    //current[i] += area * (grad_p[i][d] + ion[i].Z * p_value[i] * grad_u[d] + p_value[i] * (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2) * tmp_value[0] * grad_G2_value[d]);
                                    current[i] += area * (grad_p[i][d] + ion[i].Z * p_value[i] * grad_u[d] + p_value[i] * (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2) * multiple_value[d]);
	                                current1[i] += area * grad_p[i][d];
	                                current2[i] += area * ion[i].Z * p_value[i] * grad_u[d];
	                                //current3[i] += area * p_value[i] * (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2) * tmp_value[0] * grad_G2_value[d];
	                                current3[i] += area * p_value[i] * (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2) * multiple_value[d];
							    }
							    else {
								    if(type_flag != 2) {
									    current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL);
									    current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL);
									    current3[i] += (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2)
									    	           * PNP_quad_channel_section(e, 3, triangle_verts, 1, tmp, G_u_D_G_v, p[i], NULL);
									    current3[i] += (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2)
									    	           * PNP_quad_channel_section(e, 3, triangle_verts, 1, G_u_D_G_v, tmp, p[i], NULL);
								    }
								    else {
								    	for(j = 0; j < 4; j++) {
								    		for(k = 0; k < 3; k++) {
								    			for(l = 0; l < 3; l++) {
								    				triangle_verts[k][l] = tetrahetron_verts[(j+k)%4][l];
								    			}
								    		}
								    		current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL) / 2.0;
								    		current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL) / 2.0;
								    		current3[i] += (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2)
								    			           * PNP_quad_channel_section(e, 3, triangle_verts, 1, tmp, G_u_D_G_v, p[i], NULL) / 2.0;
								    		current3[i] += (Epsilon_s - Epsilon_m) * Pow(ion[i].a / change_unit, 3) / (2.0 * Kappa_2)
								    			           * PNP_quad_channel_section(e, 3, triangle_verts, 1, G_u_D_G_v, tmp, p[i], NULL) / 2.0;
								    	}
								    }
								    current[i] = current1[i] + current2[i] + current3[i];
							    }
                            }
                        }
                        else if (nonlinear_Epsilon == 2) {
                            sum_area += area;
                            //phgDofEval(tmp, e, lambda, tmp_value);
                            phgDofEvalGradient(TMP_u_D, e, lambda, NULL, grad_u);
                            //phgDofEvalGradient(G_u_D_G_v, e, lambda, NULL, grad_G2_value);
						    //if(current_switch) {
                                phgDofEvalGradient(multiple, e, lambda, NULL, multiple_value);
						    //}
                            for(i = 0; i < NION; i++) {
                                phgDofEval(p[i], e, lambda, p_value + i);
                                phgDofEvalGradient(p[i], e, lambda, NULL, grad_p[i]);
                            }
                            for(i = 0; i < NION; i++) {
							    if(current_switch) {
                                    //current[i] += area * (grad_p[i][d] + ion[i].Z * p_value[i] * grad_u[d] + p_value[i] * Epsilon_s * A * 1.0 / (Pow(10, 3) * NA) / 2 / (2.0 * Kappa_2) * tmp_value[0] * grad_G2_value[d]);
                                    current[i] += area * (grad_p[i][d] + ion[i].Z * p_value[i] * grad_u[d] + p_value[i] * Epsilon_s * A * 1.0 / (Pow(10, 3) * NA) / NION / (2.0 * Kappa_2) * multiple_value[d]);
                                    current1[i] += area * grad_p[i][d];
                                    current2[i] += area * ion[i].Z * p_value[i] * grad_u[d];
                                    //current3[i] += area * p_value[i] * Epsilon_s * A * 1.0 / (Pow(10, 3) * NA) / 2 / (2.0 * Kappa_2) * tmp_value[0] * grad_G2_value[d];
                                    current3[i] += area * p_value[i] * Epsilon_s * A * 1.0 / (Pow(10, 3) * NA) / NION / (2.0 * Kappa_2) * multiple_value[d];
							    }
							    else {
							    	if(type_flag != 2) {
							    		current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL);
							    		current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL);
							    		current3[i] += Epsilon_s * A * 1.0 / (Pow(10, 3) * NA) / NION / (2.0 * Kappa_2)
							    			           * PNP_quad_channel_section(e, 3, triangle_verts, 1, tmp, G_u_D_G_v, p[i], NULL);
							    			     //* PNP_quad_channel_section(e, 3, triangle_verts, 1, tmp, multiple, p[i], NULL);
							    		current3[i] += Epsilon_s * A * 1.0 / (Pow(10, 3) * NA) / NION / (2.0 * Kappa_2)
							    			           * PNP_quad_channel_section(e, 3, triangle_verts, 1, G_u_D_G_v, tmp, p[i], NULL);
							    			     //* PNP_quad_channel_section(e, 3, triangle_verts, 1, multiple, tmp, p[i], NULL);
							    	}
							    	else {
							    		for(j = 0; j < 4; j++) {
							    			for(k = 0; k < 3; k++) {
							    				for(l = 0; l < 3; l++) {
							    					triangle_verts[k][l] = tetrahetron_verts[(j+k)%4][l];
							    				}
							    			}
							    			current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL) / 2.0;
							    			current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL) / 2.0;
							    			current3[i] += Epsilon_s * A * 1.0 / (Pow(10, 3) * NA) / NION / (2.0 * Kappa_2)
							    				           * PNP_quad_channel_section(e, 3, triangle_verts, 1, tmp, G_u_D_G_v, p[i], NULL) / 2.0;
							    			current3[i] += Epsilon_s * A * 1.0 / (Pow(10, 3) * NA) / NION / (2.0 * Kappa_2)
							    				           * PNP_quad_channel_section(e, 3, triangle_verts, 1, G_u_D_G_v, tmp, p[i], NULL) / 2.0;
							    		}
							    	}
							    	current[i] = current1[i] + current2[i] + current3[i];
							    }
                            }
                        }
                    }
				    else if(MOD == 4) {
					    sum_area += area;
                        phgDofEvalGradient(TMP_u_D, e, lambda, NULL, grad_u);
                        phgDofEvalGradient(tmp_eps, e, lambda, NULL, grad_eps);
                        phgDofEvalGradient(multiple, e, lambda, NULL, multiple_value);
                        for(i = 0; i < NION; i++) {
	                        phgDofEval(p[i], e, lambda, p_value + i);
                            phgDofEvalGradient(p[i], e, lambda, NULL, grad_p[i]);
                        }
                        for(i = 0; i < NION; i++) {
                            if(current_switch) {
                                current[i] += area * (grad_p[i][d] + ion[i].Z * p_value[i] * grad_u[d] +  Kappa_2 * Pow(ion[i].Z, 2) * Mu / (2 * ion[i].a / change_unit) * p_value[i] * grad_eps[d] + Kappa_2 * Pow(ion[i].Z, 2) / 2 * (ion[i].a / change_unit) / 6 * p_value[i] * multiple_value[d]);
                                current1[i] += area * grad_p[i][d];
                                current2[i] += area * ion[i].Z * p_value[i] * grad_u[d];
                                current3[i] += area * Kappa_2 * Pow(ion[i].Z, 2) * Mu / (2 * ion[i].a / change_unit) * p_value[i] * grad_eps[d];
                                current4[i] += area * Kappa_2 * Pow(ion[i].Z, 2) / 2 * (ion[i].a / change_unit) / 6 * p_value[i] * multiple_value[d];
                            }
                            else {
                                if(type_flag != 2) {
                                    current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL);
                                    current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL);
                                	current3[i] += Kappa_2 * Pow(ion[i].Z, 2) * Mu / (2 * ion[i].a / change_unit) * PNP_quad_channel_section(e, 3, triangle_verts, 1, tmp_eps, p[i], NULL);
                                	current4[i] += Kappa_2 * Pow(ion[i].Z, 2) / 2 * (ion[i].a / change_unit) / 6 * PNP_quad_channel_section(e, 3, triangle_verts, 1, div_grad_eps, tmp, p[i], NULL);
                                }
                                else {
                                    for(j = 0; j < 4; j++) {
                                        for(k = 0; k < 3; k++) {
                                            for(l = 0; l < 3; l++) {
                                                triangle_verts[k][l] = tetrahetron_verts[(j+k)%4][l];
                                            }
                                        }
                                        current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL) / 2.0;
                                        current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL) / 2.0;
                                		current3[i] += Kappa_2 * Pow(ion[i].Z, 2) * Mu / (2 * ion[i].a / change_unit) * PNP_quad_channel_section(e, 3, triangle_verts, 1, tmp_eps, p[i], NULL) / 2.0;
                                		current4[i] += Kappa_2 * Pow(ion[i].Z, 2) / 2 * (ion[i].a / change_unit) / 6 * PNP_quad_channel_section(e, 3, triangle_verts, 1, div_grad_eps, tmp, p[i], NULL) / 2.0;
                                    }
                                }
                                current[i] = current1[i] + current2[i] + current3[i] + current4[i];
                            }
                        }
				    }
                    else{
                        sum_area += area;
                        //phgDofEvalGradient(u, e, lambda, NULL, grad_u);
                        phgDofEvalGradient(TMP_u_D, e, lambda, NULL, grad_u);
                        for(i = 0; i < NION; i++) {
                            phgDofEval(p[i], e, lambda, p_value + i);
                            phgDofEvalGradient(p[i], e, lambda, NULL, grad_p[i]);
                        }
                        for(i = 0; i < NION; i++) {
                            if(current_switch) {
                            	current[i] += area * (grad_p[i][d] + ion[i].Z * p_value[i] * grad_u[d]);
                            	current1[i] += area * grad_p[i][d];
                            	current2[i] += area * ion[i].Z * p_value[i] * grad_u[d];
                            }
                            else {
                                if(type_flag != 2) {
                                	current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL);
                                    current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL);
                                }
                                else {
                                    for(j = 0; j < 4; j++) {
                	                    for(k = 0; k < 3; k++) {
        	                                for(l = 0; l < 3; l++) {
	                                            triangle_verts[k][l] = tetrahetron_verts[(j+k)%4][l];
                                            }
                                        }
                                        current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL) / 2.0;
                                        current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL) / 2.0;
                                    }
                                }
                                current[i] = current1[i] + current2[i];
                            }
                        }
                    }
                }
            }
        }

        phgPrintf("\n\nthe section value is %f A\n", value);
        for(i = 0; i < NION; i++) {
            //p: N/cm^3 or N/m^3 -> N/A^3 | D: A^2/ps -> A^2/s | I: A ->pA | grad and quad: mesh_unit -> A
            current[i] *= -ion[i].D * chan_ratio * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;
            current1[i] *= -ion[i].D * chan_ratio * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;
            current2[i] *= -ion[i].D * chan_ratio * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;
            current3[i] *= -ion[i].D * chan_ratio * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;
            current4[i] *= -ion[i].D * chan_ratio * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;
            //current[i] *= -ion[i].D * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;

            //phgPrintf("*********8%f\n",current[i]);
            MPI_Allreduce(&current[i], &p_current[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            MPI_Allreduce(&current1[i], &p_current1[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            MPI_Allreduce(&current2[i], &p_current2[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            MPI_Allreduce(&current3[i], &p_current3[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            MPI_Allreduce(&current4[i], &p_current4[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            MPI_Allreduce(&sum_area, &total_area, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            sum_current += p_current[i];
            sum_current1 += p_current1[i];
            sum_current2 += p_current2[i];
            sum_current3 += p_current3[i];
            sum_current4 += p_current4[i];
            phgPrintf("\nion[%d]    current is %e pA", i, p_current[i]);
            phgPrintf("\nion[%d]    the first part of current is %e pA", i, p_current1[i]);
            phgPrintf("\nion[%d]    the second part of current is %e pA", i, p_current2[i]);
            //phgPrintf("\nion[%d]    the third part of current is %e pA", i, p_current3[i]);
            //phgPrintf("\nion[%d]    the forth part of current is %e pA", i, p_current4[i]);
        }
        //print current section area
        phgPrintf("\n\ntotal    current is %e pA\n", sum_current);  // 这是什么意思呢？rgshen
        phgPrintf("\n\nthe first part    current is %e pA\n", sum_current1);
        phgPrintf("\n\nthe second part    current is %e pA\n", sum_current2);
        //phgPrintf("\n\nthe third part    current is %e pA\n", sum_current3);
        //phgPrintf("\n\nthe forth part    current is %e pA\n", sum_current4);
        phgPrintf("\n\ntotal area is %e A^2\n", total_area);
        if(phgRank == 1 && out_current_file != NULL) {
            FILE *f_in;
            f_in = fopen(out_current_file, "a");
            if(f_in == NULL) {
                phgError(1, "Can not open current output file %s.\n", out_current_file);
            }
            //fprintf(f_in, "%f   %le   %le   %le   %le   %le\n", value, sum_current, sum_current1, sum_current2, sum_current3, sum_current4);
//----------------------------------
#if 1
			fprintf(f_in, "%d DOF, %d elements, %d submeshes, load imbalance: %lg\n",
			        DofGetDataCountGlobal(u), g->nleaf_global, g->nprocs, (double)g->lif);//rgshen2020-11-12 这两处打印生成一次即注释掉
			fprintf(f_in, "u_bulk(V)	total_current(pA)	the_first_part(pA)	the_second_part(pA)\n"); //rgshen2020-11-12
#endif
			fprintf(f_in, "%f   %le       %le       %le\n", u_bulk, sum_current, sum_current1, sum_current2);//rgshen2020-11-12
//--------------------------------------------
            fclose(f_in);
        }
        //phgDofFree(&TMP);
        phgDofFree(&G_u_D_G_v);
        phgDofFree(&tmp_u_D);
        phgDofFree(&TMP_u_D);
        phgDofFree(&tmp);
	    phgDofFree(&tmp_eps);
        phgDofFree(&multiple);
       	phgDofFree(&tmp_tmp);
	    phgDofFree(&div_grad_eps);
	    if(MOD == 4) {
			phgDofFree(&E);
	    }
}


/*calculate the dof radial direction value*/
/*
FLOAT Calculate_dof_circle_value(int region, DOF *dof, FLOAT ra, FLOAT rb, FLOAT tr, FLOAT value) {
//sampling step length
        FLOAT r;
        INT np = 10;
        INT Ns = np;
        FLOAT tp = 2*M_PI/np;

        INT i, j, k, count;
        GRID *g = dof->g;
        COORD x[Ns];
        FLOAT tmp[Ns];
        FLOAT sam[Ns][3];
        INT N = (INT)((rb-ra)/tr);
	FLOAT dof_value[N];
//initialize sampling points
        for(i = 0; i < np; i++){
               sam[i][0] = Cos(i*tp);
               sam[i][1] = Sin(i*tp);
               sam[i][2] = value;
        }
        r = ra;
        for(i = 0; i < N; i++){
                for(j = 0; j < Ns; j++){
                                x[j][0] = r*sam[j][0];
                                x[j][1] = r*sam[j][1];
                                x[j][2] = sam[j][2];
                }
                bzero(tmp, Ns*sizeof(*tmp));
                phgInterGridDofEval(dof, Ns, x, tmp, 0);
                if(g->rank == 0){
			dof_value[i] = 0;
                        count = 0;
                        for(k=0;k<Ns;k++){
                                //if(region[k]==solute_region) continue;
                                count++;
                                dof_value[i] += tmp[k];
                        }
                        if(count==0) phgPrintf("\nERROR!Sampling points out of domain!\n");
                        else{
                                dof_value[i] /= count;
                        }
                }
                r += tr;
        }
//        phgPrintf("====================dof Radial Distribution====================\n");
        for(i=0;i<N;i++) phgPrintf("r %lf dof_value %lf\n", ra+i*tr, dof_value[i]);
}
*/

void lambda_reg(DOF *v, SIMPLEX *ee, int bno, const FLOAT *lambda, FLOAT *value){
        *value = ee->region_mark;
}

/*calculate the dof surface value*/
FLOAT Calculate_dof_circle_value(int region, DOF *dof, FLOAT r, FLOAT value) {
//sampling step length
        INT np = 20;
        INT Ns = np;
        FLOAT tp = 2*M_PI/np;

        INT i, j, k, count;
        GRID *g = dof->g;
        COORD x[Ns];
        FLOAT tmp[Ns], reg_m[Ns];
        FLOAT sam[Ns][3];
        FLOAT dof_value;
	//DOF *reg;
        //reg = phgDofNew(g, DOF_DEFAULT, 1, "region", DofLambdaFunction);
        //phgDofSetLambdaFunction(reg, lambda_reg);


//initialize sampling points
        for(i = 0; i < np; i++) {
               sam[i][0] = Cos(i*tp);
               sam[i][1] = Sin(i*tp);
               sam[i][2] = value;
        }
	    for(j = 0; j < Ns; j++) {
		    x[j][0] = r*sam[j][0];
            x[j][1] = r*sam[j][1];
            x[j][2] = sam[j][2];
        }
        bzero(tmp, Ns*sizeof(*tmp));
        //bzero(reg_m, Ns*sizeof(*reg_m));
        phgInterGridDofEval(dof, Ns, x, tmp, 0);
        //phgInterGridDofEval(reg, Ns, x, reg_m, 0);

        if(g->rank == 0) {
			dof_value = 0;
            count = 0;
            for(k=0; k<Ns; k++) {
			    //phgPrintf("\nreg_m[%d] = %f\n", k, reg_m[k]);
            	//if(tmp[k] < 0) phgPrintf("\nERROR!!!\n");
			     //if(reg_m[k] < 1.5) {
			     	//phgPrintf("\nsolute region!!!\n");
			     //	continue;
			     //}
                count++;
			    //if(tmp[k] < 0.0) tmp[k] = 1e-6;
                dof_value += tmp[k];
            }
            if(count == 0) phgPrintf("\nERROR!Sampling points out of domain!  section value = %.2f\n", value);
            else{
                dof_value /= count;
            }
        }
	//phgDofFree(&reg);
	return dof_value;
}

/*calculate the dof section value*/
FLOAT Calculate_dof_value(int region, int d, DOF *dof, FLOAT value) {
	SIMPLEX *e;
	GRID *g = dof->g;
	FLOAT area, volume, *X, *Y;
	//FLOAT *lambda, e_lambda[4] = {0.25, 0.25, 0.25, 0.25};
	FLOAT *lambda;
	//const FLOAT *lambda = NULL;
	FLOAT sum_area = 0;
	FLOAT total_area = 0;

	COORD triangle_verts[4];
	COORD tetrahetron_verts[5];

	int i, j, k, vert;
	int type_flag;//number of p_tab
	int p_tab[NVert];//save the tab of verts of a ELEMENT which has coordinate on d direction larger than value
	int n_tab[NVert];//if smaller than value

	FLOAT e_value, dof_value, sum_value;

	dof_value = 0;
	ForAllElements(g, e) {
		if(e->region_mark == region) {
			memset(p_tab, 0, NVert*sizeof(int));
			memset(n_tab, 0, NVert*sizeof(int));
			type_flag = 0;
			i = 0;
			for(vert = 0; vert < NVert; vert++) {
				if(g->verts[e->verts[vert]][d] >= value) {
					p_tab[type_flag++] = vert;
				}
				else {
					n_tab[i++] = vert;
				}
			}
			if(type_flag == 1 || type_flag == 3) {
				memset(triangle_verts[3], 0, Dim * sizeof(FLOAT));
				for(i = 0; i < 3; i++) {
					if(type_flag == 1)
						linear_interpolation(g->verts[e->verts[p_tab[0]]],
								     g->verts[e->verts[n_tab[i]]],
								     triangle_verts[i], value, d);
					if(type_flag == 3)
						linear_interpolation(g->verts[e->verts[n_tab[0]]],
								     g->verts[e->verts[p_tab[i]]],
								     triangle_verts[i], value, d);
					for(j = 0; j < Dim; j ++) {
						triangle_verts[3][j] += triangle_verts[i][j] / 3.0;
					}
				}
				lambda = phgGeomXYZ2Lambda(g,e,triangle_verts[3][0],triangle_verts[3][1],triangle_verts[3][2]);
				X = (FLOAT *)malloc(3 * sizeof(FLOAT));
				Y = (FLOAT *)malloc(3 * sizeof(FLOAT));
				for(i = 0; i < 3; i++) {
					X[i] = triangle_verts[i][(d+1)%3];
					Y[i] = triangle_verts[i][(d+2)%3];
				}
				triangle_area(X, Y, &area);
			}
			else if(type_flag == 2) {
				memset(tetrahetron_verts[4], 0, Dim * sizeof(FLOAT));
				for(i = 0; i < 2; i++) {
					for(j = 0; j < 2; j++) {
						linear_interpolation(g->verts[e->verts[p_tab[i]]],
								     g->verts[e->verts[n_tab[j]]],
								     tetrahetron_verts[2*i+j], value, d);
						for(k = 0; k < Dim; k ++) {
							tetrahetron_verts[4][k] += tetrahetron_verts[2*i+j][k] / 4.0;
						}
					}
				}
				lambda = phgGeomXYZ2Lambda(g,e,tetrahetron_verts[4][0],tetrahetron_verts[4][1],tetrahetron_verts[4][2]);
				X = (FLOAT *)malloc(4 * sizeof(FLOAT));
				Y = (FLOAT *)malloc(4 * sizeof(FLOAT));
				for(i = 0; i < 4; i++) {
					X[i] = tetrahetron_verts[i][(d+1)%3];
					Y[i] = tetrahetron_verts[i][(d+2)%3];
				}
				tetrahetron_area(X, Y, &area);
			}
			if(type_flag > 0 && type_flag < 4) {
				sum_area += area;
//if(fabs(lambda[0])+fabs(lambda[1])+fabs(lambda[2])+fabs(lambda[3]) > 1.0001) {
//printf("sum = %lf, z = %lf\n", fabs(lambda[0])+fabs(lambda[1])+fabs(lambda[2])+fabs(lambda[3]), value);
//printf("1=%lf, 2=%lf, 3=%lf, 4=%lf, sum=%lf\n", lambda[0],lambda[1],lambda[2],lambda[3],(fabs(lambda[0])+fabs(lambda[1])+fabs(lambda[2])+fabs(lambda[3])));
//phgError(1, "larger than 1.0");
//}
//printf("1=%lf, 2=%lf, 3=%lf, 4=%lf, sum=%lf\n", lambda[0],lambda[1],lambda[2],lambda[3],(fabs(lambda[0])+fabs(lambda[1])+fabs(lambda[2])+fabs(lambda[3])));
				phgDofEval(dof, e, lambda, &e_value);

//if(e_value < -0.00001) {
//printf("dof value smaller than 0\n");
//printf(" dof value = %lf, 1=%lf, 2=%lf, 3=%lf, 4=%lf, sum=%lf, z = %lf\n", e_value, lambda[0],lambda[1],lambda[2],lambda[3],(fabs(lambda[0])+fabs(lambda[1])+fabs(lambda[2])+fabs(lambda[3])), value);
//}

				dof_value += area * e_value;
			}
		}
	}

	MPI_Allreduce(&dof_value, &sum_value, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
	MPI_Allreduce(&sum_area, &total_area, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
/*
if(sum_value < 0) {
printf("dof sum value smaller than 0\n");
printf("z = %lf\n", value);
}
*/

	//print current section area
	return sum_value / total_area;
}

FLOAT Calculate_dof_positive_value(int region, int d, DOF *dof, FLOAT value) {
        SIMPLEX *e;
        GRID *g = dof->g;
        FLOAT area, volume, *X, *Y;
        //FLOAT *lambda, e_lambda[4] = {0.25, 0.25, 0.25, 0.25};
        FLOAT *lambda;
        FLOAT sum_area = 0;
        FLOAT total_area = 0;

        COORD triangle_verts[4];
        COORD tetrahetron_verts[5];

        int i, j, k, vert;
        int type_flag;//number of p_tab
        int p_tab[NVert];//save the tab of verts of a ELEMENT which has coordinate on d direction larger than value
        int n_tab[NVert];//if smaller than value

        FLOAT e_value, dof_value, sum_value;

        dof_value = 0;
        ForAllElements(g, e) {
                if(e->region_mark == region) {
                        memset(p_tab, 0, NVert*sizeof(int));
                        memset(n_tab, 0, NVert*sizeof(int));
                        type_flag = 0;
                        i = 0;
                        for(vert = 0; vert < NVert; vert++) {
                                if(g->verts[e->verts[vert]][d] >= value) {
                                        p_tab[type_flag++] = vert;
                                }
                                else {
                                        n_tab[i++] = vert;
                                }
                        }
                        if(type_flag == 1 || type_flag == 3) {
                                memset(triangle_verts[3], 0, Dim * sizeof(FLOAT));
                                for(i = 0; i < 3; i++) {
                                        if(type_flag == 1)
                                                linear_interpolation(g->verts[e->verts[p_tab[0]]],
                                                                     g->verts[e->verts[n_tab[i]]],
                                                                     triangle_verts[i], value, d);
                                        if(type_flag == 3)
                                                linear_interpolation(g->verts[e->verts[n_tab[0]]],
                                                                     g->verts[e->verts[p_tab[i]]],
                                                                     triangle_verts[i], value, d);
                                        for(j = 0; j < Dim; j ++) {
                                                triangle_verts[3][j] += triangle_verts[i][j] / 3.0;
                                        }
                                }
                                lambda = phgGeomXYZ2Lambda(g,e,triangle_verts[3][0],triangle_verts[3][1],triangle_verts[3][2]);
                                X = (FLOAT *)malloc(3 * sizeof(FLOAT));
                                Y = (FLOAT *)malloc(3 * sizeof(FLOAT));
                                for(i = 0; i < 3; i++) {
                                        X[i] = triangle_verts[i][(d+1)%3];
                                        Y[i] = triangle_verts[i][(d+2)%3];
                                }
                                triangle_area(X, Y, &area);
                        }
                        else if(type_flag == 2) {
                                memset(tetrahetron_verts[4], 0, Dim * sizeof(FLOAT));
                                for(i = 0; i < 2; i++) {
                                        for(j = 0; j < 2; j++) {
                                                linear_interpolation(g->verts[e->verts[p_tab[i]]],
                                                                     g->verts[e->verts[n_tab[j]]],
                                                                     tetrahetron_verts[2*i+j], value, d);
                                                for(k = 0; k < Dim; k ++) {
                                                        tetrahetron_verts[4][k] += tetrahetron_verts[2*i+j][k] / 4.0;
                                                }
                                        }
                                }
                                lambda = phgGeomXYZ2Lambda(g,e,tetrahetron_verts[4][0],tetrahetron_verts[4][1],tetrahetron_verts[4][2]);
                                X = (FLOAT *)malloc(4 * sizeof(FLOAT));
                                Y = (FLOAT *)malloc(4 * sizeof(FLOAT));
                                for(i = 0; i < 4; i++) {
                                        X[i] = tetrahetron_verts[i][(d+1)%3];
                                        Y[i] = tetrahetron_verts[i][(d+2)%3];
                                }
                                tetrahetron_area(X, Y, &area);
                        }
                        if(type_flag > 0 && type_flag < 4) {
                                sum_area += area;
                                phgDofEval(dof, e, lambda, &e_value);
				if(e_value < 0.0000) {
					e_value = 1e-6;
				}
                                dof_value += area * e_value;
                        }
                }
        }
        MPI_Allreduce(&dof_value, &sum_value, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
        MPI_Allreduce(&sum_area, &total_area, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);

        //print current section area
        return sum_value / total_area;
}


static void face_measure(BTYPE BDRY, GRID *g, INT *n, FLOAT *S) {
	DOF *I = phgDofNew(g, DOF_CONSTANT, 1, "I", DofNoAction);
	SIMPLEX *e;
	int i;
	*n = 0;
	*S = 0;
	phgDofSetDataByValue(I, 1.0);

	ForAllElements(g, e) {
		if(e->region_mark == 2) {
			for(i = 0; i < NFace; i++) {
				if(e->bound_type[i] & BDRY) {
					(*n)++;
					*S += phgQuadFaceDofDotDof(e, i, I, DOF_PROJ_NONE, I, 3);
				}
			}
		}
	}
	phgDofFree(&I);
}

static void region_measure(int region, GRID *g, INT *n, FLOAT *V) {
	DOF *I = phgDofNew(g, DOF_CONSTANT, 1, "I", DofNoAction);
	SIMPLEX *e;
	*n = 0;
	*V = 0;
	ForAllElements(g, e) {
		if(e->region_mark == region) {
			(*n)++;
			*V += phgGeomGetVolume(g, e);
		}
	}
	phgDofFree(&I);
}

static void bdry_quad(DOF *p, FLOAT *value) {
	DOF *I;
	GRID *g = p->g;
	SIMPLEX *e;
	FLOAT sum = 0;
	int i;

	I = phgDofNew(g, DOF_CONSTANT, 1, "I", DofNoAction);
	phgDofSetDataByValue(I, 1.0);

	ForAllElements(g, e) {
		if(e->region_mark == solvent_region) {
			for(i = 0; i < NFace; i++) {
				if(e->bound_type[i] & reac_bdry) {
					FLOAT tmp;
					tmp = phgQuadFaceDofDotDof(e, i, p, DOF_PROJ_NONE, I, QUAD_DEFAULT);
					sum += tmp;
				}
			}
		}
	}
	phgDofFree(&I);
	*value = sum;
}

static void Estimate_Poisson_error(DOF *u, DOF **p, DOF **indicator, DOF **error) {
	GRID *g = u->g;
	SIMPLEX *e;
	DOF *jump, *residual, *grad_u, *tmp;

	residual = phgDofCopy(u, NULL, NULL, NULL);
	INT i, N = DofGetDataCount(u);
	for(i = 0; i < N; i++) {
		residual->data[i] *= Epsilon_s / Kappa_2 * Pow(change_unit, 2) / NA;
	}
	grad_u = phgDofGradient(residual, NULL, NULL, NULL);
	jump = phgQuadFaceJump(grad_u, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
	tmp = phgDofDivergence(grad_u, NULL, NULL, NULL);
	phgDofFree(&grad_u);
	phgDofCopy(tmp, &residual, NULL, NULL);
	phgDofFree(&tmp);

//	phgDofAXPBY(Epsilon_s / Kappa_2 * Pow(change_unit, 2), residual, 0.0, &tmp);
//	phgDofAXPBY(1.0, tmp, 0.0, &residual);
//
	for(i = 0; i < NION; i++) {
		phgDofAXPY(ion[i].Z / NA, p[i], &residual);
	}

	//phgPrintf("\nCreating \"%s\".\n", phgExportVTK(g, "PNP_laplace_u.vtk", residual, NULL));
	//phgPrintf("\n%f!!!!!!!!!\n", phgDofNormL2(residual));

	ForAllElements(g, e) {
//		if(e->region_mark == solute_region) {
                if(e->region_mark == solute_region || e->region_mark == membrane_region) {
			*DofElementData((*indicator), e->index) = 0.0;
			*DofElementData((*error), e->index) = 0.0;
		}
		else {
			//int N = DofGetNBas(u, e);
			FLOAT Epsilon = 0.0;
			FLOAT Kappa = 0.0;
			get_Kappa(e, &Epsilon, &Kappa);
			FLOAT h, eta = 0.0;
			FLOAT diam = phgGeomGetDiameter(g, e);
			e->mark = 0;
			for(i = 0; i < NFace; i++) {
				if (!(e->bound_type[i] & (DIRICHLET | bmlc_bdry))) {
					h = phgGeomGetFaceDiameter(g, e, i);
					eta += *DofFaceData(jump, e->faces[i]) * h;
				}
			}
			eta = eta * 0.5 + diam * diam * phgQuadDofDotDof(e, residual, residual, QUAD_DEFAULT);
			*DofElementData((*error), e->index) = (diam * diam * phgQuadDofDotDof(e, residual, residual, QUAD_DEFAULT));
			*DofElementData((*indicator), e->index) = (eta);
//printf("\n%le", *DofElementData((*indicator), e->index));
		}
	}
	phgDofFree(&jump);
	phgDofFree(&residual);
}

static void Estimate_N_P_error(INT TD, DOF *u, DOF *TMP_u, DOF **p, DOF **TMP_p, DOF **indicator, DOF **error) {
	GRID *g = u->g;
	SIMPLEX *e;
	INT i, N = DofGetDataCount(u);

	DOF *jump, *residual, *lhs, *rhs;
	DOF *flux, *TMP_flux, *grad_p, *TMP_grad_p, *grad_u, *TMP_grad_u;

	grad_p = phgDofGradient(p[which_ion], NULL, NULL, NULL);
	grad_u = phgDofGradient(u, NULL, NULL, NULL);
	dof_coef_vec(&grad_u, p[which_ion], grad_u);
	flux = phgDofCopy(grad_p, NULL, NULL, NULL);
	phgDofAXPBY(ion[which_ion].Z * ion[which_ion].D / NA, grad_u, ion[which_ion].D / NA, &flux);
	phgDofFree(&grad_p);
	phgDofFree(&grad_u);

	lhs = phgDofNew(g, DOF_DEFAULT, 1, "tmp", DofInterpolation);
	phgDofSetDataByValue(lhs, 0.0);

	jump = phgQuadFaceJump(flux, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
	rhs = phgDofDivergence(flux, NULL, NULL, NULL);
	residual = phgDofCopy(u, NULL, NULL, NULL);
	phgDofCopy(rhs, &residual, NULL, NULL);
	phgDofAXPY(-1.0, lhs, &residual);
	phgDofFree(&flux);
	phgDofFree(&lhs);
	phgDofFree(&rhs);


	ForAllElements(g, e) {
//		if(e->region_mark == solute_region) {
                if(e->region_mark == solute_region || e->region_mark == membrane_region) {
			*DofElementData((*indicator), e->index) = 0.0;
			*DofElementData((*error), e->index) = 0.0;
		}
		else {
			//int N = DofGetNBas(u, e);
			FLOAT Epsilon = 0.0;
			FLOAT Kappa = 0.0;
			get_Kappa(e, &Epsilon, &Kappa);
			FLOAT h, eta = 0.0;
			FLOAT diam = phgGeomGetDiameter(g, e);
			e->mark = 0;
			for(i = 0; i < NFace; i++) {
				if (!(e->bound_type[i] & (DIRICHLET | bmlc_bdry))) {
					h = phgGeomGetFaceDiameter(g, e, i);
					eta += *DofFaceData(jump, e->faces[i]) * h;
				}
			}
			eta = eta * 0.5 + diam * diam * phgQuadDofDotDof(e, residual, residual, QUAD_DEFAULT);
			*DofElementData((*error), e->index) = (diam * diam * phgQuadDofDotDof(e, residual, residual, QUAD_DEFAULT));
			*DofElementData((*indicator), e->index) = (eta);
//printf("\n%le", *DofElementData((*indicator), e->index));
		}
	}
	phgDofFree(&jump);
	phgDofFree(&residual);
}


/* For Calculate current */
//calculate the coordinate of the point c on the line a,b with c[d] = value
//lambda * a + (1 - lambda) * b = c
//lambda * a[d] + (1 - lambda) * b[2] = z
//lambda = (z - b[d]) / (a[d] - b[d])
//1 - lambda = (a[d] - z) / (a[d] - b[d])
//d short for direction d = 1, 2, ..., Dim;
static void linear_interpolation(COORD a, COORD b, COORD c, FLOAT value, int d) {
	if(a[d] == b[d]) {
		phgError(1, "\ninput points have same direction value");
	}
	else {
		int i;
//		c = (COORD *)malloc(sizeof(COORD));
		c[d] = value;
		for(i = 1; i < Dim; i++) {
			c[(d+i)%Dim] = (value - b[d]) / (a[d] - b[d]) * a[(d+i)%Dim]
				     + (value - a[d]) / (b[d] - a[d]) * b[(d+i)%Dim];
		}
	}
}

// (X, Y) to save all pairs of triangle vert (x[i], y[i]) i = 1, 2, 3, with X, Y belong to R^3
static void triangle_area(FLOAT *X, FLOAT *Y, FLOAT *area) {
	int i;
	*area = 0;
	for(i = 0; i < 3; i++) {
		*area += (X[i] * Y[(i+1)%3] - X[(i+1)%3] * Y[i]);
	}
	*area = Fabs(*area) * 0.5;
	free(X);
	free(Y);
}

// (X, Y) to save all pairs of tetrahetron vert (x[i], y[i]) i = 1, 2, 3, 4, with X, Y belong to R^4
static void tetrahetron_area(FLOAT *X, FLOAT *Y, FLOAT *area) {
	FLOAT *child_area = (FLOAT *)malloc(4 * sizeof(FLOAT));
	FLOAT *x, *y;
	int i, j;
	*area = 0;
	for(i = 0; i < 4; i++) {
		x = (FLOAT *)malloc(3 * sizeof(FLOAT));
		y = (FLOAT *)malloc(3 * sizeof(FLOAT));
		for(j = 0; j < 3; j++) {
			x[j] = X[(i+j)%4];
			y[j] = Y[(i+j)%4];
		}
		triangle_area(x, y, &child_area[i]);
		*area += child_area[i] * 0.5;
	}
}

static void dof_save(char *fn, int *step_n, DOF *dof, ...) {
	INT i, N = DofGetDataCount(dof);
	MPI_Comm_rank(MPI_COMM_WORLD, &i);
	va_list ap;
	int tmp[3];
	char dir_name[40] = "./PNP_dof_saves/";
	strcat(dir_name, fn);
	char filename[20] = "_000_dof.tmp";
	tmp[0] = i / 100;
	tmp[1] = (i / 10) % 10;
	tmp[2] = i % 10;
	filename[1] += tmp[0];
	filename[2] += tmp[1];
	filename[3] += tmp[2];
	strcat(dir_name, filename);

	if(*step_n >= 0) {
		FILE *file = fopen(dir_name, "wb");
		if(file == NULL) phgError(1, "cannot open %s", dir_name);
		fwrite(step_n, sizeof(int), 1, file);
		va_start(ap, dof);
		for(i = 0; i < 32; i++) {
			if(dof == NULL) break;
			fwrite(dof->data, sizeof(FLOAT), N, file);
			dof = va_arg(ap, DOF *);
		}
		va_end(ap);
		fclose(file);
		phgPrintf("time dependent dof saved in PNP_dof_saves/%s\n", fn);
	}
	else {
		FILE *file = fopen(dir_name, "w");
		if(file == NULL) phgError(1, "cannot open %s", dir_name);
		va_start(ap, dof);
		for(i = 0; i < 32; i++) {
			if(dof == NULL) break;
			int n;
			for(n = 0; n < N; n++)
			fprintf(file, "%lf\n", dof->data[n]);
			fprintf(file, "\n");
			dof = va_arg(ap, DOF *);
		}
		va_end(ap);
		fclose(file);
		phgPrintf("\nsteady-state dof loaded in %s_N_dof.data", fn);
	}
}

static void dof_load(char *fn, int *step_n, DOF *dof, ...) {
	GRID *g = dof->g;
	INT i, N = DofGetDataCount(dof);
	MPI_Comm_rank(MPI_COMM_WORLD, &i);
	va_list ap;
	char dir_name[40] = "./PNP_dof_saves/";
	strcat(dir_name, fn);
	char filename[20] = "_000_dof.tmp";
	filename[1] += i / 100;
	filename[2] += (i / 10) % 10;
	filename[3] += i % 10;
	strcat(dir_name, filename);

	if(*step_n >= 0) {
		FILE *file = fopen(dir_name, "rb");
		if(file == NULL) phgError(1, "cannot open %s", dir_name);
		fread(step_n, sizeof(int), 1, file);
		va_start(ap, dof);
		for(i = 0; i < 32; i++) {
			if(dof == NULL) break;
			fread(dof->data, sizeof(FLOAT), N, file);
			dof = va_arg(ap, DOF *);
		}
		va_end(ap);
		fclose(file);
		phgPrintf("time dependent dof loaded in PNP_dof_saves/%s\n", fn);
	}
	else {
		FILE *file = fopen(dir_name, "r");
		if(file == NULL) phgError(1, "cannot open %s", dir_name);
		va_start(ap, dof);
		for(i = 0; i < 32; i++) {
			if(dof == NULL) break;
			int n;
			for(n = 0; n < N; n++)
			fscanf(file, "%lf\n", &dof->data[n]);
			dof = va_arg(ap, DOF *);
		}
		va_end(ap);
		fclose(file);
		phgPrintf("\nsteady-state dof loaded in ./tmp/%s_N_dof.data", fn);
	}
}

FLOAT P_B_Dof_Int(int region, DOF *U) {
	int i, j, nvalues;
	FLOAT d, d0, tmp, sum;
	const FLOAT *a, *b;
	QUAD *quad;
	ELEMENT *e;
	INT N = DofGetDataCount(U);

	DOF *u = phgDofCopy(U, NULL, NULL, NULL);
	for(i = 0; i < N; i++) {
		u->data[i] = Exp(-ion[which_ion].Z * u->data[i]) - 1;
	}

	sum = 0.;
	ForAllElements(u->g, e){
		if(e->region_mark == region) {
			i = DofTypeOrder(u, e);
			if (i < 0)
				i = 1;
			quad = phgQuadGetQuad3D(i);
			d = 0.;
			b = quad->weights;
			a = phgQuadGetDofValues(e, u, quad);
			for (i = 0; i < quad->npoints; i++) {
				tmp = 0.;
				for (j = 0; j < DofDim(u); j++) {
					tmp += *(a++);
				}
				d += tmp * *(b++);
			}
			sum += d * phgGeomGetVolume(u->g, e);
		}
	}
#if USE_MPI
	MPI_Allreduce(&sum, &d0, 1, PHG_MPI_FLOAT, PHG_SUM, u->g->comm);
#else
	d0 = sum;
#endif
	d0 *= ion[which_ion].c * NA / Pow(change_mesh_unit, 3);
	phgDofFree(&u);
	return d0;
}


//-------Bernoulli function --rgshen2021-4-19
FLOAT Ber(FLOAT x) {
    FLOAT y;
    if(Fabs(x) > 1.0e-4)
        y = x/(Exp(x) - 1.0);
    else
        y = ((-1./720 * x * x + 1./12) * x -1./2) * x + 1.;
    return y;
}




//------参考书写 SMPNP 方程的 通道电流计算---2021-6-12  调试中-----暂时不对 6-13----
           //SMPNP_Calculate_current(solvent_region, current_direction, u, u_D, p, section_value, size_c, Epsilon);  //rgshne 调试中？？6-12
static void SMPNP_Calculate_current(int region, int d, DOF *u, DOF *u_D, DOF **p, FLOAT value, DOF **size_c, DOF *Epsilon) {
        FLOAT sum_current = 0;
        FLOAT p_current[NION];
        FLOAT current[NION];
        FLOAT sum_current1 = 0;
        FLOAT sum_current2 = 0;
        FLOAT sum_current3 = 0;
        FLOAT sum_current4 = 0;
        FLOAT p_current1[NION];
        FLOAT p_current2[NION];
        FLOAT p_current3[NION];
        FLOAT p_current4[NION];
        FLOAT current1[NION];
        FLOAT current2[NION];
        FLOAT current3[NION];
        FLOAT current4[NION];  //for born correct

        SIMPLEX *e;
        GRID *g = u->g;
        FLOAT area, volume, *X, *Y;
        FLOAT *lambda, e_lambda[4] = {0.25, 0.25, 0.25, 0.25};
        FLOAT sum_area = 0;
        FLOAT total_area = 0;

        COORD triangle_verts[4];
        COORD tetrahetron_verts[5];

        DOF *tmp, *tmp_eps, *E, *G_u_D_G_v, *div_grad_eps;
	    DOF *tmp_tmp, *multiple;
        DOF *TMP_u_D = phgDofCopy(u_D, NULL, NULL, NULL);
        DOF *tmp_u_D = phgDofCopy(u_D, NULL, NULL, NULL);
        tmp = phgDofNew(g, DOF_DEFAULT, 1, "tmp", DofInterpolation);
        tmp_eps = phgDofNew(g, DOF_DEFAULT, 1, "tmp_eps", DofInterpolation);
        if(!strcmp(u->type->name, "P1")) {
			G_u_D_G_v = phgDofNew(g, DOF_DG0, 1, "G_u_D_G_v", DofInterpolation);
	        div_grad_eps = phgDofNew(g, DOF_DG0, 1, "div_grad_eps", DofInterpolation);
               	tmp_tmp = phgDofNew(u->g, DOF_DG1, 1, "tmp_tmp", DofInterpolation);
               	multiple = phgDofNew(u->g, DOF_DG1, 1, "tmp_grad_u_2", DofInterpolation);
		}
        if(!strcmp(u->type->name, "P2")) {
			G_u_D_G_v = phgDofNew(g, DOF_DG1, 1, "G_u_D_G_v", DofInterpolation);
        	div_grad_eps = phgDofNew(g, DOF_DG0, 1, "div_grad_eps", DofInterpolation);
			if(born_correct) {
				tmp_tmp = phgDofNew(u->g, DOF_DG1, 1, "tmp_tmp", DofInterpolation);
                multiple = phgDofNew(u->g, DOF_DG1, 1, "tmp_grad_u_2", DofInterpolation);
			}
			else {
                tmp_tmp = phgDofNew(u->g, DOF_DG2, 1, "tmp_tmp", DofInterpolation);
                multiple = phgDofNew(u->g, DOF_DG2, 1, "tmp_grad_u_2", DofInterpolation);
			}
		}
        if(!strcmp(u->type->name, "P3")) {
			G_u_D_G_v = phgDofNew(g, DOF_DG2, 1, "G_u_D_G_v", DofInterpolation);
        	div_grad_eps = phgDofNew(g, DOF_DG1, 1, "div_grad_eps", DofInterpolation);
		    if(born_correct) {
                tmp_tmp = phgDofNew(u->g, DOF_DG2, 1, "tmp_tmp", DofInterpolation);
                multiple = phgDofNew(u->g, DOF_DG2, 1, "tmp_grad_u_2", DofInterpolation);
			}
			else {
                tmp_tmp = phgDofNew(u->g, DOF_DG3, 1, "tmp_tmp", DofInterpolation);
                multiple = phgDofNew(u->g, DOF_DG3, 1, "tmp_grad_u_2", DofInterpolation);
			}
		}

        phgDofSetDataByValue(tmp, 0.0);
        phgDofSetDataByValue(tmp_eps, 0.0);
        phgDofSetDataByValue(G_u_D_G_v, 0.0);
        phgDofSetDataByValue(div_grad_eps, 0.0);
        phgDofSetDataByValue(tmp_tmp, 0.0);
        phgDofSetDataByValue(multiple, 0.0);
        phgDofAXPBY(1.0, u, -1, &tmp_u_D);    //tmp_u_D = u - u_D
        phgDofAXPBY(1.0, u, -0.5, &TMP_u_D);  //TMP_u_D = u - 0.5*u_D
        
        int i, j, k, l, vert;
        int type_flag;//number of p_tab
        int p_tab[NVert];//save the tab of verts of a ELEMENT which has coordinate on d direction larger than value
        int n_tab[NVert];//if smaller than value
        // phgPrintf("\n%d", type_flag);
        FLOAT grad_u[Dim];
        FLOAT grad_p[NION][Dim];
        FLOAT grad_eps[Dim];
        FLOAT p_value[NION];
        FLOAT grad_G2_value[Dim];
        FLOAT tmp_value[1];
        FLOAT multiple_value[Dim];
		
		FLOAT size_c_value[NION]; //rgshen2021-6-13

        memset(current, 0, (NION)*sizeof(FLOAT));
        memset(current1, 0, (NION)*sizeof(FLOAT));
        memset(current2, 0, (NION)*sizeof(FLOAT));
        memset(current3, 0, (NION)*sizeof(FLOAT));
        memset(current4, 0, (NION)*sizeof(FLOAT));
        memset(p_current, 0, (NION)*sizeof(FLOAT));
        memset(p_current1, 0, (NION)*sizeof(FLOAT));
        memset(p_current2, 0, (NION)*sizeof(FLOAT));
        memset(p_current3, 0, (NION)*sizeof(FLOAT));
        memset(p_current4, 0, (NION)*sizeof(FLOAT));
        //phgPrintf("\n%d", type_flag);
        ForAllElements(g, e) {
			if(e->region_mark == region) {
                memset(p_tab, 0, NVert*sizeof(int));
                memset(n_tab, 0, NVert*sizeof(int));
                type_flag = 0;
                i = 0;
                //phgPrintf("@@@@@@%\n%d", NVert);
                for(vert = 0; vert < NVert; vert++) {
                    if(g->verts[e->verts[vert]][d] >= value) {
						p_tab[type_flag++] = vert;
                    }
                    else {
                        n_tab[i++] = vert;
                    }
                    //phgPrintf("\n%d", type_flag);
                }
                //phgPrintf("\n%d", type_flag);

                if(type_flag == 1 || type_flag == 3) {
                    //COORD *triangle_verts = (COORD *)malloc(4 * sizeof(COORD));
                    memset(triangle_verts[3], 0, Dim * sizeof(FLOAT));
                    for(i = 0; i < 3; i++) {
                        if(type_flag == 1)
                            linear_interpolation(g->verts[e->verts[p_tab[0]]],
                                                 g->verts[e->verts[n_tab[i]]],
                                                 triangle_verts[i], value, d);
                        if(type_flag == 3)
                            linear_interpolation(g->verts[e->verts[n_tab[0]]],
                                                 g->verts[e->verts[p_tab[i]]],
                                                 triangle_verts[i], value, d);
                        for(j = 0; j < Dim; j ++) {
                            triangle_verts[3][j] += triangle_verts[i][j] / 3.0;
                        phgPrintf("****, %e", triangle_verts[3][j]);
                        }
                    }
                    lambda = phgGeomXYZ2Lambda(g,e,triangle_verts[3][0],triangle_verts[3][1],triangle_verts[3][2]);
                    X = (FLOAT *)malloc(3 * sizeof(FLOAT));
                    Y = (FLOAT *)malloc(3 * sizeof(FLOAT));
                    for(i = 0; i < 3; i++) {
                        X[i] = triangle_verts[i][(d+1)%3];
                        Y[i] = triangle_verts[i][(d+2)%3];
                    }
                    triangle_area(X, Y, &area);
                    //free(triangle_verts);
                    //phgDofEval(flux, e, e_lambda, &value);
                    //current += area * value;
                }
                else if(type_flag == 2) {
                    //COORD *tetrahetron_verts = (COORD *)malloc(5 * sizeof(COORD));
                    memset(tetrahetron_verts[4], 0, Dim * sizeof(FLOAT));
                    for(i = 0; i < 2; i++) {
                        for(j = 0; j < 2; j++) {
                            linear_interpolation(g->verts[e->verts[p_tab[i]]],
                                                 g->verts[e->verts[n_tab[j]]],
                                                 tetrahetron_verts[2*i+j], value, d);
                            for(k = 0; k < Dim; k ++) {
                                tetrahetron_verts[4][k] += tetrahetron_verts[2*i+j][k] / 4.0;
                            }
                        }
                    }
                    lambda = phgGeomXYZ2Lambda(g,e,tetrahetron_verts[4][0],tetrahetron_verts[4][1],tetrahetron_verts[4][2]);
                    X = (FLOAT *)malloc(4 * sizeof(FLOAT));
                    Y = (FLOAT *)malloc(4 * sizeof(FLOAT));
                    for(i = 0; i < 4; i++) {
                        X[i] = tetrahetron_verts[i][(d+1)%3];
                        Y[i] = tetrahetron_verts[i][(d+2)%3];
                    }
                    tetrahetron_area(X, Y, &area);
                    //free(tetrahetron_verts);
                    //phgPrintf("\n%e",area);
                    //phgDofEval(flux, e, e_lambda, &value);
                    //current += area * value;
                }
                if(type_flag > 0 && type_flag < 4) {
					//-----------------------------------------------------------先暂时注释掉----------
                    /* //if(MOD == 0) { // 经验证，这个就是 PNP 模拟通道电流,; 参考它书写smpnp 电流计算，增加尺寸效应部分即可
                        sum_area += area;
                        //phgDofEvalGradient(u, e, lambda, NULL, grad_u);
                        phgDofEvalGradient(TMP_u_D, e, lambda, NULL, grad_u);
                        for(i = 0; i < NION; i++) {
                            phgDofEval(p[i], e, lambda, p_value + i);
                            phgDofEvalGradient(p[i], e, lambda, NULL, grad_p[i]);
                        }
                        for(i = 0; i < NION; i++) {
                            if(current_switch) {
                            	current[i] += area * (grad_p[i][d] + ion[i].Z * p_value[i] * grad_u[d]);
                            	current1[i] += area * grad_p[i][d];
                            	current2[i] += area * ion[i].Z * p_value[i] * grad_u[d];
                            }
                            else {
                                if(type_flag != 2) {
                                	current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL);
                                    current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL);
                                }
                                else {
                                    for(j = 0; j < 4; j++) {
                	                    for(k = 0; k < 3; k++) {
        	                                for(l = 0; l < 3; l++) {
	                                            triangle_verts[k][l] = tetrahetron_verts[(j+k)%4][l];
                                            }
                                        }
                                        current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL) / 2.0;
                                        current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL) / 2.0;
                                    }
                                }
                                current[i] = current1[i] + current2[i];
                            }
                        }
                    //} */
					//-------------------------------------------------------------------
					 //---------smpnp 电流计算， 在 pnp 基础上 增加尺寸效应部分即可
                        sum_area += area;
                        //phgDofEvalGradient(u, e, lambda, NULL, grad_u);
                        phgDofEvalGradient(TMP_u_D, e, lambda, NULL, grad_u);
                        for(i = 0; i < NION; i++) {
                            phgDofEval(p[i], e, lambda, p_value + i);
                            phgDofEvalGradient(p[i], e, lambda, NULL, grad_p[i]);
							
							phgDofEval(size_c[i], e, lambda, size_c_value + i); //rgshen2021-6-13
                        }
                        for(i = 0; i < NION; i++) {
                            if(current_switch) {
                            	//current[i] += area * (grad_p[i][d] + ion[i].Z * p_value[i] * grad_u[d]);
								current[i] += area * (grad_p[i][d] + size_c_value[i] * p_value[i] * (Pow(ion[0].a, 3) * grad_p[0][d] + Pow(ion[1].a, 3) * grad_p[1][d]) 
								                     + ion[i].Z * p_value[i] * grad_u[d]);
								
                            	//current1[i] += area * grad_p[i][d];
								current1[i] += area * grad_p[i][d] + size_c_value[i] * p_value[i] * (Pow(ion[0].a, 3) * grad_p[0][d] + Pow(ion[1].a, 3) * grad_p[1][d]) ;
								
                            	current2[i] += area * ion[i].Z * p_value[i] * grad_u[d];
                            }
                            else {
                                if(type_flag != 2) {
                                	current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL);  // smpnp 时 这里的 p[i] 要改动吗？？？
                                    current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL);
                                }
                                else {
                                    for(j = 0; j < 4; j++) {
                	                    for(k = 0; k < 3; k++) {
        	                                for(l = 0; l < 3; l++) {
	                                            triangle_verts[k][l] = tetrahetron_verts[(j+k)%4][l];
                                            }
                                        }
                                        current1[i] += PNP_quad_channel_section(e, 3, triangle_verts, 1, p[i], NULL) / 2.0;
                                        current2[i] += ion[i].Z * PNP_quad_channel_section(e, 3, triangle_verts, 1, TMP_u_D, p[i], NULL) / 2.0;
                                    }
                                }
                                current[i] = current1[i] + current2[i];
                            }
                        }
                    //------------------------------------------- 	
                }
            }
        }

        phgPrintf("\n\nthe section value is %f A\n", value);
        for(i = 0; i < NION; i++) {
            //p: N/cm^3 or N/m^3 -> N/A^3 | D: A^2/ps -> A^2/s | I: A ->pA | grad and quad: mesh_unit -> A
            current[i] *= -ion[i].D * chan_ratio * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;
            current1[i] *= -ion[i].D * chan_ratio * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;
            current2[i] *= -ion[i].D * chan_ratio * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;
            current3[i] *= -ion[i].D * chan_ratio * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;
            current4[i] *= -ion[i].D * chan_ratio * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;
            //current[i] *= -ion[i].D * ion[i].Z / Pow(change_unit, 3) / Pow(mesh_unit, 2) * 1.0e12 * 1.0e12 * ec;

            //phgPrintf("*********8%f\n",current[i]);
            MPI_Allreduce(&current[i], &p_current[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            MPI_Allreduce(&current1[i], &p_current1[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            MPI_Allreduce(&current2[i], &p_current2[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            MPI_Allreduce(&current3[i], &p_current3[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            MPI_Allreduce(&current4[i], &p_current4[i], 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            MPI_Allreduce(&sum_area, &total_area, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
            sum_current += p_current[i];
            sum_current1 += p_current1[i];
            sum_current2 += p_current2[i];
            sum_current3 += p_current3[i];
            sum_current4 += p_current4[i];
            phgPrintf("\nion[%d]    current is %e pA", i, p_current[i]);
            phgPrintf("\nion[%d]    the first part of current is %e pA", i, p_current1[i]);
            phgPrintf("\nion[%d]    the second part of current is %e pA", i, p_current2[i]);
            //phgPrintf("\nion[%d]    the third part of current is %e pA", i, p_current3[i]);
            //phgPrintf("\nion[%d]    the forth part of current is %e pA", i, p_current4[i]);
        }
        //print current section area
        phgPrintf("\n\ntotal    current is %e pA\n", sum_current);  // 这是什么意思呢？rgshen
        phgPrintf("\n\nthe first part    current is %e pA\n", sum_current1);
        phgPrintf("\n\nthe second part    current is %e pA\n", sum_current2);
        //phgPrintf("\n\nthe third part    current is %e pA\n", sum_current3);
        //phgPrintf("\n\nthe forth part    current is %e pA\n", sum_current4);
        phgPrintf("\n\ntotal area is %e A^2\n", total_area);
        if(phgRank == 1 && out_current_file != NULL) {
            FILE *f_in;
            f_in = fopen(out_current_file, "a");
            if(f_in == NULL) {
                phgError(1, "Can not open current output file %s.\n", out_current_file);
            }
            //fprintf(f_in, "%f   %le   %le   %le   %le   %le\n", value, sum_current, sum_current1, sum_current2, sum_current3, sum_current4);
//----------------------------------
#if 1
			fprintf(f_in, "%d DOF, %d elements, %d submeshes, load imbalance: %lg\n",
			        DofGetDataCountGlobal(u), g->nleaf_global, g->nprocs, (double)g->lif);//rgshen2020-11-12 这两处打印生成一次即注释掉
			fprintf(f_in, "u_bulk(V)	total_current(pA)	the_first_part(pA)	the_second_part(pA)\n"); //rgshen2020-11-12
#endif
			fprintf(f_in, "%f   %le       %le       %le\n", u_bulk, sum_current, sum_current1, sum_current2);//rgshen2020-11-12
//--------------------------------------------
            fclose(f_in);
        }
        //phgDofFree(&TMP);
        phgDofFree(&G_u_D_G_v);
        phgDofFree(&tmp_u_D);
        phgDofFree(&TMP_u_D);
        phgDofFree(&tmp);
	    phgDofFree(&tmp_eps);
        phgDofFree(&multiple);
       	phgDofFree(&tmp_tmp);
	    phgDofFree(&div_grad_eps);
	    
}









