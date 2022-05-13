#ifndef PNP_FUNC_H
#define PNP_FUNC_H

#include"phg.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

//if a function is not listed here, that means it remains to be test

/*functions in PNP_build_solver.c*/

static void Build_H_Solver(SOLVER *solver, DOF *H);

static void Build_Primitive_Solver(SOLVER* solver, DOF *u, DOF *size_p, DOF **p, DOF *grad_G, DOF *grad_H);

static void Build_Primitive_Solver_u_D(SOLVER* solver, DOF *u, DOF *size_p, DOF **p, DOF *grad_G, DOF *grad_H);

static void VD_Build_Primitive_Solver(SOLVER* solver, DOF *u, DOF *size_p, DOF **p, DOF *grad_G, DOF *grad_H, DOF *Epsilon);

static void VD_Build_Primitive_Solver_u_D(SOLVER* solver, DOF *u, DOF *size_p, DOF **p, DOF *grad_G, DOF *grad_H, DOF *Epsilon);

static void Build_p_Solver(SOLVER *solver, DOF *u, DOF * size_p, DOF **p, DOF *D);

static void Build_p_Solver_born(SOLVER *solver, DOF *u, DOF *size_p, DOF **p, DOF *D, DOF *Epsilon);

static void Build_p_Solver_u_D(SOLVER *solver, DOF *u, DOF *u_D, DOF *size_p, DOF **p, DOF *D);

//-------下面还在测试�?- rgshen ----
static void Benchmark_Build_p_Solver(SOLVER *solver, DOF *u, DOF *size_p, DOF **p, DOF *D);
static void Benchmark_Build_Solver(SOLVER* solver, DOF *u, DOF *size_p, DOF **p, DOF *grad_G, DOF *grad_H);

static void get_size_u_p(DOF **tmp, DOF *u, DOF **p); //rgshen2021-4-26
static void SMPNP_NP_Slotboomtransform_averaging_EAFE(SOLVER *solver, DOF *u,  DOF **p, DOF *D); 

static void get_size_c(DOF **tmp, DOF **p); //rgshen2021-6-12
static void SMPNP_NP_Slover_Primitive(SOLVER *solver, DOF *u, DOF **size_c, DOF **p, DOF *D, DOF *size_p);

static void myquad(QUAD *quad, SIMPLEX *e, int i, int j, DOF *p, FLOAT *DATA2, FLOAT *DATA5, FLOAT **DATAc, FLOAT **DATAGc, FLOAT h, FLOAT *value);
static void PNP_PRIMITIVE_NP(SOLVER *solver, DOF *u, DOF *p, DOF **c, DOF *d);  //这是乔瑜 标准FEM 求解smpnp-->
//------------------------------------------------------------------------------------------

/*functions in PNP_more.c*/

typedef void PNP_FUNC_1(FLOAT *);

typedef void PNP_FUNC_2(FLOAT *, FLOAT *);

PNP_FUNC_1 func_id, func_unit, func_exp_1, func_exp_2, func_exp_3, func_exp_4;

PNP_FUNC_2 func_p_k;

static void dof_positive(DOF *p);//if p->data < 0 set p->data = 0

static void dof_restriction(DOF *p, DOF *P);//if p->data > P->data set p->data = P->data

static void dof_direction(DOF **tmp, DOF *u, int d);//choose vector u one direction

static void dof_const(DOF **tmp, DOF *u, FLOAT k);//tmp = k * u

static void dof_coef_vec(DOF **tmp, DOF *u, DOF *v);//tmp[3] = u * v[3]

static void dof_pow(DOF **tmp, DOF *dof, FLOAT a);//tmp = pow( dof, a )

static void dof_plus(DOF **tmp, DOF *dof, ...);//tmp = dof1 + dof2 + ...

static void dof_multiple(DOF **tmp, DOF *dof, ...);//tmp = dof1 * dof2 * ...

static void Region_Quad_Dof_G_D_Bas_G_B(int region, SIMPLEX *e, DOF *dof, DOF *dof_vec, DOF *u, int i, DOF *v, int j, int order, FLOAT *value );

static void Region_Quad_Bas_Dof_G_B(int region, SIMPLEX *e, DOF *dof, DOF *u, int i, DOF *v, int j, int order, FLOAT *value);
//integrate(u_i * dof dot grad_v_j) in e and region with order

static void Region_Quad_Dof_G_B_G_B(int region, SIMPLEX *e, DOF *dof, DOF *u, int i, DOF *v, int j, int order, FLOAT *value);
//integrate(dof * grad_u_i dot grad_v_j) in e and region with order

static void Region_Quad_Dof_Bas_Bas(int region, SIMPLEX *e, DOF *dof, DOF *u, int i, DOF *v, int j, int order, FLOAT *value);
//integrate(dof * u_i * v_j) in e and region with order

static void Region_Quad_SUPG_term_Primitive(int region, SIMPLEX *e, DOF *grad_u, DOF *p, int i, int j, int order, FLOAT *value);

static void Region_Quad_SUPG_term_Slotboom(int region, SIMPLEX *e, DOF *grad_u, DOF *u, DOF *size_p, DOF *P, int i, int j, int order, FLOAT *value);

static void get_exp_u(DOF **tmp, DOF *u, PNP_FUNC_1 func);

static void get_size_p(DOF **tmp, DOF **p, PNP_FUNC_2 func);

static void p_Transform(FLOAT Alpha, DOF *p, DOF *P, DOF *u, DOF *size_p);

static void P_Transform(FLOAT Alpha, DOF *p, DOF *P, DOF *u, DOF *size_p);

static void t_Transform(FLOAT Alpha, DOF *p, DOF *T);

static void mark_region(GRID *g, int region, int level);

void get_Epsilon_all(DOF **Full_Epsilon, DOF *Epsilon, DOF **p);

static void get_Epsilon_all_test(DOF **Full_Epsilon, DOF *Epsilon);

void get_Epsilon(DOF **tmp, DOF **p);

void get_Epsilon_r(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

void get_G_u_DOT_G_v(DOF **tmp, DOF *u, DOF *v);

void get_G_D_DOT_G_D(DOF **tmp, DOF *dof);

static void See_bdry(DOF **bdry);

static void Calculate_energy(DOF *u, DOF *H);

static void Calculate_chemical_potential(DOF **tmp, DOF *u, DOF **p, DOF *Epsilon, int j);

static void Calculate_potential_energy(DOF **tmp, DOF *u, DOF *Epsilon, int j);

FLOAT PNP_quad_channel_section(ELEMENT *e, FLOAT npoint, COORD *points, INT grad_dofs, DOF *dof, ...);

static void VD_Calculate_current_u_D_born(int region, int d, DOF *u, DOF *u_D, DOF **p, FLOAT value, DOF *size_p, DOF *Epsilon);

void lambda_reg(DOF *v, SIMPLEX *ee, int bno, const FLOAT *lambda, FLOAT *value);

FLOAT Calculate_dof_circle_value(int region, DOF *dof, FLOAT r, FLOAT value);

FLOAT Calculate_dof_value(int region, int d, DOF *dof, FLOAT value);

FLOAT Calculate_dof_positive_value(int region, int d, DOF *dof, FLOAT value);

static void face_measure(BTYPE BDRY, GRID *g, INT *n, FLOAT *S);

static void region_measure(int region, GRID *g, INT *n, FLOAT *V);

static void bdry_quad(DOF *p, FLOAT *value);

static void Reactive_rate(DOF *p, DOF *u, FLOAT *value);

static void Reactive_P_rate(DOF *P, DOF *u, FLOAT *value);

static void Estimate_Poisson_error(DOF *u, DOF **p, DOF **indicator, DOF **error);

static void Estimate_N_P_error(INT TD, DOF *u, DOF *TMP_u, DOF **p, DOF **TMP_p, DOF **indicator, DOF **error);

static void linear_interpolation(COORD a, COORD b, COORD c, FLOAT value, int d);

static void triangle_area(FLOAT *X, FLOAT *Y, FLOAT *area);

static void tetrahetron_area(FLOAT *X, FLOAT *Y, FLOAT *area);

static void dof_save(char *fn, int *step_n, DOF *dof, ...);

static void dof_load(char *fn, int *step_n, DOF *dof, ...);

FLOAT P_B_Dof_Int(int region, DOF *U);

FLOAT Ber(FLOAT x); //rgshen2021-4-19

/*functions in PNP_coefficient.c*/

static void Init_coefficient();

static void change_file_path(char *path, char **file_name);

static void Read_parameter(INT read_parameter);

static void Read_pqr(const char *fn);

static void Read_bcmap(const char *fn);

static void Read_ions(const char *fn);

static void Read_cube(const char *fn);   //written by rgshen 2020-10-11

static int bc_map(int bctype);

static void initial_bulk(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void td_initial_bulk(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void func_G(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void func_grad_G(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void func_pb(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void func_ub(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void func_ub1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void func_D(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void func_scd(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void get_Kappa(SIMPLEX *e, FLOAT *epsilon, FLOAT *kappa);

static void get_Kappa1(SIMPLEX *e, FLOAT *Kappa);


/* PNP_analytic.c */  // by rgshen

static void analytic_func_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void analytic_func_n(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void analytic_func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void analytic_func_f_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void analytic_func_f_n(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

static void analytic_init();

static void analytic_print();

static void analytic_refresh();

static void analytic_finalize();

#endif
