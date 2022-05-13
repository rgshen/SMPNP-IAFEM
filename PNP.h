
/* significant definitions */

//avoiding mutiple definitions
#ifndef PNP_H
#define PNP_H



/* * * * * * * * */
/* for norms use */
/* * * * * * * * */

FLOAT u_L2, p_L2, P_L2, T_L2;

FLOAT u_H1, p_H1, P_H1, T_H1;

/* * * * * * * * */
/* end norms use */
/* * * * * * * * */

/* time dependency */
#define Steady_State        0 //consider a system in equilibrium
#define Time_Dependent      1 //system is changing over time
#define Test_Time           2 //test convergence

INT Time_type = Steady_State;

/* initialization in TD */
#define bulk_init       0 //with boundary bulk
#define func_init       1 //with concentration function
#define stst_init       2 //with steady-state result
#define load_init       3 //load data saved before
#define anal_init      -1 //with analytic function
#define mesh_init      -2 //with steady-state adaptive refinement

INT t_d_init_type = func_init;


/* * * * * * * */
/* for phg use */
/* * * * * * * */

//create solver for linear system
SOLVER *solver;

SOLVER *pc_solver;

SOLVER *pc_u, *pc_p;

//GRIDs
GRID *g;

GRID *fine_g;

GRID *coarse_g;

//DOFs
DOF *G, *H, *grad_G, *grad_H;   //initially calculated

DOF *scatter_p;         //consider fixed charge scattered as ions

//DOF *D[MNION];          //Diffusions differ in regions

DOF *epsl;          //permittivity influenced by ion species

DOF *bdry;          //boundary mark

DOF *ph;            //electrostatic potential (not dimensionless)

DOF *region_indicator;      //1 in region, 0 out of region

//dimensionless DOF
//basic     relaxation  last step   Newton's iteration
DOF *u,     *tmp_u,     *TMP_u,     *delta_u;

#define MNION1 10 // 为了解析解用--新定义的

DOF *c[MNION1],  *tmp_c[MNION1],  *TMP_c[MNION1],  *delta_c[MNION1];

//DOF *p[MNION],  *tmp_p[MNION],  *TMP_p[MNION],  *delta_p[MNION];

//DOF *P[MNION],  *tmp_P[MNION],  *TMP_P[MNION],  *delta_P[MNION];

//DOF *T[MNION],  *tmp_T[MNION],  *TMP_T[MNION],  *delta_T[MNION];

//DOF *size_p, *size_TMP_p;   //size effect

DOF *indicator, *error, *td_tmp; //for refinement use

DOF *lub_u;
/* * * * * * * */
/* end phg use */
/* * * * * * * */


/* * * * * * * * */
/* for count use */
/* * * * * * * * */
int lu_count = 0;       //count step for converge-lubrication
int ni_count = 0;       //count step for Newton's iteration
int ss_count = 0;       //count iteration step for steady-state
int mn_count = 0;       //count step for monolithic steady-state
int td_count = 0;       //count time step for time-dependent
int tl_count = 0;       //count step in every time level
int pr_count = 0;       //count VTK printing
int re_count = 0;       //count refine time
int break_flag = 0;     //break flag in counting

/* * * * * * * * */
/* end count use */
/* * * * * * * * */







/* * * * * * * * */
/* analytic test */
/* * * * * * * * */

INT analytic_test = 0;
DOF *analytic_u;
DOF *analytic_grad_u;
DOF *analytic_p;
DOF *analytic_n;
DOF *analytic_grad_p;
DOF *analytic_grad_n;
DOF *analytic_f_p;
DOF *analytic_f_n;
DOF *anal_err_u;
DOF *anal_err_p;
DOF *anal_err_n;


#endif	/* PNP_H */








