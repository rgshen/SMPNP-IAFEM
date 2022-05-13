#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"phg.h"

#include"PNP.h"



FLOAT K1, K2; //rgshen  

FLOAT asize_ratio = 6.022140857e-4;   //!SMPNP --Slotboom 变量 无量纲化（转成 "埃米" 尺度 而引入的比例系数）2021-5-28

#define CGSE 			0
#define SI 			1

INT SOU = SI;

#define Primitive		2

INT Method = Primitive;

//-----------rgshen20201201----------------------------
/* different transformation */
//#define Primitive       0 //c = p ---已经用liuxj 单独原始方法
#define Slotboom        1 //c = p exp(-z u)
#define Logarithm1      2 //p = p_bulk exp(c)
#define Logarithm2      3 //p = p_bulk exp(c - z u)

INT Transform_type = Primitive;
int tmp_type;           //for temporary use
//--------------------------------------------------

INT MOD = 0; //rgshen

#define SUPG 			1
#define PRFB 			2
#define EAFE         3

INT Stabilize = 0;

#define biomolecular		0
#define ion_channel		1
#define nanotube		2

INT model_type = biomolecular;

#define x_axis			0
#define y_axis			1
#define z_axis			2

INT current_direction = z_axis;

const int MNATOMS = 200000; /* The upper bound of atoms numbers */
int NATOMS = 0;
const int MNBCMAP = 10;/* The upper bound of bc_mark numbers*/
int NBCMAP = 0;
//const int MNION = 10;/* The upper bound of ions numbers */
#define MNION 10
int NION = 3;

int NPRO = 0;//number of product species --rgshen <--xujj


INT Example = -1;
INT energy_calculation = 0;
INT face_calculation = 0;
INT region_calculation = 0;
INT current_calculation = 0;
INT chemical_potential_calculation = 0;
INT potential_energy_calculation = 0;
INT dof_section_value_calculation = 0;
INT dof_circle_value_calculation = 0;
INT current_switch = 0; // switch the order of quad function
INT test_current = 0;   //test the order
INT voltage_test = 0;
INT born_test = 0;
INT born_correct = 0;

FLOAT section_value = 0.00;
FLOAT delta_value = 1.00;

#define r2(x,y,z,a,b,c) (Sqrt((x-a)*(x-a)+(y-b)*(y-b)+(z-c)*(z-c)))

/* Initialize constants */

/* the unit of coordinate in mesh is e-10 */

//char *fn_ions = "../PNP_mesh/1MAG.ions";
//char *fn_pqr = "../PNP_mesh/connexin.pqr";
//char *fn_mesh = "../PNP_mesh/connexin.mesh";
//char *fn_pqr = "../PNP_mesh/1MAG.pqr";
//char *fn_mesh = "../PNP_mesh/1MAG_2.mesh";
//char *fn_bcmap = "../PNP_mesh/1MAG.bc_dirichlet";
//char *fn_bcmap = "../PNP_program/PNP_mesh/1MAG.bc_neumann";

//char *fn_ions = "3LDC_channel.ions";
//char *fn_pqr = "3LDC_1.pqr";
//char *fn_mesh = "../PNP_mesh/3LDC_channel2.mesh";
//char *fn_mesh = "../PNP_mesh/1bl8_channel_nano.mesh";
//char *fn_mesh = "../PNP_mesh/1bl8_channel_tms.mesh";

char *fn_mesh = "../PNP_mesh/1bl8_tu.mesh";
char *fn_pqr = "../PNP_mesh/1bl8.pqr";  
//char *fn_ions = "1bl8_3ions.ions";
char *fn_bcmap = "1bl8_channel.bc_dirichlet";

//char *fn_mesh = "3LDC_protein.mesh";
//char *fn_mesh = "k_channel_2.mesh";
//char *fn_mesh = "3LDC_protein+k_1.mesh";
//char *fn_mesh = "3LDC_protein+k_2.mesh";
//char *fn_mesh = "3LDC_protein+k_3.mesh";
//char *fn_bcmap = "3LDC_channel.bc_dirichlet";
//char *fn_bcmap = "3LDC_protein.bc_dirichlet";
//char *fn_mesh = "../PNP_mesh/2JK4.mesh";
//char *fn_pqr = "../PNP_mesh/2JK4.pqr";
//char *fn_ions = "1bl8_2ions.ions";

char *fn_ions = "cylinder_nanopore.ions";
//char *fn_mesh = "cylinder_nanopore_R4.mesh";
//char *fn_mesh = "../PNP_mesh/cylinder_nanopore_R2.mesh";  // 3 = 3 + 9
//char *fn_mesh = "../PNP_mesh/new_cylinder1_refine_new_bd.mesh";// z = -z
//char *fn_bcmap = "cylinder_nanopore.bc_dirichlet";

//char *fn_ions = "circular_cone.ions";
//char *fn_mesh = "circular_cone_2_5.mesh";
//char *fn_bcmap = "circular_cone.bc_dirichlet";

//char *fn_ions = "cylinder_nanopore.ions";
//char *fn_mesh = "../PNP_mesh/cylinder_cone_2_30.mesh";
//char *fn_bcmap = "cylinder_cone_2_30.bc_dirichlet";


//---------瑞刚加的 -
//char *fn_mesh = "./mesh/cube4.dat";   // by rgshen for benchmark model to test convergence order
                                      //  立方体盒子 ，一致加密工作
/* * * * * * * * * * * */
/* for input and output*/
/* * * * * * * * * * * */

//INT read_parameter = 0;
char *fn_parameter = NULL;  //control parameter
char *fn_input_path = NULL; //used in transformation of input path
char *fn_output_path = NULL;    //used in transformation of output path
//char *tmp_file = NULL;      //used in combinition of path string
INT fn_dof_type = 1;        //0-save dof in binary file
                //1-save dof in normal file

//output
char *fn_dof = NULL;        //save SS dof
char *fn_td_dof = NULL;     //save TD dof
char *fn_current = NULL;    //save current
char *dof_name;
char *td_dof_name;
char *current_name;
//INT off_vtk = 0;        //turn off steady-state VTK print
INT t_d_off_vtk = 1;        //turn off time-dependent(TD) VTK print
INT ss_step_vtk = 0;        //print VTK every steady-state step
INT t_d_pri_int = 100;      //interval of VTK print in TD method
INT mono_pri_int = 0;       //interval of VTK print in monolithic method
INT reac_pri_int = -1;      //interval of reaction rate

/* * * * * * */
/* end input */
/* * * * * * */

/* * * * * * * * * */
/* for loop maxnum */   //<--xujj
/* * * * * * * * * */

INT mono_maxit = 20;
INT prim_maxit = 500;
INT slot_maxit = 50;
INT loga_maxit = 20;
INT time_maxit = 1000;
INT nonlinear_maxit = 50;

/* * * * * * * * * */
/* end loop maxnum */
/* * * * * * * * * */

/* * * * * * */

INT read_parameter = 0;
char *tmp_file = NULL;
char *tmp_path;
char read_tmp[5][30];

char *dof_file = NULL;
char *out_dof_section_file = NULL;
char *out_dof_circle_file = NULL;
char *out_current_file = NULL;

FLOAT CGSE_unit = 1e8;/* 1cm = 1e8A */
FLOAT SI_unit = 1e10;/* 1m = 1e10A */
FLOAT change_mesh_unit;
FLOAT mesh_unit = 1.0;/* default 1 A */
FLOAT change_unit;

int which_ion = 0;
int which_bc;

//INT solute_region = 1;
//INT solvent_region = 2;
INT solute_region = 120;
INT solvent_region = 100;

INT membrane_region = 110;

//----rgshen added20201201
 //int bc_mark[MNBCMAP] = {0}; //tab used in *.mesh
 //int u_bc_tab[MNBCMAP] = {0};
//int p_bc_tab[MNION][MNBCMAP] = {{0}};
//---------------------------------------------------


int bc_mark[10] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
int bc_tab[10] = {0};
BTYPE bmlc_bdry = 0;//biomolecule boundary
BTYPE reac_bdry = 0;//reactive boundary
BTYPE diri_bdry = 0;//all Dirichlet boundary
BTYPE scd_bdry = 0;//surface_charge boundary
BTYPE u_bdry = 0;//only u boundary
BTYPE c_bdry = 0;//channel boundary: H = G on c_bdry and Neumann for u

BTYPE sfcg_bdry = 0; //surface charge boundary--rgshen <--xujj


int count = 0;//counting the loop times in every method

//FLOAT biomolecular_size = 400;//ache
//FLOAT bulk_radius = 200;//ache
FLOAT biomolecular_size = 200;	//atom
FLOAT bulk_radius = 100;	//atom
FLOAT bulk_length = 10;

/*for nanopore*/
//FLOAT mesh_top = 50; 		//nanopore
//FLOAT mesh_bottom = 50; 	//nanopore
//FLOAT channel_top = 25; 	//nanopore
//FLOAT channel_bottom = 25; 	//nanopore
FLOAT scd_top = 25;
FLOAT scd_bottom = -25;
INT V_scd = 0;

/*for GA channel*/
//FLOAT mesh_top = 46.6; 		
//FLOAT mesh_bottom = 53.4; 	
//FLOAT channel_top = 7.0; 	
//FLOAT channel_bottom = 14.0; 	

/*for 2JK4*/
//FLOAT channel_top = 20.0; 	
//FLOAT channel_bottom = 10.0; 	

/*for connexin*/
//FLOAT mesh_top = 138; 		
//FLOAT mesh_bottom = 62; 	
//FLOAT channel_top = 55.0; 	
//FLOAT channel_bottom = -20.0; 	

/*for k+ channel box = 200*200*200 */
//FLOAT mesh_top = 101.217; 	
//FLOAT mesh_bottom = 98.783; 	
FLOAT mesh_top = 100.0; 	
FLOAT mesh_bottom = 100.0; 	
FLOAT channel_top = 15.0; 	
FLOAT channel_bottom = 15.0; //因为code里默认是负数，所以要输入正值 rgshen	

FLOAT buffer_length = 5;
FLOAT buffer_length_eps = 5;
FLOAT chan_ratio = 1.0 / 18.0;
INT V_D = 0;
INT V_N = 9;
FLOAT s_0 = 10.0;
//FLOAT a_1 = Epsilon_c / Epsilon_s;
//FLOAT z_0 = (channel_top + channel_top - buffer_length_eps) / 2.0;
FLOAT delta_z_E = 1.0;
FLOAT delta_z_D = 1.0;
INT Epsilon_type = 1;
INT V_D_type = 0;

INT all_level = 0;
INT region_level = -1, refine_region = 1;
INT adaptive_level = 1, adaptive_refine = -1, prior_refine = 0;
INT refine_time = 0;
INT input_data = 0;
INT loop_maxit = 500;
INT mem_pri_int = 0;

INT off_vtk = 0;/* set 1, if don't want vtk output with in steady-state method */

FLOAT u_bulk = 0.00;
FLOAT u_bulk1 = 0.00;
INT u_init = 1;

FLOAT init_bulk[MNION] = {0.0}; // rgshen <--xujj


FLOAT bulk_ratio = 0.0;/* ratio of density of reaction product to that of reactive ion*/
FLOAT diff_ratio = 1.0;/* ratio of diffusion of reaction product to that of reactive ion*/

FLOAT Lambda = 1;
FLOAT grow_Lambda = 1.0;
FLOAT initial_Alpha = 0.2;

FLOAT Ken = 10;


/* * * * * * * * * * * */ 
/* for reaction system */   //沈瑞刚加的
/* * * * * * * * * * * */

INT rr_bdry = 0;        //enable Robin reactive boundary
FLOAT rr_Alpha = 0.0;       //reactive rate Alpha
INT simple_reaction = 0;    //calculate on rr_bdry without gradient
INT rp_tab = 1;         //for reaction product TD calculation
                //0-substrate and product -> nonreactive ion
                //1-consider the product as steady state;
                //2-ordinary(problem)
INT lubrication = 0;        //when lubrication > 0
                //lub_ratio increase every step
                //when lubrication < 0
                //lub_ratio increase everytime its
FLOAT lub_ratio = 1.0;
FLOAT reac_bulk_ratio = 0.0;    //density ratio of product to substrate
FLOAT chan_bulk_ratio = 1.0;    //density ratio of outer
//FLOAT diff_ratio = 1.0;     //diffusion ratio of product to substrate

/* * * * * * * * * * * */
/* end reaction system */
/* * * * * * * * * * * */

FLOAT *tsc;         /* rgshen -original transformed results of c*/
//FLOAT num, delta;   /* rgshen -num is the maximum radius of sv func, delta is the unit step*/
FLOAT NUM = 70;
FLOAT delta = 0.05;


const FLOAT R = 1.0; /* inner radius of sphere model*/
const FLOAT NA = 6.022e+23; /* Avogadro constant */
const FLOAT T = 298; /* temperature */
const FLOAT Epsilon_m = 2.0; /* relative permittivity in \Omega_m */
const FLOAT Epsilon_s = 78.0; /* relative permittivity in \Omega_s */
FLOAT Epsilon_c = 78.0; /* relative permittivity in channel*/
FLOAT Epsilon_vacuum = 0.0; /* vacuum permittivity */
FLOAT Theta = 0.8; /* theta in element refinement */
FLOAT Alpha = 0.950;
FLOAT ec, Beta, Kappa_2, kB, kcal, Is, p_bulk;

FLOAT a0 = 3.1; /* size of H2O */
FLOAT a[MNION] = {0.0};
FLOAT K_i[MNION] = {0.0};
FLOAT bulk[MNION] = {0.0};
//FLOAT bulk_t[MNION] = {0.0};  //for bulk concerntration on the top 
//FLOAT bulk_b[MNION] = {0.0};  //for bulk concerntration on the bottom
FLOAT tol = 1.0e-6;

FLOAT A = 0.27;  // for variable dielectic
FLOAT Mu = 0.1; // for born term

FLOAT net_charge = 0;
FLOAT initial_tol = 0;
INT surface_charge = 0;
FLOAT scd = -0.025; //surface charge density
//FLOAT scd = 1.0 / 4.0 / M_PI;
FLOAT tmp_scd = 0; //store scd
INT nonlinear_Epsilon = 0;
INT solve_u_D = 0;
FLOAT PNP_grow = 0;
INT PNP_grow_n = 1;
FLOAT pqr_grow = 0;
INT pqr_grow_n = 1;

int exflag = 0;
typedef struct atoms_{
	FLOAT x, y, z;
	FLOAT Z;
	FLOAT r;
}atoms;

atoms *atom;

typedef struct ions_{
	FLOAT Z;//N
	FLOAT c;//mol/L
	//FLOAT c_t;//mol/L
	//FLOAT c_b;//mol/L
	FLOAT D;//A^2/ps
	FLOAT a;//A
	BOOLEAN reactive;
}ions;

ions ion[MNION]={{+1, 0.04, 0.196, 0, FALSE}, {-1, 0.05, 0.203, 0, FALSE}, {+1, 0.01, 0.203, 0, TRUE}};
//ions ion[MNION]={{+1, 0.04, 0.04, 0.196, 0, FALSE}, {-1, 0.05, 0.05, 0.203, 0, FALSE}, {+1, 0.01, 0.01, 0.203, 0, TRUE}};


static void Init_coefficient() {
//----------rgshen added for smpnp2021-4-28
  
   //num = 70.0;
   //delta = 0.05;

  
    //FLOAT *tsc; 
    int i;
    if(!(tsc = (FLOAT *)phgAlloc(NION * sizeof(FLOAT))))
        phgError(1, "Error in memory allocation for tsc.\n");
    FLOAT tm, tmp;
	tm = 0.0;
    for(i = 0; i < NION; i++){
        tm += Pow(a[i], 3) * ion[i].c;
    }
    tmp = 1 - tm * NA * 1e-27;
    if(tmp < 0.0) tmp = -tmp;
    for(i = 0; i < NION; i++){
        tsc[i] = ion[i].c * Pow(tmp, -K_i[i]);
    }
//--------------------------------------------


	if(current_calculation) {
		if(!SOU) {
			phgPrintf("Current must be calculated in SI, set SOU 1\n");
			SOU = 1;
		}
	}
	if(SOU == CGSE) {
		ec = 4.8032424e-10; /* single charge */
		kB = 1.38064852e-16; /* Boltzmann constant */
		Beta = 1 / (kB * T);
		Kappa_2 = 4 * M_PI * Beta * Pow(ec, 2);
		change_unit = CGSE_unit / mesh_unit;
		change_mesh_unit = 1.0e9 / mesh_unit;/* dm -> mesh_unit */
		kcal = 4.184e10;
	}
	if(SOU == SI) {
		ec = 1.6021892e-19; /* single charge */
		kB = 1.38064852e-23; /* Boltzmann constant */
		Epsilon_vacuum = 8.854187817e-12; /* vacuum permittivity */
		Beta = 1 / (kB * T);
		Kappa_2 = Beta * Pow(ec, 2) / Epsilon_vacuum;
		change_unit = SI_unit / mesh_unit;
		change_mesh_unit = 1.0e9 / mesh_unit;
		kcal = 4.184e3;
	}
	if(tmp_file != NULL) {
		//to combine path and file name
		tmp_path = (char *)malloc(64 * sizeof(char));
		int length = strlen(tmp_file);
		strcpy(tmp_path, tmp_file);
		if(tmp_path[length - 1]!='/') {
			tmp_path[length] = '/';
			tmp_path[length + 1] = '\0';
		}
		change_file_path(tmp_path, &fn_mesh);
		change_file_path(tmp_path, &fn_pqr);
		change_file_path(tmp_path, &fn_bcmap);
		change_file_path(tmp_path, &fn_ions);
	}
}

static void change_file_path(char *path, char **file_name) {
	if(tmp_file != NULL) {
		char *tmp = (char *)malloc(64 * sizeof(char));
		strcpy(tmp, path);
		strcat(tmp, *file_name);
		*file_name = tmp;
	}
}

static void Read_parameter(INT read_parameter) {
	if(read_parameter) {
		char *fn_para = "parameters.txt";
		FILE *f_in = fopen(fn_para, "rt");
		if(f_in == NULL) {
			phgError(1, "Can not open file %s.\n", fn_para);
		}
		free(dof_file);
		fscanf(f_in, "%*s%s%*s", read_tmp[0]);
		dof_file = (char *)malloc(64 * sizeof(char));
		strcpy(dof_file, read_tmp[0]);
		free(fn_mesh);
		fscanf(f_in, "%*s%s%*s", read_tmp[1]);
		fn_mesh = (char *)malloc(64 * sizeof(char));
		strcpy(fn_mesh, read_tmp[1]);
		free(fn_pqr);
		fscanf(f_in, "%*s%s%*s", read_tmp[2]);
		fn_pqr = (char *)malloc(64 * sizeof(char));
		strcpy(fn_pqr, read_tmp[2]);
		free(fn_bcmap);
		fscanf(f_in, "%*s%s%*s", read_tmp[3]);
		fn_bcmap = (char *)malloc(64 * sizeof(char));
		strcpy(fn_bcmap, read_tmp[3]);
		free(fn_ions);
		fscanf(f_in, "%*s%s%*s", read_tmp[4]);
		fn_ions = (char *)malloc(64 * sizeof(char));
		strcpy(fn_ions, read_tmp[4]);
		fscanf(f_in, "%*s%d%*s", &Method);
		fscanf(f_in, "%*s%lf%*s", &u_bulk);
		fscanf(f_in, "%*s%d%*s", &surface_charge);
		fscanf(f_in, "%*s%lf%*s", &scd);
		fclose(f_in);
		if(Method == 1 || Method == 3) Alpha = 1.0;
		if(Method == 2 || Method == 4) Alpha = 0.5;
		//int i = 0;
		//while(fscanf(f_in, "%lf", &ion[i++].c) != EOF);
	}
}

static void Read_pqr(const char *fn) {
	FILE *f_in;
	char pqr_atom[10];
	f_in = fopen(fn, "rt");
	if(f_in == NULL) {
		phgError(1, "Can not open file %s.\n", fn);
	}

	atom = (atoms *)malloc(MNATOMS * sizeof(*atom));
	NATOMS = 0;

	while(fscanf(f_in, "%s", pqr_atom) != EOF) {
		if(strcmp(pqr_atom, "ATOM") != 0) continue;
		fscanf(f_in, "%*s%*s%*s%*s%lf%lf%lf%lf%lf", &atom[NATOMS].x, &atom[NATOMS].y, &atom[NATOMS].z, &atom[NATOMS].Z, &atom[NATOMS].r);
		net_charge += atom[NATOMS].Z;
		NATOMS++;
	}
	fclose(f_in);
	phgPrintf("\nRead pqr file\n");
	phgPrintf("NET_CHARGE = %.4e e\n", net_charge);
	net_charge = 0;
}

static void Read_bcmap(const char *fn) {
	FILE *f_in;
	char bc_name[20];
	f_in = fopen(fn, "rt");
	if(f_in == NULL) {
		phgPrintf("No bc_map file\n");
//		phgError(1, "Can not open file %s.\n", fn);
	}
	else {
		phgPrintf("\nRead bc_map file\n");
		NBCMAP = 0;
		while(NBCMAP < MNBCMAP) {
			if(fscanf(f_in, "%s%d%d", bc_name, &bc_mark[NBCMAP], &bc_tab[NBCMAP]) != EOF) {
				if(bc_mark[NBCMAP] != -1) {
					phgPrintf("%s = %3.d	", bc_name, bc_mark[NBCMAP]);
					if(bc_tab[NBCMAP] == 0) {
						phgPrintf("Dirichlet	boundary");
						diri_bdry += BDRY_USER0 * pow(2, NBCMAP);
					}
					else if(bc_tab[NBCMAP] == 1) {
						phgPrintf("Dirichlet	boundary(u = 0)");
						diri_bdry += BDRY_USER0 * pow(2, NBCMAP);
					}
					else if(bc_tab[NBCMAP] == 2) {
						phgPrintf("Dirichlet	boundary(p = 0)");
						diri_bdry += BDRY_USER0 * pow(2, NBCMAP);
					}
					else if(bc_tab[NBCMAP] == 3) {
						phgPrintf("Biomolecular	boundary(nonreactive)");
						bmlc_bdry += BDRY_USER0 * pow(2, NBCMAP);
					}
					else if(bc_tab[NBCMAP] == 4) {
						phgPrintf("Biomolecular	boundary(reactive)");
						bmlc_bdry += BDRY_USER0 * pow(2, NBCMAP);
						reac_bdry += BDRY_USER0 * pow(2, NBCMAP);
					}
					else if(bc_tab[NBCMAP] == 5) {
						phgPrintf("Extraordinary	boundary(Neum for p, Diri for u != 0)");
						u_bdry += BDRY_USER0 * pow(2, NBCMAP);
					}
					else if(bc_tab[NBCMAP] == 6) {
						phgPrintf("Extraordinary	boundary(Neum for p, Diri for u == 0)");
						u_bdry += BDRY_USER0 * pow(2, NBCMAP);
					}
					else if(bc_tab[NBCMAP] == 7) {
						phgPrintf("Channel	boundary(Neum for u, Diri for H = G)");
						c_bdry += BDRY_USER0 * pow(2, NBCMAP);
					}
					else if(surface_charge && bc_tab[NBCMAP] == 8) {
						phgPrintf("Surface Charge	boundary");
						scd_bdry += BDRY_USER0 * pow(2, NBCMAP);
						bmlc_bdry += BDRY_USER0 * pow(2, NBCMAP);
					}
					else u_bdry += BDRY_USER0 * pow(2, NBCMAP);
					phgPrintf("\n");
					NBCMAP++;
				}
			}
			else
				break;
		}
		fclose(f_in);
		//phgPrintf("input %d kind of boundaries\n", NBCMAP);
	}
}

static void Read_ions(const char *fn) {
	FILE *f_in;
	f_in = fopen(fn, "rt");
	if(f_in == NULL) {
		phgPrintf("No ions file\n");
//		phgError(1, "Can not open file %s.\n", fn);
	}
	else {
		phgPrintf("\nRead ions file\n");
		NION = 0;
		while(NION < MNION) {
			if(fscanf(f_in, "%lf%lf%lf%lf%d", &ion[NION].Z, &ion[NION].c, &ion[NION].D, &ion[NION].a, &ion[NION].reactive) != EOF){
			//if(fscanf(f_in, "%lf%lf%lf%lf%lf%d", &ion[NION].Z, &ion[NION].c_t, &ion[NION].c_b, &ion[NION].D, &ion[NION].a, &ion[NION].reactive) != EOF){
				NION++;
			}
			else
				break;
		}
		fclose(f_in);
	}
}



// The folllowing written by rgshen for benchmark model to test convergence 2020-10-11   待完善
static void Read_cube(const char *fn) {
	FILE *f_in;
	f_in = fopen(fn, "rt");
	if(f_in == NULL) {
		phgPrintf("No ions file\n");
		phgError(1, "Can not open file %s.\n", fn);
	}
    else {
		phgPrintf("\nRead ions file\n");
	}
    fclose(f_in);
}




static int bc_map(int bctype) {
	if(bctype == bc_mark[0]) {
		if(bc_tab[0] >= 0 && bc_tab[0] <= 2)
			return (BDRY_USER0|DIRICHLET);
		else
			return BDRY_USER0;
	}
	else if(bctype == bc_mark[1]) {
		if(bc_tab[1] >= 0 && bc_tab[1] <= 2)
			return (BDRY_USER1|DIRICHLET);
		else
			return BDRY_USER1;
	}
	else if(bctype == bc_mark[2]) {
		if(bc_tab[2] >= 0 && bc_tab[2] <= 2)
			return (BDRY_USER2|DIRICHLET);
		else
			return BDRY_USER2;
	}
	else if(bctype == bc_mark[3]) {
		if(bc_tab[3] >= 0 && bc_tab[3] <= 2)
			return (BDRY_USER3|DIRICHLET);
		else
			return BDRY_USER3;
	}
	else if(bctype == bc_mark[4]) {
		if(bc_tab[4] >= 0 && bc_tab[4] <= 2)
			return (BDRY_USER4|DIRICHLET);
		else
			return BDRY_USER4;
	}
	else if(bctype == bc_mark[5]) {
		if(bc_tab[5] >= 0 && bc_tab[5] <= 2)
			return (BDRY_USER5|DIRICHLET);
		else
			return BDRY_USER5;
	}
	else if(bctype == bc_mark[6]) {
		if(bc_tab[6] >= 0 && bc_tab[6] <= 2)
			return (BDRY_USER6|DIRICHLET);
		else
			return BDRY_USER6;
	}
	else if(bctype == bc_mark[7]) {
		if(bc_tab[7] >= 0 && bc_tab[7] <= 2)
			return (BDRY_USER7|DIRICHLET);
		else
			return BDRY_USER7;
	}
	else if(bctype == bc_mark[8]) {
		if(bc_tab[8] >= 0 && bc_tab[8] <= 2)
			return (BDRY_USER8|DIRICHLET);
		else
			return BDRY_USER8;
	}
	else if(bctype == bc_mark[9]) {
		if(bc_tab[9] >= 0 && bc_tab[9] <= 2)
			return (BDRY_USER9|DIRICHLET);
		else
			return BDRY_USER9;
	}
	else return UNDEFINED;				  /**< belongs to a boundary but has no bdry type */

}

static void initial_bulk(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
	*value = bulk[which_ion];
}

static void td_initial_bulk(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
	if(model_type == biomolecular) {
		FLOAT r = r2(x,y,z,0,0,0);
		if(r < (biomolecular_size - bulk_radius)) *value = 0;
		else *value = bulk[which_ion] * (r - (biomolecular_size - bulk_radius)) / bulk_radius;
	}
	else if(model_type == ion_channel || model_type == nanotube) {
		*value = 0;
		if(mesh_top != 0 && mesh_bottom != 0)
			if(mesh_top + mesh_bottom < 2 * bulk_length)
				phgError(1, "\nbulk length is larger than half length of channel\n");
		if(mesh_top > 0 && z > (mesh_top - bulk_length)) {
			*value = bulk[which_ion] * (z - (mesh_top - bulk_length)) / bulk_length;
		}
		else if(mesh_bottom > 0 && z < - (mesh_bottom - bulk_length)) {
			*value = bulk[which_ion] * (- z - (mesh_bottom - bulk_length)) / bulk_length;
		}
	}
}

static void func_G(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
	FLOAT dist;
	int i;
	*value = 0;
	for(i = 0; i < NATOMS; i++) {
 		dist = r2(x, y, z, atom[i].x, atom[i].y, atom[i].z);
		if(dist < 0.05) continue;
 		*value += (atom[i].Z/dist);
	}
	//cm -> A: 1/r -> 1/x * 1.0e8 
//	*value *= (1.0e8 * Beta * ec * ec / Epsilon_m);
	*value *= change_unit * Kappa_2 / 4 / M_PI / Epsilon_m;
}

static void func_grad_G(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
	FLOAT ret[3], r;
	int i;
	ret[0] = ret[1] = ret[2] = 0.0;
	for (i = 0; i < NATOMS; i++) {
		r = r2(x, y, z, atom[i].x, atom[i].y, atom[i].z);
		if(r < 0.05) continue;
		ret[0] += -1.0 * atom[i].Z * (x - atom[i].x) / Pow(r,3);
		ret[1] += -1.0 * atom[i].Z * (y - atom[i].y) / Pow(r,3);
		ret[2] += -1.0 * atom[i].Z * (z - atom[i].z) / Pow(r,3);
	}
	value[0] = ret[0] * change_unit * Kappa_2 / 4 / M_PI / Epsilon_m;
	value[1] = ret[1] * change_unit * Kappa_2 / 4 / M_PI / Epsilon_m;
	value[2] = ret[2] * change_unit * Kappa_2 / 4 / M_PI / Epsilon_m;
}

static void func_pb(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
	*value = bulk[which_ion];
	if(PNP_grow < 0 && count < Fabs(PNP_grow)) {
		*value *= count / -(int)PNP_grow;
	}
	else if(PNP_grow > 0 && PNP_grow_n <= PNP_grow) {
		*value = *value * PNP_grow_n / (int)PNP_grow;
	}
}


//static void func_pb_t(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
//	*value = bulk_t[which_ion];
//}

//static void func_pb_b(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
//	*value = bulk_b[which_ion];
//}


static void func_ub(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
	if(u_bulk != 0) {
		*value = (Beta * ec) * (u_bulk);
	}
	else
		*value = 0;
	if(PNP_grow < 0 && count < Fabs(PNP_grow)) {
		*value *= count / -PNP_grow;
	}
	else if(PNP_grow > 0 && PNP_grow_n <= PNP_grow) {
		*value = *value * PNP_grow_n / (int)PNP_grow;
	}
}


static void func_ub1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
        if(u_bulk1 != 0) {
		*value = (Beta * ec) * (u_bulk1);
        }
        else
                *value = 0;
        if(PNP_grow < 0 && count < Fabs(PNP_grow)) {
                *value *= count / -PNP_grow;
        }
        else if(PNP_grow > 0 && PNP_grow_n <= PNP_grow) {
                *value = *value * PNP_grow_n / (int)PNP_grow;
        }
}


static void func_D(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
	FLOAT f_z;
	if(V_D == 1) {
		*value = ion[which_ion].D;
	}
	else if(V_D == 2) {
		if (V_D_type == 0) {
			if(z <= -(channel_bottom + buffer_length) || z >= (channel_top + buffer_length) ) {
				*value = ion[which_ion].D;
			}
			else if(z >= -channel_bottom && z <= channel_top) {
				*value = ion[which_ion].D * chan_ratio;
			}
			else if(z >= -(channel_bottom + buffer_length) && z <= - channel_bottom) {
				f_z = V_N * Pow((- z - channel_bottom) / buffer_length, V_N + 1)
			    	- (V_N + 1) * Pow((- z - channel_bottom) / buffer_length, V_N);
				//*value = ion[which_ion].D * (1.0 - (1.0 - chan_ratio) * f_z);
				*value = ion[which_ion].D * (chan_ratio - (1.0 - chan_ratio) * f_z);
			}
			else if(z <= (channel_top + buffer_length) && z >= channel_top) {
				f_z = V_N * Pow((z - channel_top) / buffer_length, V_N + 1)
			    	- (V_N + 1) * Pow((z - channel_top) / buffer_length, V_N);
				//*value = ion[which_ion].D * (1.0 - (1.0 - chan_ratio) * f_z);
				*value = ion[which_ion].D * (chan_ratio - (1.0 - chan_ratio) * f_z);
			}
		}
		if (V_D_type == 1) {
	                if(z <= -channel_bottom || z >= channel_top ) {
        	                *value = ion[which_ion].D;
                	}
                	else if(z >= -(channel_bottom - buffer_length) && z <= (channel_top - buffer_length)) {
        	                *value = ion[which_ion].D * chan_ratio;
              		}
	                else if(z >= - channel_bottom && z <= -(channel_bottom - buffer_length)) {
        	                *value = (ion[which_ion].D - ion[which_ion].D * chan_ratio) / (2*atan(s_0)) * atan(2 * s_0 / (-buffer_length) * z - s_0 * (-(channel_bottom - buffer_length) - channel_bottom) / (-buffer_length)) + (ion[which_ion].D + ion[which_ion].D * chan_ratio) / 2.0;
                	}
                	else if(z >= (channel_top - buffer_length) && z <= channel_top) {
                        	*value = (ion[which_ion].D - ion[which_ion].D * chan_ratio) / (2*atan(s_0)) * atan(2 * s_0 / buffer_length * z - s_0 * (channel_top - buffer_length + channel_top) / buffer_length) + (ion[which_ion].D + ion[which_ion].D * chan_ratio) / 2.0;
                	}
		}
		if (V_D_type == 2) {
		        FLOAT a_1 = chan_ratio;
		        FLOAT z_0 = (channel_top + channel_top - buffer_length) / 2.0;
	                *value = ion[which_ion].D * (a_1 + (1.0 - a_1) / (1.0 + exp(-(abs(z) - z_0) / delta_z_D)));
		}
	}
}

static void func_scd(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
	if(z >= scd_bottom && z <= scd_top) *value = scd;
	else *value = 0.0;
}

static void get_Kappa(SIMPLEX *e, FLOAT *Epsilon, FLOAT *Kappa) {
//	if(e->region_mark == solute_region) {
        if(e->region_mark == solute_region || e->region_mark == membrane_region){ 
		*Epsilon = Epsilon_m;
		*Kappa = 0;
	}
	else if(e->region_mark == solvent_region) {
		*Epsilon = Epsilon_s;
		*Kappa = Kappa_2;
	}
	else{
		phgError(1, "UNKNOWN REGION MARK!\n");
	}
} 
/*for variable dielectric*/
static void get_Kappa1(SIMPLEX *e, FLOAT *Kappa) {
//        if(e->region_mark == solute_region) { 
        if(e->region_mark == solute_region || e->region_mark == membrane_region){ 
                *Kappa = 0;
        }
        else if(e->region_mark == solvent_region) {
                *Kappa = Kappa_2;
        }
        else{
                phgError(1, "UNKNOWN REGION MARK!\n");
        }
}

static void initialize_options(INT mem_max) {
/* command line options preprocess */
	//about reading file
	phgOptionsRegisterFilename("tmp_file", "tmp_file", &tmp_file);
	phgOptionsRegisterFilename("fn_mesh", "mesh_file", &fn_mesh);
	phgOptionsRegisterFilename("fn_pqr", "fixed_ion_file", &fn_pqr);
	phgOptionsRegisterFilename("fn_bcmap", "boundary_mapping", &fn_bcmap);
	phgOptionsRegisterFilename("fn_ions", "different_ions", &fn_ions);
	phgOptionsRegisterFilename("dof_file", "save_dof", &dof_file);
        phgOptionsRegisterFilename("out_dof_section_file", "dof_section_value_out", &out_dof_section_file);
        phgOptionsRegisterFilename("out_dof_circle_file", "dof_circle_value_out", &out_dof_circle_file);
        phgOptionsRegisterFilename("out_current_file", "current_value_out", &out_current_file);

	//about many other coefficients
	phgOptionsRegisterFloat("tol", "Tolerance", &tol);
	phgOptionsRegisterFloat("Mu", "Mu", &Mu);
	phgOptionsRegisterFloat("initial_tol", "Monolithic_Initialization_Tolerance", &initial_tol);
	phgOptionsRegisterFloat("scd", "surface_charge_density", &scd);
	phgOptionsRegisterFloat("PNP_grow", "for_convergence", &PNP_grow);
	phgOptionsRegisterFloat("pqr_grow", "for_convergence", &pqr_grow);
	phgOptionsRegisterFloat("Epsilon_c", "the_dielectric_within_channel", &Epsilon_c);
	phgOptionsRegisterFloat("s_0", "the_smoother_of_Epsilon_type1", &s_0);
	phgOptionsRegisterFloat("delta_z_E", "the_smoother_of_Epsilon_type2", &delta_z_E);
	phgOptionsRegisterFloat("delta_z_D", "the_smoother_of_V_D_type2", &delta_z_D);
	phgOptionsRegisterFloat("Ken", "coefficient", &Ken);
	phgOptionsRegisterFloat("Alpha", "coefficient", &Alpha);
	phgOptionsRegisterFloat("Theta", "coefficient", &Theta);
	phgOptionsRegisterFloat("Lambda", "coefficient", &Lambda);
	phgOptionsRegisterFloat("initial_Alpha", "coefficient", &initial_Alpha);
	phgOptionsRegisterFloat("u_bulk", "u_boundary_bulk", &u_bulk);
	phgOptionsRegisterFloat("u_bulk1", "u_boundary1_bulk", &u_bulk1);
	phgOptionsRegisterFloat("section_value","value_of_section", &section_value);
	phgOptionsRegisterFloat("delta_value","delta_value_of_section", &delta_value);
	phgOptionsRegisterFloat("mesh_unit", "unit_of_mesh", &mesh_unit);
	phgOptionsRegisterFloat("biomolecular_size", "radius_of_biomolecular", &biomolecular_size);
	phgOptionsRegisterFloat("bulk_radius", "radius_of_bulk_distribution", &bulk_radius);
	phgOptionsRegisterFloat("mesh_top", "top_of_mesh", &mesh_top);
	phgOptionsRegisterFloat("mesh_bottom", "bottom_of_mesh", &mesh_bottom);
	phgOptionsRegisterFloat("channel_top", "top_of_channel", &channel_top);
	phgOptionsRegisterFloat("channel_bottom", "bottom_of_channel", &channel_bottom);
	phgOptionsRegisterFloat("buffer_length", "length_of_buffer_for_V_D", &buffer_length);
	phgOptionsRegisterFloat("buffer_length_eps", "length_of_buffer_for_Eps", &buffer_length_eps);
	phgOptionsRegisterFloat("chan_ratio", "ratio_of_diffusion_in_channel_to_outside", &chan_ratio);
	phgOptionsRegisterFloat("bulk_length", "length_of_bulk_distribution", &bulk_length);
	phgOptionsRegisterFloat("bulk_ratio", "ratio_of_bulk", &bulk_ratio);
	phgOptionsRegisterFloat("diff_ratio", "ratio_of_diffusion", &diff_ratio);

	phgOptionsRegisterInt("read_parameter", "enable_parameter", &read_parameter);
	phgOptionsRegisterInt("loop_maxit", "times_of_loop", &loop_maxit);
	phgOptionsRegisterInt("off_vtk", "turn_off_vtk_output", &off_vtk);
	phgOptionsRegisterInt("mem_pri_int", "memory_print_interval", &mem_pri_int);
	phgOptionsRegisterInt("V_scd", "variable_surface_charge_density", &V_scd);
	phgOptionsRegisterInt("V_D", "variable_diffusion_coefficient", &V_D);
	phgOptionsRegisterInt("V_N", "level_of_func_for_variable_diffusion", &V_N);
	phgOptionsRegisterInt("Example", "examples_to_test", &Example);
	phgOptionsRegisterInt("input_data", "inpute_data_for_3_ions", &input_data);
	phgOptionsRegisterInt("surface_charge", "charged_surface_boundary", &surface_charge);
	    phgOptionsRegisterInt("scd", "charged_surface_giving", &scd);  //rgshen 2021-1-10

	phgOptionsRegisterInt("all_level", "all_refine_level", &all_level);
	phgOptionsRegisterInt("solute_region", "region_of_solute", &solute_region);
	phgOptionsRegisterInt("solvent_region", "region_of_solvent", &solvent_region);
        phgOptionsRegisterInt("membrane_region", "region_of_membrane", &membrane_region);
	phgOptionsRegisterInt("region_level", "region_refine_level", &region_level);
	phgOptionsRegisterInt("refine_region", "region_to_refine", &refine_region);
	phgOptionsRegisterInt("adaptive_level", "adaptive_refine_level", &adaptive_level);
	phgOptionsRegisterInt("adaptive_refine", "refine_adaptively", &adaptive_refine);
	phgOptionsRegisterInt("prior_refine", "refine_before_converged", &prior_refine);
	phgOptionsRegisterInt("energy_calculation", "calculate_energy", &energy_calculation);
	phgOptionsRegisterInt("face_calculation", "calculate_face", &face_calculation);
	phgOptionsRegisterInt("region_calculation", "calculate_region", &region_calculation);
	phgOptionsRegisterInt("current_calculation", "calculate_current", &current_calculation);
	phgOptionsRegisterInt("chemical_potential_calculation", "calculate_chemical_potential", &chemical_potential_calculation);
	phgOptionsRegisterInt("potential_energy_calculation", "calculate_potential_energy", &potential_energy_calculation);
	phgOptionsRegisterInt("dof_section_value_calculation", "calculate_dof_section_value", &dof_section_value_calculation);
	phgOptionsRegisterInt("dof_circle_value_calculation", "calculate_dof_circle_value", &dof_circle_value_calculation);
	phgOptionsRegisterInt("current_direction", "direction_of_current", &current_direction);
	phgOptionsRegisterInt("model_type", "biomolecular_or_ionchannel", &model_type);
	phgOptionsRegisterInt("u_init", "initialization_of_u_bulk", &u_init);
	phgOptionsRegisterInt("SOU", "system_of_units", &SOU);
	phgOptionsRegisterInt("solve_u_D", "solve_u_D", &solve_u_D);
        phgOptionsRegisterInt("MOD", "different_mode", &MOD); //rgshen
		//phgOptionsRegisterInt("a0", "solvent_size", &a0); //rgshen2021-5-28
        phgOptionsRegisterInt("nonlinear_Epsilon", "different_epsilon", &nonlinear_Epsilon);
        phgOptionsRegisterInt("Epsilon_type", "different_epsilon_type", &Epsilon_type);
        phgOptionsRegisterInt("V_D_type", "different_diffusion_type", &V_D_type);
	phgOptionsRegisterInt("Method", "Method_to_solve", &Method);
	phgOptionsRegisterInt("Stabilize", "Stabilization_Scheme", &Stabilize);
	phgOptionsRegisterInt("mem_max", "Maximum memory (MB)", &mem_max);
	phgOptionsRegisterInt("current_switch", "advanced_current", &current_switch);
	phgOptionsRegisterInt("test_current", "test_current_order", &test_current);
	phgOptionsRegisterInt("voltage_test", "test_voltage", &voltage_test);
	phgOptionsRegisterInt("born_test", "test_born_energy", &born_test);
	phgOptionsRegisterInt("born_correct", "correct_born_energy", &born_correct);
	phgOptionsPreset("-mem_max 5120");
}
