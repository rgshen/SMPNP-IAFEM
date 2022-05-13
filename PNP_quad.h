
//element quad

//----without gradient
/* FLOAT PNP_Quad_1_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *u, int n, int order); */
#define	PNP_Quad_1_D_Bas(e, func, dof1, u, n, order) \
	PNP_Quad_5_D_Bas(e, func, dof1, NULL, NULL, NULL, NULL, u, n, order)
/* FLOAT PNP_Quad_2_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int n, int order); */
#define	PNP_Quad_2_D_Bas(e, func, dof1, dof2, u, n, order) \
	PNP_Quad_5_D_Bas(e, func, dof1, dof2, NULL, NULL, NULL, u, n, order)
/* FLOAT PNP_Quad_3_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int n, int order); */
#define	PNP_Quad_3_D_Bas(e, func, dof1, dof2, dof3, u, n, order) \
	PNP_Quad_5_D_Bas(e, func, dof1, dof2, dof3, NULL, NULL, u, n, order)
/* FLOAT PNP_Quad_4_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int n, int order); */
#define	PNP_Quad_4_D_Bas(e, func, dof1, dof2, dof3, dof4, u, n, order) \
	PNP_Quad_5_D_Bas(e, func, dof1, dof2, dof3, dof4, NULL, u, n, order)
FLOAT PNP_Quad_5_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int n, int order);

//----with gradient
/* FLOAT PNP_Quad_1_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *u, DOF *v, int n, int order); */
#define	PNP_Quad_1_D_G_D_G_B(e, func, dof1, u, v, n, order) \
	PNP_Quad_5_D_G_D_G_B(e, func, dof1, NULL, NULL, NULL, NULL, u, v, n, order)
/* FLOAT PNP_Quad_2_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, DOF *v, int n, int order); */
#define	PNP_Quad_2_D_G_D_G_B(e, func, dof1, dof2, u, v, n, order) \
	PNP_Quad_5_D_G_D_G_B(e, func, dof1, dof2, NULL, NULL, NULL, u, v, n, order)
/* FLOAT PNP_Quad_3_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, DOF *v, int n, int order); */
#define	PNP_Quad_3_D_G_D_G_B(e, func, dof1, dof2, dof3, u, v, n, order) \
	PNP_Quad_5_D_G_D_G_B(e, func, dof1, dof2, dof3, NULL, NULL, u, v, n, order)
/* FLOAT PNP_Quad_4_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, DOF *v, int n, int order); */
#define	PNP_Quad_4_D_G_D_G_B(e, func, dof1, dof2, dof3, dof4, u, v, n, order) \
	PNP_Quad_5_D_G_D_G_B(e, func, dof1, dof2, dof3, dof4, NULL, u, v, n, order)
FLOAT PNP_Quad_5_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, DOF *v, int n, int order);

//2 bas

//----without gradient bas
/* FLOAT PNP_Quad_1_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_1_D_Bas_Bas(e, func, dof1, u, m, v, n, order) \
	PNP_Quad_5_D_Bas_Bas(e, func, dof1, NULL, NULL, NULL, NULL, u, m, v, n, order)
/* FLOAT PNP_Quad_2_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_2_D_Bas_Bas(e, func, dof1, dof2, u, m, v, n, order) \
	PNP_Quad_5_D_Bas_Bas(e, func, dof1, dof2, NULL, NULL, NULL, u, m, v, n, order)
/* FLOAT PNP_Quad_3_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_3_D_Bas_Bas(e, func, dof1, dof2, dof3, u, m, v, n, order) \
	PNP_Quad_5_D_Bas_Bas(e, func, dof1, dof2, dof3, NULL, NULL, u, m, v, n, order)
/* FLOAT PNP_Quad_4_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_4_D_Bas_Bas(e, func, dof1, dof2, dof3, dof4, u, m, v, n, order) \
	PNP_Quad_5_D_Bas_Bas(e, func, dof1, dof2, dof3, dof4, NULL, u, m, v, n, order)
FLOAT PNP_Quad_5_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int m, DOF *v, int n, int order);

//----with 2 gradient bas
/* FLOAT PNP_Quad_1_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_1_D_G_B_G_B(e, func, dof1, u, m, v, n, order) \
	PNP_Quad_5_D_G_B_G_B(e, func, dof1, NULL, NULL, NULL, NULL, u, m, v, n, order)
/* FLOAT PNP_Quad_2_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_2_D_G_B_G_B(e, func, dof1, dof2, u, m, v, n, order) \
	PNP_Quad_5_D_G_B_G_B(e, func, dof1, dof2, NULL, NULL, NULL, u, m, v, n, order)
/* FLOAT PNP_Quad_3_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_3_D_G_B_G_B(e, func, dof1, dof2, dof3, u, m, v, n, order) \
	PNP_Quad_5_D_G_B_G_B(e, func, dof1, dof2, dof3, NULL, NULL, u, m, v, n, order)
/* FLOAT PNP_Quad_4_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int m, DOF *v, int n, int order); */
#define	PNP_Quad_4_D_G_B_G_B(e, func, dof1, dof2, dof3, dof4, u, m, v, n, order) \
	PNP_Quad_5_D_G_B_G_B(e, func, dof1, dof2, dof3, dof4, NULL, u, m, v, n, order)
FLOAT PNP_Quad_5_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int m, DOF *v, int n, int order);

//----with 1 gradient bas
/* FLOAT PNP_Quad_1_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *u, int m, DOF *gdof, DOF *v, int n, int order); */
#define	PNP_Quad_1_D_Bas_G_D_G_B(e, func, dof1, u, m, gdof, v, n, order) \
	PNP_Quad_4_D_Bas_G_D_G_B(e, func, dof1, NULL, NULL, NULL, u, m, gdof, v, n, order)
/* FLOAT PNP_Quad_2_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int m, DOF *gdof, DOF *v, int n, int order); */
#define	PNP_Quad_2_D_Bas_G_D_G_B(e, func, dof1, dof2, u, m, gdof, v, n, order) \
	PNP_Quad_4_D_Bas_G_D_G_B(e, func, dof1, dof2, NULL, NULL, u, m, gdof, v, n, order)
/* FLOAT PNP_Quad_3_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int m, DOF *gdof, DOF *v, int n, int order); */
#define	PNP_Quad_3_D_Bas_G_D_G_B(e, func, dof1, dof2, dof3, u, m, gdof, v, n, order) \
	PNP_Quad_4_D_Bas_G_D_G_B(e, func, dof1, dof2, dof3, NULL, u, m, gdof, v, n, order)
FLOAT PNP_Quad_4_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int m, DOF *gdof, DOF *v, int n, int order);
//dim of every dof is 1
//at least 1 dof
FLOAT PNP_Quad_5_D_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int n, int order) {
	int i1, i2 ,i3, i4, i5, j1;
	int i;
	FLOAT tmp;
	FLOAT d;
	const FLOAT *v1, *v2, *v3, *v4, *v5, *x1, *w;
	QUAD *quad;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
					if(dof5 != NULL) {
						i5 = DofTypeOrder(dof5, e);
					}
					else {
						i5 = 0;
					}
				}
				else {
					i4 = i5 = 0;
				}
			}
			else {
				i3 = i4 = i5 = 0;
			}
		}
		else {
			i2 = i3 = i4 = i5 = 0;
		}
		j1 = DofTypeOrder(u, e);
		order = i1 + i2 + i3 + i4 + i5 + j1;
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	if(dof2 != NULL) {
		v2 = phgQuadGetDofValues(e, dof2, quad);
		if(dof3 != NULL) {
			v3 = phgQuadGetDofValues(e, dof3, quad);
			if(dof4 != NULL) {
				v4 = phgQuadGetDofValues(e, dof4, quad);
				if(dof5 != NULL) {
					v5 = phgQuadGetDofValues(e, dof5, quad);
				}
			}
		}
	}
	x1 = phgQuadGetBasisValues(e, u, n, quad);

	if(func == NULL)
		func = func_id;

	d = 0.;
	if(dof2 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(x1++) * *(w++);
		}
	}
	else if(dof3 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(x1++) * *(w++);
		}
	}
	else if(dof4 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(x1++) * *(w++);
		}
	}
	else if(dof5 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(x1++) * *(w++);
		}
	}
	else {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(v5++) * *(x1++) * *(w++);
		}
	}

	return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_5_D_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, DOF *v, int n, int order) {
	int i1, i2 ,i3, i4, i5, j1 ,j2;
	int i, j;
	FLOAT tmp;
	FLOAT d, d0;
	const FLOAT *v1, *v2, *v3, *v4, *v5, *x1, *x2, *w;
	QUAD *quad;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
					if(dof5 != NULL) {
						i5 = DofTypeOrder(dof5, e);
					}
					else {
						i5 = 0;
					}
				}
				else {
					i4 = i5 = 0;
				}
			}
			else {
				i3 = i4 = i5 = 0;
			}
		}
		else {
			i2 = i3 = i4 = i5 = 0;
		}
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(v, e);
		order = i1 + i2 + i3 + i4 + i5 + j1 + j2;
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	if(dof2 != NULL) {
		v2 = phgQuadGetDofValues(e, dof2, quad);
		if(dof3 != NULL) {
			v3 = phgQuadGetDofValues(e, dof3, quad);
			if(dof4 != NULL) {
				v4 = phgQuadGetDofValues(e, dof4, quad);
				if(dof5 != NULL) {
					v5 = phgQuadGetDofValues(e, dof5, quad);
				}
			}
		}
	}
	x1 = phgQuadGetDofValues(e, u, quad);
	x2 = phgQuadGetBasisGradient(e, v, n, quad);

	if(func == NULL)
		func = func_id;

	d = 0.;
	if(dof2 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * d0 * *(w++);
		}
	}
	else if(dof3 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * d0 * *(w++);
		}
	}
	else if(dof4 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * d0 * *(w++);
		}
	}
	else if(dof5 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * d0 * *(w++);
		}
	}
	else {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(v5++) * d0 * *(w++);
		}
	}

	return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_5_D_Bas_Bas(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int m, DOF *v, int n, int order) {
	int i1, i2 ,i3, i4, i5, j1, j2;
	int i;
	FLOAT tmp;
	FLOAT d;
	const FLOAT *v1, *v2, *v3, *v4, *v5, *x1, *x2, *w;
	QUAD *quad;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
					if(dof5 != NULL) {
						i5 = DofTypeOrder(dof5, e);
					}
					else {
						i5 = 0;
					}
				}
				else {
					i4 = i5 = 0;
				}
			}
			else {
				i3 = i4 = i5 = 0;
			}
		}
		else {
			i2 = i3 = i4 = i5 = 0;
		}
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(v, e);
		order = i1 + i2 + i3 + i4 + i5 + j1 + j2;
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	if(dof2 != NULL) {
		v2 = phgQuadGetDofValues(e, dof2, quad);
		if(dof3 != NULL) {
			v3 = phgQuadGetDofValues(e, dof3, quad);
			if(dof4 != NULL) {
				v4 = phgQuadGetDofValues(e, dof4, quad);
				if(dof5 != NULL) {
					v5 = phgQuadGetDofValues(e, dof5, quad);
				}
			}
		}
	}
	x1 = phgQuadGetBasisValues(e, u, m, quad);
	x2 = phgQuadGetBasisValues(e, v, n, quad);
	
	if(func == NULL)
		func = func_id;

	d = 0.;
	if(dof2 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(x1++) * *(x2++) * *(w++);
		}
	}
	else if(dof3 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(x1++) * *(x2++) * *(w++);
		}
	}
	else if(dof4 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(x1++) * *(x2++) * *(w++);
		}
	}
	else if(dof5 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(x1++) * *(x2++) * *(w++);
		}
	}
	else {
		for(i = 0; i < quad->npoints; i++) {
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(v5++) * *(x1++) * *(x2++) * *(w++);
		}
	}

	return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_5_D_G_B_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int m, DOF *v, int n, int order) {
	int i1, i2 ,i3, i4, i5, j1 ,j2;
	int i, j;
	FLOAT tmp;
	FLOAT d, d0;
	const FLOAT *v1, *v2, *v3, *v4, *v5, *x1, *x2, *w;
	QUAD *quad;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
					if(dof5 != NULL) {
						i5 = DofTypeOrder(dof5, e);
					}
					else {
						i5 = 0;
					}
				}
				else {
					i4 = i5 = 0;
				}
			}
			else {
				i3 = i4 = i5 = 0;
			}
		}
		else {
			i2 = i3 = i4 = i5 = 0;
		}
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(v, e);
		order = i1 + i2 + i3 + i4 + i5 + j1 + j2;
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	if(dof2 != NULL) {
		v2 = phgQuadGetDofValues(e, dof2, quad);
		if(dof3 != NULL) {
			v3 = phgQuadGetDofValues(e, dof3, quad);
			if(dof4 != NULL) {
				v4 = phgQuadGetDofValues(e, dof4, quad);
				if(dof5 != NULL) {
					v5 = phgQuadGetDofValues(e, dof5, quad);
				}
			}
		}
	}
	x1 = phgQuadGetBasisGradient(e, u, m, quad);
	x2 = phgQuadGetBasisGradient(e, v, n, quad);

	if(func == NULL)
		func = func_id;

	d = 0.;
	if(dof2 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * d0 * *(w++);
		}
	}
	else if(dof3 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * d0 * *(w++);
		}
	}
	else if(dof4 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * d0 * *(w++);
		}
	}
	else if(dof5 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * d0 * *(w++);
		}
	}
	else {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x1++) * *(x2++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += tmp * *(v2++) * *(v3++) * *(v4++) * *(v5++) * d0 * *(w++);
		}
	}

	return d * phgGeomGetVolume(u->g, e);
}

FLOAT PNP_Quad_4_D_Bas_G_D_G_B(ELEMENT *e, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int m, DOF *gdof, DOF *v, int n, int order) {
	int i1, i2 ,i3, i4, j1 ,j2, j3;
	int i, j;
	FLOAT tmp;
	FLOAT d, d0;
	const FLOAT *v1, *v2, *v3, *v4, *x1, *x2, *x3, *w;
	QUAD *quad;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
				}
				else {
					i4 = 0;
				}
			}
			else {
				i3 = i4 = 0;
			}
		}
		else {
			i2 = i3 = i4 = 0;
		}
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(gdof, e);
		j3 = DofTypeOrder(v, e);
		order = i1 + i2 + i3 + i4 + j1 + j2 + j3;
	}
	quad = phgQuadGetQuad3D(order);

	w = quad->weights;

	v1 = phgQuadGetDofValues(e, dof1, quad);
	if(dof2 != NULL) {
		v2 = phgQuadGetDofValues(e, dof2, quad);
		if(dof3 != NULL) {
			v3 = phgQuadGetDofValues(e, dof3, quad);
			if(dof4 != NULL) {
				v4 = phgQuadGetDofValues(e, dof4, quad);
			}
		}
	}
	x1 = phgQuadGetBasisValues(e, u, m, quad);
	x2 = phgQuadGetDofValues(e, gdof, quad);
	x3 = phgQuadGetBasisGradient(e, v, n, quad);

	if(func == NULL)
		func = func_id;

	d = 0.;
	if(dof2 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x2++) * *(x3++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += *(x1++) * tmp * d0 * *(w++);
		}
	}
	else if(dof3 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x2++) * *(x3++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += *(x1++) * tmp * *(v2++) * d0 * *(w++);
		}
	}
	else if(dof4 == NULL) {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x2++) * *(x3++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += *(x1++) * tmp * *(v2++) * *(v3++) * d0 * *(w++);
		}
	}
	else {
		for(i = 0; i < quad->npoints; i++) {
			d0 = 0.;
			for(j = 0; j < Dim; j++) {
				d0 += *(x2++) * *(x3++);
			}
			tmp = *(v1++);
			func(&tmp);
			d += *(x1++) * tmp * *(v2++) * *(v3++) * *(v4++) * d0 * *(w++);
		}
	}

	return d * phgGeomGetVolume(u->g, e);
}

//face quad

/* FLOAT PNP_Quad_Face_1_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *u, int n, DOF *v, int m, int order); */
#define PNP_Quad_Face_1_D_Bas_Bas(e, face, func, dof1, u, n, v, m, order) \
	PNP_Quad_Face_5_D_Bas_Bas(e, face, func, dof1, NULL, NULL, NULL, NULL, u, n, v, m, order) 
/* FLOAT PNP_Quad_Face_2_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int n, DOF *v, int m, int order); */
#define PNP_Quad_Face_2_D_Bas_Bas(e, face, func, dof1, dof2, u, n, v, m, order) \
	PNP_Quad_Face_5_D_Bas_Bas(e, face, func, dof1, dof2, NULL, NULL, NULL, u, n, v, m, order) 
/* FLOAT PNP_Quad_Face_3_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int n, DOF *v, int m, int order); */
#define PNP_Quad_Face_3_D_Bas_Bas(e, face, func, dof1, dof2, dof3, u, n, v, m, order) \
	PNP_Quad_Face_5_D_Bas_Bas(e, face, func, dof1, dof2, dof3, NULL, NULL, u, n, v, m, order) 
/* FLOAT PNP_Quad_Face_4_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int n, DOF *v, int m, int order); */
#define PNP_Quad_Face_4_D_Bas_Bas(e, face, func, dof1, dof2, dof3, dof4, u, n, v, m, order) \
	PNP_Quad_Face_5_D_Bas_Bas(e, face, func, dof1, dof2, dof3, dof4, NULL, u, n, v, m, order) 
FLOAT PNP_Quad_Face_5_D_Bas_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int n, DOF *v, int m, int order) {
	int i, i1, i2, i3, i4, i5, j1, j2, v0, v1, v2;
	const FLOAT *x1, *x2, *x3, *x4, *x5;
	FLOAT tmp, d, d0, lambda[Dim + 1];
	FLOAT *buffer;
	const FLOAT *bas0, *bas1, *p, *w;
	QUAD *quad;
	DOF_TYPE *type_u, *type_v;

	assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
	assert(face >= 0 && face < NFace);

	type_u = u->type;
	type_v = v->type;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		x1 = phgAlloc(sizeof(FLOAT));
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			x2 = phgAlloc(sizeof(FLOAT));
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				x3 = phgAlloc(sizeof(FLOAT));
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
					x4 = phgAlloc(sizeof(FLOAT));
					if(dof5 != NULL) {
						i5 = DofTypeOrder(dof5, e);
						x5 = phgAlloc(sizeof(FLOAT));
					}
					else {
						i5 = 0;
					}
				}
				else {
					i4 = i5 = 0;
				}
			}
			else {
				i3 = i4 = i5 = 0;
			}
		}
		else {
			i2 = i3 = i4 = i5 = 0;
		}
		j1 = DofTypeOrder(u, e);
		j2 = DofTypeOrder(v, e);
		order = i1 + i2 + i3 + i4 + i5 + j1 + j2;
	}
	quad = phgQuadGetQuad2D(order);

	v0 = GetFaceVertex(face, 0);
	v1 = GetFaceVertex(face, 1);
	v2 = GetFaceVertex(face, 2);
	lambda[face] = 0.;

	if (n != m || u != v)
		buffer = phgAlloc(sizeof(*buffer));
	else
		buffer = NULL;
	p = quad->points;
	w = quad->weights;

	if(func == NULL)
		func = func_id;


	d = 0.;
	for (i = 0; i < quad->npoints; i++) {
	        lambda[v0] = *(p++);
	        lambda[v1] = *(p++);
	        lambda[v2] = *(p++);
		/* Note: type->BasFuncs returns an internal static buffer, its contents
		 * are changed by the next call if type_u==type_v */
		phgDofEval(dof1, e, lambda, x1);
		if(dof2 != NULL) {
			phgDofEval(dof2, e, lambda, x2);
			if(dof3 != NULL) {
				phgDofEval(dof3, e, lambda, x3);
				if(dof4 != NULL) {
					phgDofEval(dof4, e, lambda, x4);
					if(dof5 != NULL) {
						phgDofEval(dof5, e, lambda, x5);
					}
				}
			}
		}
		bas0 = type_u->BasFuncs(u, e, n, n + 1, lambda);
		if (n == m && u == v) {
			bas1 = bas0;
		}
		else {
			memcpy(buffer, bas0, sizeof(*buffer));
			bas0 = buffer;
			bas1 = type_v->BasFuncs(v, e, m, m + 1, lambda);
		}
		if(dof2 == NULL) {
			tmp = *x1;
			func(&tmp);
			d += tmp * *(bas0++) * *(bas1++) * *(w++);
		}
		else if(dof3 == NULL) {
			tmp = *x1;
			func(&tmp);
			d += tmp * *x2 * *(bas0++) * *(bas1++) * *(w++);
		}
		else if(dof4 == NULL) {
			tmp = *x1;
			func(&tmp);
			d += tmp * *x2 * *x3 * *(bas0++) * *(bas1++) * *(w++);
		}
		else if(dof5 == NULL) {
			tmp = *x1;
			func(&tmp);
			d += tmp * *x2 * *x3 * *x4 * *(bas0++) * *(bas1++) * *(w++);
		}
		else {
			tmp = *x1;
			func(&tmp);
			d += tmp * *x2 * *x3 * *x4 * *x5 * *(bas0++) * *(bas1++) * *(w++);
		}
	}

	if (n != m || u != v)
		phgFree(buffer);

	return d * phgGeomGetFaceArea(u->g, e, face);
}

/* FLOAT PNP_Quad_Face_1_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *u, int n, int order); */
#define PNP_Quad_Face_1_D_Bas(e, face, func, dof1, u, n, order) \
	PNP_Quad_Face_5_D_Bas(e, face, func, dof1, NULL, NULL, NULL, NULL, u, n, order) 
/* FLOAT PNP_Quad_Face_2_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *u, int n, int order); */
#define PNP_Quad_Face_2_D_Bas(e, face, func, dof1, dof2, u, n, order) \
	PNP_Quad_Face_5_D_Bas(e, face, func, dof1, dof2, NULL, NULL, NULL, u, n, order) 
/* FLOAT PNP_Quad_Face_3_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *u, int n, int order); */
#define PNP_Quad_Face_3_D_Bas(e, face, func, dof1, dof2, dof3, u, n, order) \
	PNP_Quad_Face_5_D_Bas(e, face, func, dof1, dof2, dof3, NULL, NULL, u, n, order) 
/* FLOAT PNP_Quad_Face_4_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *u, int n, int order); */
#define PNP_Quad_Face_4_D_Bas(e, face, func, dof1, dof2, dof3, dof4, u, n, order) \
	PNP_Quad_Face_5_D_Bas(e, face, func, dof1, dof2, dof3, dof4, NULL, u, n, order) 
FLOAT PNP_Quad_Face_5_D_Bas(ELEMENT *e, int face, PNP_FUNC_1 func, DOF *dof1, DOF *dof2, DOF *dof3, DOF *dof4, DOF *dof5, DOF *u, int n, int order) {
	int i, i1, i2, i3, i4, i5, j, v0, v1, v2;
	const FLOAT *x1, *x2, *x3, *x4, *x5;
	FLOAT tmp, d, lambda[Dim + 1];
	const FLOAT *bas, *p, *w;
	QUAD *quad;
	DOF_TYPE *type_u;

	assert(face >= 0 && face < NFace);

	type_u = u->type;

	if (order < 0) {
		i1 = DofTypeOrder(dof1, e);
		x1 = phgAlloc(sizeof(FLOAT));
		if(dof2 != NULL) {
			i2 = DofTypeOrder(dof2, e);
			x2 = phgAlloc(sizeof(FLOAT));
			if(dof3 != NULL) {
				i3 = DofTypeOrder(dof3, e);
				x3 = phgAlloc(sizeof(FLOAT));
				if(dof4 != NULL) {
					i4 = DofTypeOrder(dof4, e);
					x4 = phgAlloc(sizeof(FLOAT));
					if(dof5 != NULL) {
						i5 = DofTypeOrder(dof5, e);
						x5 = phgAlloc(sizeof(FLOAT));
					}
					else {
						i5 = 0;
					}
				}
				else {
					i4 = i5 = 0;
				}
			}
			else {
				i3 = i4 = i5 = 0;
			}
		}
		else {
			i2 = i3 = i4 = i5 = 0;
		}
		j = DofTypeOrder(u, e);
		order = i1 + i2 + i3 + i4 + i5 + j;
	}
	quad = phgQuadGetQuad2D(order);

	v0 = GetFaceVertex(face, 0);
	v1 = GetFaceVertex(face, 1);
	v2 = GetFaceVertex(face, 2);
	lambda[face] = 0.;

	p = quad->points;
	w = quad->weights;

	if(func == NULL)
		func = func_id;


	d = 0.;
	for (i = 0; i < quad->npoints; i++) {
	        lambda[v0] = *(p++);
	        lambda[v1] = *(p++);
	        lambda[v2] = *(p++);
		/* Note: type->BasFuncs returns an internal static buffer, its contents
		 * are changed by the next call if type_u==type_v */
		phgDofEval(dof1, e, lambda, x1);
		if(dof2 != NULL) {
			phgDofEval(dof2, e, lambda, x2);
			if(dof3 != NULL) {
				phgDofEval(dof3, e, lambda, x3);
				if(dof4 != NULL) {
					phgDofEval(dof4, e, lambda, x4);
					if(dof5 != NULL) {
						phgDofEval(dof5, e, lambda, x5);
					}
				}
			}
		}
		bas = type_u->BasFuncs(u, e, n, n + 1, lambda);
		if(dof2 == NULL) {
			tmp = *x1;
			func(&tmp);
			d += tmp * *bas * *(w++);
		}
		else if(dof3 == NULL) {
			tmp = *x1;
			func(&tmp);
			d += tmp * *x2 * *bas * *(w++);
		}
		else if(dof4 == NULL) {
			tmp = *x1;
			func(&tmp);
			d += tmp * *x2 * *x3 * *bas * *(w++);
		}
		else if(dof5 == NULL) {
			tmp = *x1;
			func(&tmp);
			d += tmp * *x2 * *x3 * *x4 * *bas * *(w++);
		}
		else {
			tmp = *x1;
			func(&tmp);
			d += tmp * *x2 * *x3 * *x4 * *x5 * *bas * *(w++);
		}
	}

	return d * phgGeomGetFaceArea(u->g, e, face);
}

