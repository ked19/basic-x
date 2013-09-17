#ifndef _GSL_MTXOP_H
#define _GSL_MTXOP_H

#include "denseMatrix.h"
#include "matrixOperation.h"
#include <gsl/gsl_linalg.h>

void CopyTo(const gsl_matrix &gsl, Mtx &mtx);
void CopyTo(const Mtx &mtx, gsl_matrix &gsl);
void CopyTo(const gsl_vector *gsl, Mtx &mtx);
void CopyTo(const Mtx &mtx, gsl_vector &gsl);

//*************************************************************************************************

class Inv
{
public:
	Inv();
	~Inv();

	void Gen(Mtx &mtxInv, Mtx &mtxOrg);

	//static void Test();

private:
	gsl_matrix *m_pGMtxOrg;
	gsl_matrix *m_pGMtxInv;
	gsl_permutation *m_pP;
	unsigned m_dim;
};

class Solve
{
public:
	Solve();
	~Solve();

	void Gen(Mtx &mtxX, Mtx &mtxA, Mtx &mtxB);

private:
	gsl_matrix *m_pGMtxOrg;
	gsl_permutation *m_pP;
	gsl_vector *m_pB;
	gsl_vector *m_pX;
	unsigned m_dim;
};

//*************************************************************************************************

class G_gslMtxOp
{
public:
	Inv				inv;
	Solve			sol;

private:
};

extern G_gslMtxOp gslMtxOp;

#endif