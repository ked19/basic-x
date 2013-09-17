#include "gslMtxOp.h"

G_gslMtxOp gslMtxOp;

//*************************************************************************************************

void CopyTo(const gsl_matrix &gsl, Mtx &mtx)
{
	Vect2D<unsigned> dim = mtx.GetDim();
	assert(gsl.size1==dim.m_y && gsl.size2==dim.m_x);

	for(unsigned y=0; y<dim.m_y; y++)
		for(unsigned x=0; x<dim.m_x; x++)
			mtx.CellRef(x, y) = (DATA)gsl_matrix_get(&gsl, y, x);
}

void CopyTo(const Mtx &mtx, gsl_matrix &gsl)
{
	Vect2D<unsigned> dim = mtx.GetDim();
	assert(gsl.size1==dim.m_y && gsl.size2==dim.m_x);

	for(unsigned y=0; y<dim.m_y; y++)
		for(unsigned x=0; x<dim.m_x; x++)
			gsl_matrix_set(&gsl, y, x, mtx.CellVal(x, y));
}

//***********************************************

void CopyTo(const Mtx &mtx, gsl_vector &gsl)
{
	Vect2D<unsigned> dim = mtx.GetDim();
	//assert(gsl.size1==dim.m_y && gsl.size2==dim.m_x);
	assert(dim.m_x == 1);

	for(unsigned y=0; y<dim.m_y; y++)
	{
		gsl_vector_set(&gsl, y, mtx.CellVal(0, y));
	}
}

void CopyTo(const gsl_vector &gsl, Mtx &mtx)
{
	Vect2D<unsigned> dim = mtx.GetDim();
	//assert(gsl.size1==dim.m_y && gsl.size2==dim.m_x);
	assert(dim.m_x == 1);

	for(unsigned y=0; y<dim.m_y; y++)
	{
		mtx.CellRef(0, y) = (DATA)gsl_vector_get(&gsl, y);
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

Inv::Inv()
	:m_dim(0), m_pGMtxOrg(0), m_pGMtxInv(0), m_pP(0)
{}

Inv::~Inv()
{
	if(m_pGMtxOrg)
	{
		gsl_matrix_free(m_pGMtxOrg);
		gsl_matrix_free(m_pGMtxInv);
		gsl_permutation_free(m_pP);
	}
}

void Inv::Gen(Mtx &mtxInv, Mtx &mtxOrg)
{
	Vect2D<unsigned> dimOrg = mtxOrg.GetDim();
	Vect2D<unsigned> dimInv = mtxInv.GetDim();
	assert(dimInv.m_x==dimOrg.m_x && dimInv.m_y==dimOrg.m_y);
	assert(dimOrg.m_x == dimOrg.m_y);

	if(m_dim != dimOrg.m_x)
	{
		if(m_pGMtxOrg)
		{
			gsl_matrix_free(m_pGMtxOrg);
			gsl_matrix_free(m_pGMtxInv);
			gsl_permutation_free(m_pP);
		}

		m_dim = dimOrg.m_x;
		m_pGMtxOrg = gsl_matrix_alloc(m_dim, m_dim);
		m_pGMtxInv = gsl_matrix_alloc(m_dim, m_dim);
		m_pP = gsl_permutation_alloc(m_dim);
	}
	else {}

	CopyTo(mtxOrg, *m_pGMtxOrg);

	int sign;
	gsl_linalg_LU_decomp(m_pGMtxOrg, m_pP, &sign);
	gsl_linalg_LU_invert(m_pGMtxOrg, m_pP, m_pGMtxInv);
	
	CopyTo(*m_pGMtxInv, mtxInv);
}

/*
void Inv::Test()
{
	cout << "start Inv testing" << endl;

	const unsigned dim = 2;
	DATA aaA[dim][dim] = {{ 4, 3}, {3,  2}};
	DATA aaB[dim][dim] = {{-2, 3}, {3, -4}};
	
	Mtx mtxIn(dim, dim);
	Mtx mtxOut(dim, dim);
	for(unsigned y=0; y<dim; y++)
	{
		for(unsigned x=0; x<dim; x++)
		{
			mtxIn.CellRef(x, y)  = aaA[y][x];
		}
	}

	mtxOp.inv.Gen(mtxOut, mtxIn);
	mtxOp.out << mtxIn;
	cout << endl;
	mtxOp.out << mtxOut;

	for(unsigned y=0; y<dim; y++)
	{
		for(unsigned x=0; x<dim; x++)
		{
			assert(myMath.IsEqual(mtxOut.CellVal(x, y), aaB[y][x]));
		}
	}

	cout << "Inv test ok" << endl << endl;
}
*/

//*************************************************************************************************
//
//*************************************************************************************************

Solve::Solve()
	:m_dim(0), m_pGMtxOrg(0), m_pP(0), m_pB(0), m_pX(0)
{}

Solve::~Solve()
{}

void Solve::Gen(Mtx &mtxX, Mtx &mtxA, Mtx &mtxB)
{
	Vect2D<unsigned> dimOrg = mtxA.GetDim();

	if(m_dim != dimOrg.m_x)
	{
		if(m_pGMtxOrg)
		{
			gsl_matrix_free(m_pGMtxOrg);
			gsl_permutation_free(m_pP);
			gsl_vector_free(m_pB);
			gsl_vector_free(m_pX);
		}

		m_dim = dimOrg.m_x;
		m_pGMtxOrg = gsl_matrix_alloc(m_dim, m_dim);
		m_pP = gsl_permutation_alloc(m_dim);
		m_pB = gsl_vector_alloc(m_dim);
		m_pX = gsl_vector_alloc(m_dim);
	}
	else {}

	CopyTo(mtxA, *m_pGMtxOrg);
	CopyTo(mtxB, *m_pB);

	int sign;
	gsl_linalg_LU_decomp(m_pGMtxOrg, m_pP, &sign);
	gsl_linalg_LU_solve(m_pGMtxOrg, m_pP, m_pB, m_pX);
	
	CopyTo(*m_pX, mtxX);
}