#include "denseMatrix.h"
//#include "matrixOperation.h"

Mtx2D::Mtx2D(const unsigned xDim, const unsigned yDim)
	:m_xDim(xDim), m_yDim(yDim), m_lineSize(xDim)
{
	m_pMemUsedNum = new unsigned;
	*m_pMemUsedNum = 1;

	unsigned cellNum = m_xDim * m_yDim;
	m_pMem = new DATA[cellNum];
	m_pCell = m_pMem;
}

Mtx2D::Mtx2D(Vect2D<unsigned> dim)
	:m_xDim(dim.m_x), m_yDim(dim.m_y), m_lineSize(dim.m_x)
{
	m_pMemUsedNum = new unsigned;
	*m_pMemUsedNum = 1;

	unsigned cellNum = m_xDim * m_yDim;
	m_pMem = new DATA[cellNum];
	m_pCell = m_pMem;
}

Mtx2D::Mtx2D(const Mtx2D &mtxOrg)
	:m_xDim(mtxOrg.m_xDim), m_yDim(mtxOrg.m_yDim), m_lineSize(mtxOrg.m_lineSize), m_pMem(mtxOrg.m_pMem)
{
	m_pMemUsedNum = mtxOrg.m_pMemUsedNum;
	++(*m_pMemUsedNum);

	m_pCell = mtxOrg.m_pCell;
}


Mtx2D::Mtx2D(const Mtx2D &mtxOrg, unsigned xOffset, unsigned yOffset, unsigned xDim, unsigned yDim)
	:m_xDim(xDim), m_yDim(yDim), m_lineSize(mtxOrg.m_lineSize), m_pMem(mtxOrg.m_pMem)
{
	m_pMemUsedNum = mtxOrg.m_pMemUsedNum;
	++(*m_pMemUsedNum);

	unsigned offset = yOffset*m_lineSize + xOffset;
	m_pCell = mtxOrg.m_pCell + offset;
}

Mtx2D::~Mtx2D()
{
	m_pCell = 0;

	if(*m_pMemUsedNum == 1)
	{
		delete []m_pMem;
	}
	else
	{
		--(*m_pMemUsedNum);
		m_pMemUsedNum = 0;
		m_pMem = 0;
	}
}

//*************************************************************************************************

bool Mtx2D::IsXInside(int x) const
{
	if(x>=0 && x<(int)m_xDim)	return true;
	else						return false;
}
bool Mtx2D::IsYInside(int y) const
{
	if(y>=0 && y<(int)m_yDim)	return true;
	else						return false;
}
bool Mtx2D::IsInside(int x, int y) const
{
	if(x>=0 && x<(int)m_xDim && y>=0 && y<(int)m_yDim)	return true;
	else												return false;
}

Vect2D<unsigned> Mtx2D::GetDim() const
{
	Vect2D<unsigned> dim(m_xDim, m_yDim);
	return dim;
}

DATA& Mtx2D::CellRef(unsigned x, unsigned y)
{
	unsigned idx = y*m_lineSize + x;
	return m_pCell[idx];
}
DATA Mtx2D::CellVal(unsigned x, unsigned y) const
{
	unsigned idx = y*m_lineSize + x;
	return 	m_pCell[idx];
}

void Mtx2D::CopyFrom(const Mtx2D &mtxFrom)
{
	assert(m_xDim==mtxFrom.m_xDim && m_yDim==mtxFrom.m_yDim);

	for(unsigned y=0; y<m_yDim; y++)
	{
		for(unsigned x=0; x<m_xDim; x++)
		{
			CellRef(x, y) = mtxFrom.CellVal(x, y);
		}
	}
}

void Mtx2D::CopyFrom(float *pData)
{
	unsigned loc = 0;
	for(unsigned y=0; y<m_yDim; y++)
	{
		for(unsigned x=0; x<m_xDim; x++)
		{
			CellRef(x, y) = pData[loc++];
		}
	}
}

void Mtx2D::CopyFrom(unsigned char *pData)
{
	unsigned loc = 0;
	for(unsigned y=0; y<m_yDim; y++) {
		for(unsigned x=0; x<m_xDim; x++) {
			CellRef(x, y) = pData[loc++];
		}
	}
}

void Mtx2D::CopyTo(DATA *pData) const
{
	unsigned mSize = m_xDim * m_yDim;
	unsigned loc = 0;
	for(unsigned y=0; y<m_yDim; y++)
	{
		for(unsigned x=0; x<m_xDim; x++)
		{
			pData[loc++] = CellVal(x, y);
		}
	}
}

void Mtx2D::CopyTo(float *pData) const
{
	unsigned mSize = m_xDim * m_yDim;
	unsigned loc = 0;
	for(unsigned y=0; y<m_yDim; y++)
	{
		for(unsigned x=0; x<m_xDim; x++)
		{
			pData[loc++] = (float)CellVal(x, y);
		}
	}
}

void Mtx2D::CopyTo(unsigned char *pData) const
{
	unsigned mSize = m_xDim * m_yDim;
	unsigned loc = 0;
	for(unsigned y = 0; y < m_yDim; y++) {
		for(unsigned x = 0; x < m_xDim; x++) {
			pData[loc++] = (unsigned char)CellVal(x, y);
		}
	}
}

/*
void Mtx2D::Test()
{
	cout << "start denseMatrix testing.." << endl; 

	const int xDim = 8;
	const int yDim = 4;
	float dataArray[yDim][xDim] = 
	{{ 0,  1,  2,  3,  4,  5,  6,  7}, 
	 { 8,  9, 10, 11, 12, 13, 14, 15}, 
	 {16, 17, 18, 19, 20, 21, 22, 23},
	 {24, 25, 26, 27, 28, 29, 30, 31}};

	mtxOp.out.SetWidth(3);

	Mtx mtx(xDim, yDim);
	for(unsigned y=0; y<yDim; y++)
	{
		mtx.CellRef(0, y) = dataArray[y][0];

		for(unsigned x=1; x<xDim; x++)
		{
			mtx.CellRef(x, y) = dataArray[y][x];
		}
	}
	cout << "matrix of input data:" << endl;
	mtxOp.out << mtx;

	Vect2D<unsigned> dimMtx = mtx.GetDim();
	assert(dimMtx.m_x == xDim);
	assert(dimMtx.m_y == yDim);
	for(unsigned y=0; y<dimMtx.m_y; y++)
	{
		for(unsigned x=0; x<dimMtx.m_x; x++)
		{
			assert(mtx.CellVal(x, y) == dataArray[y][x]);
		}
	}
	cout << "input / output ok" << endl << endl;
	
	//************************************************************************************************

	const unsigned xOffset = 2, yOffset = 1;
	const unsigned xSubDim = 3, ySubDim = 2;
	Mtx mtxSub(mtx, xOffset, yOffset, xSubDim, ySubDim);
	cout << "sub matrix with offset:(" << xOffset << ',' << yOffset << ')' << 
		                     "  dim:(" << xSubDim << ',' << ySubDim << ')' << endl;
	mtxOp.out << mtxSub;
	Vect2D<unsigned> dimSub = mtxSub.GetDim();
	assert(dimSub.m_x == xSubDim);
	assert(dimSub.m_y == ySubDim);
	for(unsigned y=0; y<dimSub.m_y; y++)
	{
		for(unsigned x=0; x<dimSub.m_x; x++)
		{
			assert(mtxSub.CellVal(x, y) == mtx.CellVal(x+xOffset, y+yOffset));
		}
	}

	cout << "set sub matrix(0,0) = 100" << endl;
	mtxSub.CellRef(0, 0) = 100.F;
	mtxOp.out << mtx;
	assert(mtx.CellVal(xOffset, yOffset) == 100.F);
	cout << "sub matrix ok" << endl << endl;

	//************************************************************************************************

	const unsigned x2Offset = 1, y2Offset = 0;
	const unsigned x2SubDim = 2, y2SubDim = 2;
	Mtx mtx2Sub(mtxSub, x2Offset, y2Offset, x2SubDim, y2SubDim);
	cout << "double sub matrix with offset:(" << x2Offset << ',' << y2Offset << ')' << 
									"  dim:(" << x2SubDim << ',' << y2SubDim << ')' << endl;
	mtxOp.out << mtx2Sub;
	Vect2D<unsigned> dim2Sub = mtx2Sub.GetDim();
	assert(dim2Sub.m_x == x2SubDim);
	assert(dim2Sub.m_y == y2SubDim);
	for(unsigned y=0; y<dim2Sub.m_y; y++)
	{
		for(unsigned x=0; x<dim2Sub.m_x; x++)
		{
			assert(mtx2Sub.CellVal(x, y) == mtxSub.CellVal(x+x2Offset, y+y2Offset));
		}
	}

	cout << "set double sub matrix(0,0) = 100" << endl;
	mtx2Sub.CellRef(0, 0) = 100.F;
	mtxOp.out << mtx;
	assert(mtx.CellVal(xOffset+x2Offset, yOffset+y2Offset) == 100.F);
	cout << "double sub matrix ok" << endl << endl;


	cout << "denseMatrix testing ok" << endl << endl;
}
*/

//*************************************************************************************************
//
//*************************************************************************************************