#ifndef _DENSE_MATRIX_H
#define _DENSE_MATRIX_H

#include "define.h"

#include <vector>
#include <memory>
#include <cassert>
#include <limits>

using namespace std;

class Mtx2D
{
public:
	Mtx2D(const unsigned xDim, const unsigned yDim);
	Mtx2D(Vect2D<unsigned> dim);
	Mtx2D(const Mtx2D &mtxOrg);
	Mtx2D(const Mtx2D &mtxOrg, unsigned xOffset, unsigned yOffset, unsigned xDim, unsigned yDim);
	~Mtx2D();

	bool IsXInside(int x) const;
	bool IsYInside(int y) const;
	bool IsInside(int x, int y) const;
	Vect2D<unsigned> GetDim() const;
	
	DATA& CellRef(unsigned x, unsigned y);
	DATA  CellVal(unsigned x, unsigned y) const;

	void CopyFrom(const Mtx2D &mtxFrom);
	void CopyFrom(float *pData);
	void CopyFrom(unsigned char *pData);
	void CopyTo(DATA *pData) const;
	void CopyTo(float *pData) const;
	void CopyTo(unsigned char* pData) const;

	static void Test();

private:
	const unsigned m_xDim;
	const unsigned m_yDim;
	const unsigned m_lineSize;

	unsigned *m_pMemUsedNum;
	DATA *m_pCell;
	DATA *m_pMem;
};

typedef Mtx2D Mtx;

//*************************************************************************************************



#endif