#ifndef _LAYER_3D_H
#define _LAYER_3D_H

#include "layer.h"
#include <vector>

using namespace std;

class Layer3D
{
public:
	Layer3D(unsigned xDim, unsigned yDim, unsigned zDim, unsigned cDim);
	Layer3D(const Layer3D &lay3Org);
	Layer3D(const Layer3D &lay3Org, unsigned xOffset, unsigned yOffset, unsigned zOffset, unsigned cOffset, 
		unsigned xDim, unsigned yDim, unsigned zDim, unsigned cDim);
	//Layer3D(const Mtx &mtx);
	//Layer3D(const Mtx *ptrMtx, unsigned num);
	~Layer3D();

	bool IsXInside(int x) const;
	bool IsYInside(int y) const;
	bool IsZInside(int z) const;
	bool IsCInside(int c) const;
	bool IsInside(int x, int y, int z, int c) const;
	Vect4D<unsigned> GetDim() const;
	
	DATA& CellRef(unsigned x, unsigned y, unsigned z, unsigned c);
	DATA  CellVal(unsigned x, unsigned y, unsigned z, unsigned c) const;

	void CopyFrom(const Layer3D &lyr3From);
	//void CopyTo_zFirst(float *pTo, bool bNormal=false) const;
	//void CopyTo_zLast(float *pTo) const;

	Layer* GetLyr(unsigned c) const;
	//Mtx* GetMtx(unsigned c) const;
	//Layer2D* GetLayer2D() const;

private:
	vector<Layer2D*> m_vpLyr; 

protected:
	const unsigned m_xDim;
	const unsigned m_yDim;
	const unsigned m_zDim;
	const unsigned m_cDim;
};

//typedef Layer2D		Layer;
//typedef	Layer		MyImg

//*************************************************************************************************

class CplxLyr2D
{
public:
	CplxLyr2D(unsigned xDim, unsigned yDim, unsigned zDim);
	CplxLyr2D(const CplxLyr2D &cLyrOrg, unsigned xOffset=0, unsigned yOffset=0, unsigned zOffset=0, unsigned xDim=0, unsigned yDim=0, unsigned zDim=0);
	~CplxLyr2D();

	bool IsXInside(int x) const;
	bool IsYInside(int y) const;
	bool IsZInside(int z) const;
	bool IsInside(int x, int y, int z) const;
	Vect3D<unsigned> GetDim() const;

	Vect2D<DATA&> CellRef(unsigned x, unsigned y, unsigned z);
	Vect2D<DATA>  CellVal(unsigned x, unsigned y, unsigned z) const;

	Layer* GetRealLyr() const;
	Layer* GetImgLyr() const;

	//static void Test();

private:
	const unsigned m_xDim;
	const unsigned m_yDim;
	const unsigned m_zDim;

	Layer3D *m_pLyr3;
};

typedef CplxLyr2D CplxLyr;

#endif