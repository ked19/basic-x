#ifndef _LAYER_H
#define _LAYER_H

#include "denseMatrix.h"
#include <vector>

using namespace std;

class Layer2D
{
public:
	Layer2D(const unsigned xDim, const unsigned yDim, const unsigned cDim);
	Layer2D(const Layer2D &layOrg);
	Layer2D(const Layer2D &layOrg, unsigned xOffset, unsigned yOffset, unsigned cOffset, unsigned xDim, unsigned yDim, unsigned cDim);
	Layer2D(const Mtx &mtx);
	Layer2D(const Mtx *ptrMtx, unsigned num);
	~Layer2D();

	bool IsXInside(int x) const;
	bool IsYInside(int y) const;
	bool IsCInside(int c) const;
	bool IsInside(int x, int y, int c) const;
	void ClipBound(int &xL, unsigned &xR, int &yB, unsigned &yT, int &zN, unsigned &zF);
	Vect3D<unsigned> GetDim() const;
	
	DATA& CellRef(unsigned x, unsigned y, unsigned c);
	DATA  CellVal(unsigned x, unsigned y, unsigned c) const;

	void CopyFrom(const Layer2D &lyrFrom);
	void CopyFrom_zFirst(float *pFrom);
	void CopyTo_zFirst(float *pTo, bool bNormal=false) const;
	void CopyTo_zLast(float *pTo) const;

	Mtx* GetMtx(unsigned c) const;
	//Layer2D* GetLayer2D() const;

	Layer2D GetRgb();

private:
	vector<Mtx*> m_vpMtx; 

protected:
	const unsigned m_xDim;
	const unsigned m_yDim;
	const unsigned m_cDim;
};

typedef Layer2D		Layer;
//typedef	Layer		MyImg

//*************************************************************************************************

/*
#include "layer.h"
class Layer2D;
*/

class CplxMtx2D
{
public:
	CplxMtx2D(unsigned xDim, unsigned yDim);
	CplxMtx2D(const CplxMtx2D &cMtxOrg, unsigned xOffset=0, unsigned yOffset=0, unsigned xDim=0, unsigned yDim=0);
	~CplxMtx2D();

	bool IsXInside(int x) const;
	bool IsYInside(int y) const;
	bool IsInside(int x, int y) const;
	Vect2D<unsigned> GetDim() const;

	Vect2D<DATA&> CellRef(unsigned x, unsigned y);
	Vect2D<DATA>  CellVal(unsigned x, unsigned y) const;

	Mtx* GetRealMtx() const;
	Mtx* GetImgMtx() const;

	static void Test();

private:
	const unsigned m_xDim;
	const unsigned m_yDim;

	Layer2D *m_pLay;
};

typedef CplxMtx2D CplxMtx;

//*************************************************************************************************

class MyImg: public Layer2D
{
public:
	MyImg(unsigned xDim, unsigned yDim, unsigned cDim);
	MyImg(const MyImg &imgOrg);
	MyImg(const MyImg &imgOrg, unsigned xOffset, unsigned yOffset, unsigned xDim, unsigned yDim);
	MyImg(const Mtx &mtx);
	MyImg(const Layer &lyr);
	~MyImg();

	MyImg* ConvertGray();

private:
};

//*************************************************************************************************

class VolumeData
{
public:
	VolumeData(Layer &lyrData, Vect3D<DATA> step);
	VolumeData(Vect3D<unsigned> dim, Vect3D<DATA> step);
	~VolumeData();

	Vect3D<unsigned> GetDim() const;
	Vect3D<DATA> GetStep() const;
	DATA& CellRef(unsigned x, unsigned y, unsigned c);
	DATA  CellVal(unsigned x, unsigned y, unsigned c) const;
	Layer& GetLyrRef();
	const Layer& GetLyrVal() const;

private:
	Vect3D<DATA> m_vStep;
	Layer m_lyrData;
};

#endif