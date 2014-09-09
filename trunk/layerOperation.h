#ifndef _LAYER_OPERATION_H
#define _LAYER_OPERATION_H

#include "denseMatrix.h"
#include "define.h"
#include "matrixOperation.h"
#include "myMath.h"
#include "layer3D.h"
//#include "volumeData.h"

#include <fstream>
#include <fftw3.h>
#include <cmath>
#include <cassert>
#include <sstream>

using namespace std;

class ClampLyr
{
public:
	ClampLyr();
	~ClampLyr();

	void Gen(Layer &lyr, DATA thrd, bool bRmLow, DATA repV);

};

class BlendLyr
{
public:
	BlendLyr();
	~BlendLyr();

	void Gen(Layer &lyrOut, const Mtx &mtxIn, DATA ratio) const;
};

class SaveLyr
{
public:
	SaveLyr();
	~SaveLyr();

	void Gen(const Layer &lyr, string f, bool bB=true) const;
};
class SaveImgLyr
{
public:
	SaveImgLyr();
	~SaveImgLyr();

	void Gen(const Layer & lyr, string f, string ext) const;

//private:
//	stringstream ss;
//	string sN;
};

class ZeroLyr
{
public:
	ZeroLyr();
	~ZeroLyr();

	void Gen(Layer &lyr) const;

private:
	Zero m_zero;
};

class OneLyr
{
public:
	OneLyr();
	~OneLyr();

	void Gen(Layer &lyr) const;	
};

class RndLyr
{
public:
	RndLyr();
	~RndLyr();

	void Gen(Layer &lyr) const;
};

class Gauss3DLyr
{
public:
	Gauss3DLyr();
	~Gauss3DLyr();

	void Gen(Layer &lyr, DATA s, bool bNormal=true) const;
	void Gen(VolumeData &vol, DATA s, bool bNormal=true) const;
};

class Norm3DLyr
{
public:
	Norm3DLyr();
	~Norm3DLyr();

	DATA Gen(Layer &lyr, bool bNor=false) const;	
};

//*************************************************************************************************

// (0,0,0)
// dxdx  dxdy  dxdz
//
// dydx  dydy  dydz 
//
// dzdx  dzdy  dzdz
class HessianLyr
{
public:
	HessianLyr();
	~HessianLyr();

	void Gen(Mtx2D &mtxH, const Layer &lyr, unsigned x, unsigned y, unsigned c) const;

private:
};

//*************************************************************************************************

class DerivateLyr
{
public:
	DerivateLyr();
	~DerivateLyr();

	void Gen(Mtx &mtxDer, const Layer &lyr, unsigned x, unsigned y, unsigned c) const;
	void Gen(Layer &lyrDer, const Mtx &mtx, bool bNormal=true) const;

private:
};

//*************************************************************************************************

class Orient
{
public:
	Orient();
	~Orient();

	void Gen(Layer &lyrOrnt, const Mtx &mtx);

private:
};

//*************************************************************************************************

class ScaleDimMtxLyr
{
public:
	ScaleDimMtxLyr();
	~ScaleDimMtxLyr();

	enum METHOD {NEAREST, BILINEAR, BLOCK_AVG};
	void Gen(Layer &lyrOut, const Layer &lyrIn, METHOD m=BILINEAR);
};

class ScaleDimLyr
{
public:
	ScaleDimLyr();
	~ScaleDimLyr();

	enum METHOD	{NEAREST, TRLINEAR, BLOCK_AVG};
	void Gen(Layer &lyrOut, const Layer &lyrIn, METHOD m=TRLINEAR, Layer *pLyrBuf=0) const;

private:
	void Gen_nearest(Layer &lyrOut, const Layer &lyrIn) const;
	void Gen_bilinear(Layer &lyrOut, const Layer &lyrIn) const;
	//void Gen_blockAvg(Layer &lyrOut, const Mtx &mIn, Mtx *pmBuf) const;
};

//*************************************************************************************************

class Rgb2Hsv
{
public:
	Rgb2Hsv();
	~Rgb2Hsv();

	void Gen(Layer &lyrHsv, const Layer &lyrRgb, bool bNormH) const;
};

class Rgb2Ycbcr
{
public:
	Rgb2Ycbcr();
	~Rgb2Ycbcr();

	void Gen(Layer &lyrYcbcr, const Layer &lyrRgb) const;
};

class Gray
{
public:
	Gray();
	~Gray();

	void Gen(Mtx &mtxG, const Layer &lyrIn) const;
};

//*************************************************************************************************

class HisLyr
{
public:
	HisLyr();
	~HisLyr();

	//void Gen(Layer &lyrHis, const Layer &lyrIn, DATA aMaxV[], DATA aMinV[]) const;
	void Gen(Mtx &mtxHis, const Layer &lyrIn, DATA minV, DATA maxV) const;
	void Gen(Mtx &mtxHis, const Layer &lyrX, const Layer &lyrY, 
		DATA minX, DATA maxX, DATA minY, DATA maxY) const;
	
	//void SetBound(DATA aMaxV[], DATA aMinV[]);
	//void Add(Layer &lyrHis, DATA aVal[]) const;

private:
	DATA m_aMaxV[3];
	DATA m_aMinV[3];
	DATA m_aVLen[3];
	DATA m_aInvVLen[3];
};

//*************************************************************************************************

class MulLyr
{
public:
	MulLyr();
	~MulLyr();

	//void Gen(Layer &lyrR, const Layer &mtxA, const Mtx &mtxB) const;
	void Gen(Layer &lyr, DATA scl) const;

private:
};

class SubLyr
{
public:
	SubLyr();
	~SubLyr();

	void Gen(Layer &lyrA, const Layer &lyrB);
	void Gen(Layer &lyr, DATA v);
};

class AddLyr
{
public:
	AddLyr();
	~AddLyr();

	void Gen(Layer &lyr, DATA v);
	void Gen(Layer &lyr, const Mtx &mtx);
};

class CellMultiplyLyr
{
public:
	CellMultiplyLyr();
	~CellMultiplyLyr();

	void Gen(CplxLyr &clTo, const CplxLyr &clA, const CplxLyr &clB);
	void Gen(Layer &lyrTo, const Layer &lyrA, const Layer &lyrB);
};

//*************************************************************************************************

class AvgLyr
{
public:
	AvgLyr();
	~AvgLyr();

	DATA Gen(const Layer &lyr) const;

private:
};

class VarLyr
{
public:
	VarLyr();
	~VarLyr();

	DATA Gen(const Layer &lyr) const;

private:
};

class RngLyr
{
public:
	RngLyr();
	~RngLyr();

	void Gen(DATA &vMin, DATA &vMax, const Layer &lyr) const;
private:
};

//*************************************************************************************************

class FFTLyr
{
public:
	FFTLyr();
	~FFTLyr();

	void Gen(CplxLyr &cLyrTo, Layer &lyrFrom, bool bDividNum=false);
	void Gen(Layer &lyrTo, CplxLyr &cLyrFrom, bool bDividNum=false);
	//void Gen(CplxMtx &cMtxTo, CplxMtx &cMtxFrom, bool bDividNum=false);

	//static void Test();

private:
	enum T_CASE {R2C, C2R, C2C};
	void SetPlan(unsigned xDim, unsigned yDim, unsigned zDim, T_CASE tCase);

	void NewFMem(unsigned size);
	void NewTMem(unsigned size);
	void NewCFMem(unsigned size);
	void NewCTMem(unsigned size);

	T_CASE m_tCase;
	unsigned m_planXDim;
	unsigned m_planYDim;
	unsigned m_planZDim;
	fftw_plan m_plan;
	DATA *m_pFrom;					unsigned m_fSize;
	DATA *m_pTo;					unsigned m_tSize;
	fftw_complex *m_pCplxFrom;		unsigned m_cfSize;
	fftw_complex *m_pCplxTo;		unsigned m_ctSize;
};

class ConvLyr
{
public:
	ConvLyr();
	~ConvLyr();

	//void Gen(Layer &lyrOut, const Layer &lyrA, const Layer &lyrB);
	void Gen(Layer &lyrBase, const Layer &lryKerl);

private:
	void SetPlan();

	FFTLyr m_fft;
	FFTLyr m_ifft;
};

class DoGLyr
{
public:
	DoGLyr();
	~DoGLyr();

	void Gen(Layer &lyrIn, Layer &lyrGauss, DATA s, DATA sScl=sqrt(1.6F));
	void Gen(VolumeData &volIn, Layer &lyrGauss, DATA s, DATA sScl=sqrt(1.6F));
};

class SplattingLyr
{
public:
	SplattingLyr();
	~SplattingLyr();

	void New(Layer &lyrIn);
	void AddPoint(Layer &lyrIn, DATA val, const Layer &lyrW, unsigned x, unsigned y, unsigned z);
	void Gen(Layer &lyrIn);

private:
	Layer *m_pLyrCount;
};

//*************************************************************************************************

class G_LayerOp
{
public:
	SaveLyr				save;
	SaveImgLyr  		saveImg;		 		
	HessianLyr			Hess;
	ZeroLyr				zero;
	OneLyr				one;
	RndLyr 				rnd;
	Norm3DLyr 			norm;
	Gauss3DLyr			Gauss3D;
	CellMultiplyLyr		cellX;
	DerivateLyr			der;
	Orient				ornt;
	ScaleDimMtxLyr		scaleDimMtx;
	ScaleDimLyr			scaleDim;
	Rgb2Hsv				rgb2hsv;
	Rgb2Ycbcr			rgb2ycbcr;
	Gray				gray;
	HisLyr				his;
	MulLyr				mul;
	SubLyr				sub;
	AddLyr 				add;
	AvgLyr				avg;
	VarLyr				var;
	RngLyr				rng;
	ConvLyr				conv;
	DoGLyr				DoG;
	SplattingLyr		splt;
	BlendLyr			blend;
	ClampLyr 			clamp;

private:
};

extern G_LayerOp lyrOp;

#endif