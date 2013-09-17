#ifndef _MATRIX_OPERATION_H
#define _MATRIX_OPERATION_H

#include "denseMatrix.h"
#include "define.h"
#include "myMath.h"
#include "imageIO.h"

#include <fftw3.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <list>

using namespace std;

class Output
{
public:
	Output();
	~Output();

	Output& operator<<(const DATA &data);
	Output& operator<<(const Mtx &mtx);
	Output& operator<<(const CplxMtx &cMtx);

	void SetWidth(unsigned width);
	void SetPrecision(unsigned precision);

private:
	unsigned m_width;
	unsigned m_precision;
};

//*************************************************************************************************

class Zero
{
public:
	Zero();
	~Zero();

	void Gen(Mtx &mtx) const;

private:
};

class One
{
public:
	One();
	~One();

	void Gen(Mtx &mtx) const;
};

class I
{
public:
	I();
	~I();

	void Gen(Mtx &mtx) const;
	bool Is(const Mtx &mtx) const;

private:
};

class Norm
{
public:
	Norm();
	~Norm();

	DATA Gen(Mtx &mtx, bool bNor=false) const;

private:
};

class Abs
{
public:
	Abs();
	~Abs();

	void Gen(Mtx &mtx) const;

private:
};

//*************************************************************************************************

class CellMultiply
{
public:
	CellMultiply();
	~CellMultiply();

	void Gen(CplxMtx &cmTo, const CplxMtx &cmA, const CplxMtx &cmB);
	void Gen(Mtx &mtxTo, const Mtx &mtxA, const Mtx &mtxB);

private:
};

//*************************************************************************************************

class Gauss2D
{
public:
	Gauss2D();
	~Gauss2D();

	void Gen(Mtx &mtx, DATA devX, DATA devY, bool bNormal, DATA ang=0) const;

	static void Test();

private:
};

//*************************************************************************************************

class FFT
{
public:
	FFT();
	~FFT();

	void Gen(CplxMtx &cMtxTo, Mtx &mtxFrom, bool bDividNum=false);
	void Gen(Mtx &mtxTo, CplxMtx &cMtxFrom, bool bDividNum=false);
	void Gen(CplxMtx &cMtxTo, CplxMtx &cMtxFrom, bool bDividNum=false);

	static void Test();

private:
	enum T_CASE {R2C, C2R, C2C};
	void SetPlan(unsigned xDim, unsigned yDim, T_CASE tCase);

	void NewFMem(unsigned size);
	void NewTMem(unsigned size);
	void NewCFMem(unsigned size);
	void NewCTMem(unsigned size);

	T_CASE m_tCase;
	unsigned m_planXDim;
	unsigned m_planYDim;
	fftw_plan m_plan;
	DATA *m_pFrom;					unsigned m_fSize;
	DATA *m_pTo;					unsigned m_tSize;
	fftw_complex *m_pCplxFrom;		unsigned m_cfSize;
	fftw_complex *m_pCplxTo;		unsigned m_ctSize;
};

class Conv
{
public:
	Conv();
	~Conv();

	void Gen(Mtx &mtxBase, const Mtx &mtxKerl);

	static void Test();

private:
	void SetPlan();

	FFT m_fft;
	FFT m_ifft;
};

//*************************************************************************************************

class Sum
{
public:
	Sum();
	~Sum();

	void Gen(Mtx &mOut, const Mtx &mIn) const;
	DATA Gen(const Mtx &mIn) const;

private:;
};

class Avg
{
public:
	Avg();
	~Avg();

	DATA Gen(const Mtx &mIn) const;

private:
};

class Dev
{
public:
	Dev();
	~Dev();

	DATA Gen(const Mtx &mIn, bool bDivN=true) const;
};

class Rng
{
public:
	Rng();
	~Rng();

	void Gen(DATA &vMin, DATA &vMax, const Mtx &m) const;
};
//*************************************************************************************************

class ScaleDim
{
public:
	ScaleDim();
	~ScaleDim();

	enum METHOD {NEAREST, BILINEAR, BLOCK_AVG};
	void Gen(Mtx &mOut, const Mtx &mIn, METHOD m=BILINEAR, Mtx *pmBuf=0) const;

private:
	void Gen_nearest(Mtx &mOut, const Mtx &mIn) const;
	void Gen_bilinear(Mtx &mOut, const Mtx &mIn) const;
	void Gen_blockAvg(Mtx &mOut, const Mtx &mIn, Mtx *pmBuf) const;
};

//*************************************************************************************************

class Add
{
public:
	Add();
	~Add();

	void Gen(Mtx &mtxBase, const Mtx &adder) const;

private:
};

class Sub
{
public:
	Sub();
	~Sub();

	void Gen(Mtx &mtxBase, const Mtx &suber) const;
	void Gen(DATA cNum, Mtx &mtxSuber) const;

private:
};

class Mul
{
public:
	Mul();
	~Mul();

	void Gen(Mtx &mtxR, const Mtx &mtxA, const Mtx &mtxB) const;
	void Gen(Mtx &mtx, DATA scl) const;

	static void Test();

private:
};


//*************************************************************************************************

class Hessian
{
public:
	Hessian();
	~Hessian();

	void Gen(Mtx2D &mtxH, const Mtx &mtx, unsigned x, unsigned y) const;

private:
};

//*************************************************************************************************

class Derivate
{
public:
	Derivate();
	~Derivate();

	void Gen(Mtx &mtxDer, const Mtx &mtx, unsigned x, unsigned y, bool bNormal=true) const;

private:
};

//*************************************************************************************************

class ValCount
{
public:
	DATA m_val;
	DATA m_count;
};

//*************************************************************************************************

class Median
{
public:
	Median();
	~Median();

	DATA Gen(DATA *pData, unsigned size);
};

//*************************************************************************************************

class His
{
public:
	His();
	~His();

	void Gen(Mtx &mtxHis, Mtx &mtxIn, DATA maxV, DATA minV);

	void SetBound(DATA maxV, DATA minV);
	void Add(Mtx &mtxHis, DATA val) const;

private:
	DATA m_maxV;
	DATA m_minV;
	DATA m_vLen;
	DATA m_invVLen;
};

class HisEqual
{
public:
	HisEqual();
	~HisEqual();

	void Gen(Mtx &mtx, unsigned levlNum, bool bReval=true, Mtx *pMSkip=0);
	void Reval(Mtx &mtx, Mtx *pMSkip=0);

private:
	unsigned m_levelNum;
	unsigned m_memSize;
	unsigned *m_pCount;
};

//*************************************************************************************************

class Distance
{
public:
	Distance();
	~Distance();

	enum METHOD {SAD, SSD};
	DATA Gen(Mtx &mtxA, Mtx &mtxB, METHOD m=SSD);
};

//*************************************************************************************************

class Thrd
{
public:
	Thrd();
	~Thrd();

	void Gen(Mtx &mtx, DATA thrd);
};

class Otsu
{
public:
	Otsu();
	~Otsu();

	DATA Gen(Mtx &mtx, DATA &thrd, unsigned lev=256, DATA valMax=255.F, DATA valMin=0);
};

//*************************************************************************************************

class MedFlt
{
public:
	MedFlt();
	~MedFlt();

	void Gen(Mtx &mtxOut, Mtx &mtxIn, unsigned winLen);
};

//*************************************************************************************************

class Dilate
{
public:
	Dilate();
	~Dilate();

	void Gen(Mtx &mtxOut, Mtx &mtxIn, unsigned len);
};
class Erose
{
public:
	Erose();
	~Erose();

	void Gen(Mtx &mtxOut, Mtx &mtxIn, unsigned len);
};

//*************************************************************************************************

class DoG
{
public:
	DoG();
	~DoG();

	void Gen(Mtx &mtxIn, Mtx* apMtxG[], DATA s, DATA wSub=0.98F, DATA sScl=sqrt(1.6));
	void GenNPR(Mtx &mtxIn, Mtx* apMtxG[], DATA s, DATA coeSharp);
};

//*************************************************************************************************

class CrossCorr
{
public:
	CrossCorr();
	~CrossCorr();

	DATA Gen(const Mtx &mtxIn, Mtx &mtxMsk, unsigned xOff, unsigned yOff);
};

//*************************************************************************************************

void Rotate2D(DATA aOut[], const DATA aIn[], DATA ang, const DATA aCnt[]);
DATA Length2D(DATA x, DATA y);

//*************************************************************************************************

class G_MtxOp
{
public:
	Output			out;
	Zero			zero;
	One				one;
	I				I;
	Abs				abs;
	Norm			norm;
	CellMultiply	cellX;
	Gauss2D			Gauss;
	FFT				fft;
	FFT				ifft;
	Conv			conv;
	Sum				sum;
	Avg				avg;
	Dev				dev;
	Rng				rng;
	CrossCorr		crossCorr;
	ScaleDim		scaleDim;
	Add				add;
	Sub				sub;		
	Mul				mul;
	//Inv				inv;
	Hessian			Hess;
	Derivate		der;
	Median			median;
	MedFlt			medFlt;
	His				his;
	HisEqual		hisEqual;
	Distance		dis;
	Thrd			thrd;
	Otsu			Otsu;
	Dilate			dilate;
	Erose			erose;
	DoG				DoG;

private:
};

extern G_MtxOp mtxOp;

#endif