#include "matrixOperation.h"

G_MtxOp mtxOp;

//*************************************************************************************************

Output::Output()
{
	m_width		= (unsigned)cout.width();
	m_precision = (unsigned)cout.precision(); 
}

Output::~Output()
{}

Output& Output::operator<<(const DATA &data) 
{
	cout.width(m_width);
	cout.precision(m_precision);
	cout << data;
	return *this;
}

Output& Output::operator<<(const Mtx &mtx)
{
	cout.width(m_width);
	cout.precision(m_precision);
	Vect2D<unsigned> dim = mtx.GetDim();
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			cout << mtx.CellVal(x, y);
			cout << " ";
		}  
		cout << endl;
	}  
	return *this;
}

Output& Output::operator <<(const CplxMtx &cMtx)
{
	cout.width(m_width);
	cout.precision(m_precision);
	Vect2D<unsigned> dim = cMtx.GetDim();
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			Vect2D<DATA> val = cMtx.CellVal(x, y);
			cout << val.m_x;	
			if(val.m_y>=0)	cout << "+";
			else			cout << "-";
			cout << fabs(val.m_y) << "i"; 
			cout << " ";
		}
		cout << endl;
	}
	return *this;
}

void Output::SetWidth(unsigned width)
{
	m_width = width;
}
void Output::SetPrecision(unsigned precision)
{
	m_precision = precision;
}

//*************************************************************************************************
//
//*************************************************************************************************

Zero::Zero()
{}

Zero::~Zero()
{}

void Zero::Gen(Mtx &mtx) const
{
	Vect2D<unsigned> dim = mtx.GetDim();
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			mtx.CellRef(x, y) = 0;
		}
	}
}

One::One()
{}

One::~One()
{}

void One::Gen(Mtx &mtx) const
{
	Vect2D<unsigned> dim = mtx.GetDim();
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			mtx.CellRef(x, y) = 1.F;
		}
	}
}

I::I()
{}

I::~I()
{}

void I::Gen(Mtx &mtx) const
{
	mtxOp.zero.Gen(mtx);

	Vect2D<unsigned> dim = mtx.GetDim();
	assert(dim.m_x == dim.m_y);
	for(unsigned i=0; i<dim.m_x; i++)
	{
		mtx.CellRef(i, i) = 1.F;
	}
}

bool I::Is(const Mtx &mtx) const
{
	Vect2D<unsigned> dim = mtx.GetDim();
	if(dim.m_x != dim.m_y)		return false;
	else{}

	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			DATA val = (x==y)? 1.F: 0;
			if(myMath.IsEqual(mtx.CellVal(x, y), val))		return false;
			else {}
		}
	}
	return true;
}

Norm::Norm()
{}

Norm::~Norm()
{}

DATA Norm::Gen(Mtx &m, bool bNor) const
{
	Vect2D<unsigned> dim = m.GetDim();

	DATA len = 0;
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			len += m.CellVal(x, y) * m.CellVal(x, y);
		}
	}
	len = sqrt(len);

	if(bNor)
	{
		if(myMath.IsEqual(len, 0))
		{
			return len;
		}
		else
		{
			for(unsigned y=0; y<dim.m_y; y++)
			{
				for(unsigned x=0; x<dim.m_x; x++)
				{
					m.CellRef(x, y) /= len;
				}
			}
			return len;
		}
	}
	else 
	{
		return len;
	}
}

Abs::Abs()
{}

Abs::~Abs()
{}

void Abs::Gen(Mtx &mtx) const
{
	Vect2D<unsigned> dim = mtx.GetDim();

	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			mtx.CellRef(x, y) = fabs(mtx.CellVal(x, y));
		}
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

CellMultiply::CellMultiply()
{}

CellMultiply::~CellMultiply()
{}

void CellMultiply::Gen(CplxMtx &cmTo, const CplxMtx &cmA, const CplxMtx &cmB)
{
	Vect2D<unsigned> dimA  = cmA.GetDim();
	Vect2D<unsigned> dimB  = cmB.GetDim();
	Vect2D<unsigned> dimTo = cmTo.GetDim();
	assert(dimA.m_x==dimB.m_x && dimA.m_y==dimB.m_y);
	assert(dimTo.m_x==dimA.m_x && dimTo.m_y==dimA.m_y);

	for(unsigned y=0; y<dimA.m_y; y++)
	{
		for(unsigned x=0; x<dimA.m_x; x++)
		{
			Vect2D<DATA> valA = cmA.CellVal(x, y);
			Vect2D<DATA> valB = cmB.CellVal(x, y);
			Vect2D<DATA&> valTo = cmTo.CellRef(x, y);
			valTo.m_x = valA.m_x*valB.m_x - valA.m_y*valB.m_y;
			valTo.m_y = valA.m_x*valB.m_y + valA.m_y*valB.m_x;
		}
	}
}

void CellMultiply::Gen(Mtx &mtxTo, const Mtx &mtxA, const Mtx &mtxB)
{
	Vect2D<unsigned> dimA  = mtxA.GetDim();
	Vect2D<unsigned> dimB  = mtxB.GetDim();
	Vect2D<unsigned> dimTo = mtxTo.GetDim();
	assert(dimA.m_x==dimB.m_x && dimA.m_y==dimB.m_y);
	assert(dimTo.m_x==dimA.m_x && dimTo.m_y==dimA.m_y);

	for(unsigned y=0; y<dimA.m_y; y++)
	{
		for(unsigned x=0; x<dimA.m_x; x++)
		{
			mtxTo.CellRef(x, y) = mtxA.CellVal(x, y) * mtxB.CellVal(x, y);
		}
	}
}

void CellMultiply::Gen(Mtx &mtxTo, const Mtx &mtxB)
{
	Vect2D<unsigned> dimTo = mtxTo.GetDim();
	Vect2D<unsigned> dimB  = mtxB.GetDim();
	MyAssert(dimTo.m_x == dimB.m_x && 
			 dimTo.m_y == dimB.m_y);

	for(unsigned y = 0; y < dimTo.m_y; y++) {
		for(unsigned x = 0; x < dimTo.m_x; x++) {
			mtxTo.CellRef(x, y) *= mtxB.CellVal(x, y);
		}
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

Gauss2D::Gauss2D()
{}

Gauss2D::~Gauss2D()
{}

DATA ComputeGVal_0(int x, int y, DATA halfInvVarX, DATA halfInvVarY, DATA ang) 
{
	// g = e^-(XX/2 + YY/2)
	//
	// X = x / dx
	//
	// Y = y / dy

	DATA expVal = -(x*x*halfInvVarX + y*y*halfInvVarY);
	return exp(expVal);
}

DATA ComputeGVal_ang(int x, int y, DATA halfInvVarX, DATA halfInvVarY, DATA ang)
{
	// g = e^-(axx + 2bxy + cyy)
	//
	// a = [(cosA/dx)^2 + (sinA/dy)^2] / 2
	//
	// b = [-sin2A/dx^2 + sin2A/dy^2] / 4
	//
	// c = [(sinA/dx)^2 + (cosA/dy)^2] / 2

	DATA cosA = (DATA)cos(ang);		DATA powCosA = cosA * cosA;
	DATA sinA = (DATA)sin(ang);		DATA powSinA = sinA * sinA;
	DATA sin2A = (DATA)sin(2.F*ang);
	
	DATA xxCoef = powCosA*halfInvVarX + powSinA*halfInvVarY;
	DATA yyCoef = powSinA*halfInvVarX + powCosA*halfInvVarY;
	DATA xyCoef = -sin2A*halfInvVarX + sin2A*halfInvVarY;

	DATA expVal = -(x*x*xxCoef + y*y*yyCoef + x*y*xyCoef);
	return exp(expVal);
}

void Gauss2D::Gen(Mtx &mtx, DATA devX, DATA devY, bool bNormal, DATA ang) const
{
	assert(devX != 0);
	assert(devY != 0);

	Vect2D<unsigned> dim = mtx.GetDim();
	int yT = (dim.m_y - 1) / 2;		int yB = -yT;
	int xR = (dim.m_x - 1) / 2;		int xL = -xR;

	DATA varX = devX * devX;	DATA halfInvVarX = myMath.IsEqual(varX, 0)?	0:	0.5F / varX;
	DATA varY = devY * devY;	DATA halfInvVarY = myMath.IsEqual(varY, 0)?	0:	0.5F / varY;

	DATA (*ComputeGVal)(int, int, DATA, DATA, DATA);
	ComputeGVal = myMath.IsEqual(ang, 0)? ComputeGVal_0: ComputeGVal_ang;

	Zero zero;
	zero.Gen(mtx);
	for(int y=yB, yLoc=0; y<=yT; y++, yLoc++)
	{
		for(int x=xL, xLoc=0; x<=xR; x++, xLoc++)
		{
			mtx.CellRef(xLoc, yLoc) = ComputeGVal(x, y, halfInvVarX, halfInvVarY, ang);
		}
	}

	if(bNormal)
	{
		DATA sum = 0;
		DATA normalRatio = 1.F / (2.F * PI * devX * devY);
		for(unsigned y=0; y<dim.m_y; y++)
		{
			for(unsigned x=0; x<dim.m_x; x++)
			{
				mtx.CellRef(x, y) *= normalRatio;
				sum += mtx.CellRef(x, y);
			}
		}
		cout << "sum of Gauss window is: " << sum << endl;
		//assert(myMath.IsEqual(1.f, sum));
	}
	else
	{}
}

void Gauss2D::Test()
{
	cout << "start Gauss2D testing" << endl;

	const unsigned width = 5;
	cout.width(width);
	mtxOp.out.SetWidth(width);

	double aRealGVal[] = {0.13534, 0.60653, 1.00000, 0.60653, 0.13534};  // dev = 1.0
	const unsigned len = sizeof(aRealGVal) / sizeof(double);
	cout << "real value:\n";
	for(unsigned i=0; i<len; i++)	cout << aRealGVal[i] << " ";
	cout << endl;

	const DATA dev = 1.F;
	const DATA d2 = 0.00001F;
	Mtx vectX(5, 1), vectX2(6, 1);
	Mtx vectY(1, 5), vectY2(1, 6);
	Gauss2D gGen;
	gGen.Gen(vectX,  dev, d2,  false);
	gGen.Gen(vectX2, dev, d2,  false);
	gGen.Gen(vectY,	 d2,  dev, false);
	gGen.Gen(vectY2, d2,  dev, false);

	cout << "vect  x:\n";	mtxOp.out << vectX;	cout << endl;
	cout << "vect x2:\n";	mtxOp.out << vectX2;	cout << endl;
	cout << "vect  y:\n";	mtxOp.out << vectY;	cout << endl;
	cout << "vect y2:\n";	mtxOp.out << vectY2;	cout << endl;

	cout << "pass? (n/y):";
	char ny;
	cin >> ny;
	assert(ny == 'y');
	cout << endl;

	cout << "test normal distribution" << endl;
	Mtx m(10, 10);
	gGen.Gen(m, 1.F, 1.F, true);
	cout << "normal distribution ok" << endl << endl;

	cout << "Gauss2D testing ok" << endl << endl;
}

//*************************************************************************************************
//
//*************************************************************************************************

FFT::FFT()
	:m_pFrom(0), m_pTo(0),   m_pCplxFrom(0), m_pCplxTo(0) 
	,m_fSize(0), m_tSize(0), m_cfSize(0),    m_ctSize(0)
	,m_tCase(R2C), m_planXDim(0), m_planYDim(0), m_plan(0)
{}

FFT::~FFT()
{
	fftw_destroy_plan(m_plan);
	fftw_free(m_pFrom);	
	fftw_free(m_pTo);		
	fftw_free(m_pCplxFrom);
	fftw_free(m_pCplxTo);
}

void FFT::NewFMem(unsigned size)
{
	if(m_fSize >= size)
	{
		return;
	}
	else
	{
		fftw_free(m_pFrom);
		m_fSize = size;
		m_pFrom = (double*)fftw_malloc(sizeof(double) * m_fSize);
	}
}
void FFT::NewTMem(unsigned size)
{
	if(m_tSize >= size)
	{
		return;
	}
	else
	{
		fftw_free(m_pTo);
		m_tSize = size;
		m_pTo = (double*)fftw_malloc(sizeof(double) * m_tSize);
	}
}
void FFT::NewCFMem(unsigned size)
{
	if(m_cfSize >= size)
	{
		return;
	}
	else
	{
		fftw_free(m_pCplxFrom);
		m_cfSize = size;
		m_pCplxFrom = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * m_cfSize);
	}
}
void FFT::NewCTMem(unsigned size)
{
	if(m_ctSize >= size)
	{
		return;
	}	
	else
	{
		fftw_free(m_pCplxTo);
		m_ctSize = size;
		m_pCplxTo = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * m_ctSize);
	}
}

void FFT::SetPlan(unsigned xDim, unsigned yDim, T_CASE tCase)
{
	if(xDim==m_planXDim && yDim==m_planYDim && tCase ==m_tCase)	return;
	else {}

	fftw_destroy_plan(m_plan);

	unsigned rSize = xDim * yDim;
	unsigned cSize = (xDim/2 + 1) * yDim; 
	if(tCase == C2C)
	{
		NewCFMem(rSize);
		NewCTMem(rSize);
		m_plan = fftw_plan_dft_2d(yDim, xDim, m_pCplxFrom, m_pCplxTo, FFTW_FORWARD, FFTW_ESTIMATE);
	}
	else if(tCase == C2R)
	{
		NewCFMem(cSize);
		NewTMem (rSize);
		m_plan = fftw_plan_dft_c2r_2d(yDim, xDim, m_pCplxFrom, m_pTo, FFTW_ESTIMATE);
	}
	else if(tCase == R2C)
	{
		NewFMem (rSize);
		NewCTMem(cSize);
		m_plan = fftw_plan_dft_r2c_2d(yDim, xDim, m_pFrom, m_pCplxTo, FFTW_ESTIMATE);
	}
	else assert(0);
	m_tCase = tCase;
	m_planXDim = xDim;
	m_planYDim = yDim;
}

void FFT::Gen(CplxMtx &cMtxTo, CplxMtx &cMtxFrom, bool bDividNum)
{
	Vect2D<unsigned> dimFrom = cMtxFrom.GetDim(); 
	Vect2D<unsigned> dimTo   = cMtxTo.GetDim();
	assert(dimFrom.m_x==dimTo.m_x && dimFrom.m_y==dimTo.m_y);

	SetPlan(dimFrom.m_x, dimFrom.m_y, C2C);

	unsigned idx = 0;
	for(unsigned y=0; y<m_planYDim; y++)
	{
		for(unsigned x=0; x<m_planXDim; x++)
		{
			Vect2D<DATA> valFrom = cMtxFrom.CellVal(x, y);
			m_pCplxFrom[idx][0] = valFrom.m_x;
			m_pCplxFrom[idx][1] = valFrom.m_y;
		}
	}

	fftw_execute(m_plan);

	idx = 0;
	for(unsigned y=0; y<m_planYDim; y++)
	{
		for(unsigned x=0; x<m_planXDim; x++)
		{
			Vect2D<DATA&> valTo = cMtxTo.CellRef(x, y);
			valTo.m_x = m_pCplxTo[idx][0];
			valTo.m_y = m_pCplxTo[idx][1];
			++idx;
		}
	}

	if(bDividNum)
	{
		unsigned size = m_planXDim * m_planYDim;
		DATA invN = 1.F / size;
		for(unsigned y=0; y<m_planYDim; y++)
		{
			for(unsigned x=0; x<m_planXDim; x++)
			{
				Vect2D<DATA&> val = cMtxTo.CellRef(x, y);
				val.m_x *= invN;
				val.m_y *= invN;
			}
		}
	}
	else
	{}
}

void FFT::Gen(CplxMtx &cMtxTo, Mtx &mtxFrom, bool bDividNum)
{
	Vect2D<unsigned> dimFrom = mtxFrom.GetDim(); 
	Vect2D<unsigned> dimTo   = cMtxTo.GetDim();
	unsigned planCXDim = dimFrom.m_x/2 + 1;
	unsigned planCYDim = dimFrom.m_y;
	assert(dimTo.m_x==planCXDim && dimTo.m_y==planCYDim);

	SetPlan(dimFrom.m_x, dimFrom.m_y, R2C);

	unsigned idx = 0;
	for(unsigned y=0; y<m_planYDim; y++)
	{
		for(unsigned x=0; x<m_planXDim; x++)
		{
			m_pFrom[idx] = mtxFrom.CellVal(x, y);
			++idx;
		}
	}

	fftw_execute(m_plan);

	idx = 0;
	for(unsigned y=0; y<planCYDim; y++)
	{
		for(unsigned x=0; x<planCXDim; x++)
		{
			Vect2D<DATA&> valTo = cMtxTo.CellRef(x, y);
			valTo.m_x = m_pCplxTo[idx][0] ;
			valTo.m_y = m_pCplxTo[idx][1];
			++idx;
		}
	}

	if(bDividNum)
	{
		unsigned size = m_planXDim * m_planYDim;
		DATA invN = 1.F / size;
		for(unsigned y=0; y<planCYDim; y++)
		{
			for(unsigned x=0; x<planCXDim; x++)
			{
				Vect2D<DATA&> val = cMtxTo.CellRef(x, y);
				val.m_x *= invN;
				val.m_y *= invN;
			}
		}
	}
	else
	{}
}

void FFT::Gen(Mtx &mtxTo, CplxMtx &cMtxFrom, bool bDividNum)
{
	Vect2D<unsigned> dimFrom = cMtxFrom.GetDim(); 
	Vect2D<unsigned> dimTo   = mtxTo.GetDim();
	unsigned planCXDim = dimTo.m_x/2 + 1;
	unsigned planCYDim = dimTo.m_y;
	assert(dimFrom.m_x==planCXDim  && dimFrom.m_y==planCYDim);

	SetPlan(dimTo.m_x, dimTo.m_y, C2R);

	unsigned idx = 0;
	for(unsigned y=0; y<planCYDim; y++)
	{
		for(unsigned x=0; x<planCXDim; x++)
		{
			Vect2D<DATA> valFrom = cMtxFrom.CellVal(x, y);
			m_pCplxFrom[idx][0] = valFrom.m_x;
			m_pCplxFrom[idx][1] = valFrom.m_y;
			++idx;
		}
	}

	fftw_execute(m_plan);

	idx = 0;
	for(unsigned y=0; y<m_planYDim; y++)
	{
		for(unsigned x=0; x<m_planXDim; x++)
		{
			mtxTo.CellRef(x, y) = m_pTo[idx];
			++idx;
		}
	}

	if(bDividNum)
	{
		unsigned size = m_planXDim * m_planYDim;
		DATA invN = 1.F / size;
		for(unsigned y=0; y<m_planYDim; y++)
		{
			for(unsigned x=0; x<m_planXDim; x++)
			{
				mtxTo.CellRef(x, y) *= invN;
			}
		}
	}
	else
	{}
}

/*
void FFT::Test()
{
	cout << "start FFT testing" << endl;
	mtxOp.out.SetWidth(4);
	mtxOp.out.SetPrecision(4);

	cout << "test simple sequence" << endl;
	DATA aVal[] = {1.F, 2.F, 3.F, 4.F, 5.F};
	unsigned inNum = sizeof(aVal) / sizeof(DATA);
	Mtx mtxIn(inNum, 1);
	for(unsigned x=0; x<inNum; x++)
	{
		mtxIn.CellRef(x, 0) = aVal[x];
	}

	CplxMtx cMtx(inNum/2+1, 1);
	mtxOp.fft.Gen(cMtx, mtxIn);
	
	mtxOp.ifft.Gen(mtxIn, cMtx, true);
	for(unsigned x=0; x<inNum; x++)
	{
		assert(myMath.IsEqual(mtxIn.CellVal(x, 0), aVal[x]));
	}	
	cout << "sequence ok" << endl << endl;

	//************************************************************************************************

	cout << "test image convertion" << endl;
	string inName  = TEST_DATA_DIR	+ "\\bmpImg.bmp"; 
	string outName = OUTPUT_DIR		+ "\\fftOut.bmp";
	string tmpName = OUTPUT_DIR		+ "\\tmpOut.bmp";

	MyImg *pImg = imgIO.Read(inName);
	Mtx *pMtxOrg = pImg->GetMtx(0);
	delete pImg;

	Vect2D<unsigned> dim = pMtxOrg->GetDim();
	CplxMtx *pCplxMtx = new CplxMtx(dim.m_x/2+1, dim.m_y);
	mtxOp.fft.Gen(*pCplxMtx, *pMtxOrg);
	delete pMtxOrg;

	Mtx *pMtxNew = new Mtx(dim.m_x, dim.m_y);
	mtxOp.ifft.Gen(*pMtxNew, *pCplxMtx, true);
	delete pCplxMtx;

	MyImg *pImgOut = new MyImg(*pMtxNew);
	delete pMtxNew;
	assert(imgIO.IsBmpFormat(*pImgOut));
	imgIO.Write(outName, *pImgOut);
	delete pImgOut;
	cout << "image conversion ok" << endl << endl;

	cout << "FFT testing ok" << endl << endl;
}
*/

//*************************************************************************************************
//
//*************************************************************************************************

Conv::Conv()
{}

Conv::~Conv()
{}

void Conv::Gen(Mtx &mtxBase, const Mtx &mtxKerl)
{
	Vect2D<unsigned> dimBase = mtxBase.GetDim();
	Vect2D<unsigned> dimKerl = mtxKerl.GetDim();
	Vect2D<unsigned> dimFft(dimBase.m_x + dimKerl.m_x - 1, 
						    dimBase.m_y + dimKerl.m_y - 1);

	cout << "convolve " << dimBase.m_x << "x" << dimBase.m_y << " with "
		 << dimKerl.m_x << "x" << dimKerl.m_y << endl;
	
	unsigned cXDim = dimFft.m_x/2 + 1;
	unsigned cYDim = dimFft.m_y;
	Mtx mBasePad(dimFft.m_x, dimFft.m_y);	mtxOp.zero.Gen(mBasePad);	
	Mtx mKerlPad(dimFft.m_x, dimFft.m_y);	mtxOp.zero.Gen(mKerlPad);	
	Mtx mBaseV(mBasePad, 0, 0, dimBase.m_x, dimBase.m_y);		mBaseV.CopyFrom(mtxBase);
	Mtx mKerlV(mKerlPad, 0, 0, dimKerl.m_x, dimKerl.m_y);		mKerlV.CopyFrom(mtxKerl);

	CplxMtx cmBase(cXDim, cYDim);			m_fft.Gen(cmBase, mBasePad);
	CplxMtx cmKerl(cXDim, cYDim);			m_fft.Gen(cmKerl, mKerlPad);

	CplxMtx cmFft(cXDim, cYDim);
	mtxOp.cellX.Gen(cmFft, cmBase, cmKerl);

	Mtx cFft(dimFft.m_x, dimFft.m_y);
	m_ifft.Gen(cFft, cmFft, true);
	
	unsigned kerlCntX = dimKerl.m_x / 2;
	unsigned kerlCntY = dimKerl.m_y / 2;
	mtxBase.CopyFrom(Mtx(cFft, kerlCntX, kerlCntY, dimBase.m_x, dimBase.m_y));

	cout << "convolution end" << endl << endl;
}

void Conv::Test()
{
	cout << "start Conv testing" << endl;
	mtxOp.out.SetWidth(4);
	mtxOp.out.SetPrecision(4);

	cout << "test simple sequence" << endl;
	DATA aBase[] = {1.F, 2.F, 3.F, 4.F, 5.F, 6.F};
	DATA aKerl[] = {1.F, 1.F, 1.F};
	DATA aReal[] = {3.F, 6.F, 9.F, 12.F, 15.F, 11.F};

	unsigned baseNum = sizeof(aBase) / sizeof(DATA);
	unsigned kerlNum = sizeof(aKerl) / sizeof(DATA);

	Mtx mtxBase(baseNum, 1);
	for(unsigned x=0; x<baseNum; x++)
	{
		mtxBase.CellRef(x, 0) = aBase[x];
	}
	Mtx mtxKerl(kerlNum, 1);
	for(unsigned x=0; x<kerlNum; x++)
	{
		mtxKerl.CellRef(x, 0) = aKerl[x];
	}

	mtxOp.conv.Gen(mtxBase, mtxKerl);

	for(unsigned x=0; x<baseNum; x++)
	{
		assert(myMath.IsEqual(mtxBase.CellVal(x, 0), aReal[x]));
	}
	cout << "sequence ok" << endl << endl;


	cout << "Conv test ok" << endl << endl;
}

//*************************************************************************************************
//
//*************************************************************************************************

Sum::Sum()
{}

Sum::~Sum()
{}

void Sum::Gen(Mtx &mOut, const Mtx &mIn) const
{
	Vect2D<unsigned> dim = mIn.GetDim();
	
	for(unsigned y=0; y<dim.m_y; y++)
	{
		mOut.CellRef(0, y) = mIn.CellVal(0, y);
		for(unsigned x=1; x<dim.m_x; x++)
		{
			mOut.CellRef(x, y) = mOut.CellVal(x-1, y) + mIn.CellVal(x, y);
		}
	}

	for(unsigned x=0; x<dim.m_x; x++)
	{
		for(unsigned y=1; y<dim.m_y; y++)
		{
			mOut.CellRef(x, y) += mOut.CellRef(x, y-1);
		}
	}
}

DATA Sum::Gen(const Mtx &mIn) const
{
	Vect2D<unsigned> dim = mIn.GetDim();

	DATA s = 0;
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			s += mIn.CellVal(x, y);
		}
	}
	return s;
}

//*************************************************************************************************

Avg::Avg()
{}

Avg::~Avg()
{}

DATA Avg::Gen(const Mtx &mIn) const
{
	Vect2D<unsigned> dim = mIn.GetDim();
	unsigned size = dim.m_x * dim.m_y;
	
	DATA s = mtxOp.sum.Gen(mIn);
	DATA a = s / size;
	return a;
}

//*************************************************************************************************

Dev::Dev()
{}

Dev::~Dev()
{}

DATA Dev::Gen(const Mtx &mIn, bool bDivN) const
{
	DATA avg = mtxOp.avg.Gen(mIn);	
	Vect2D<unsigned> dim = mIn.GetDim();
	
	DATA d = 0;
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			DATA diff = mIn.CellVal(x, y) - avg;
			d += diff*diff;
		}
	}

	if(bDivN)
	{
		unsigned num = dim.m_x * dim.m_y;
		d /= num;
	}
	else {}
	d = sqrt(d);
	return d;
}

//*************************************************************************************************

Rng::Rng()
{}

Rng::~Rng()
{}

void Rng::Gen(DATA &vMin, DATA &vMax, const Mtx &m) const
{
	Vect2D<unsigned> dim = m.GetDim();
	vMin = 1e10;
	vMax = -1e10;
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			DATA v = m.CellVal(x, y);
			if(v < vMin)
			{
				vMin = v;
			}
			else {}
			if(v > vMax)
			{
				vMax = v;
			}
			else {}
		}
	}
}

ScaleVal::ScaleVal()
{}

ScaleVal::~ScaleVal()
{}

void ScaleVal::Gen(Mtx &m, DATA min, DATA max) const
{
	DATA mMin = 1e10;
	DATA mMax = -1e10;
	mtxOp.rng.Gen(mMin, mMax, m);
	MyAssert(mMax > mMin);

	DATA scl = (max - min) / (mMax - mMin);
	Vect2D<unsigned> dim = m.GetDim();
	for (unsigned y = 0; y < dim.m_y; y++) {
		for (unsigned x = 0; x < dim.m_x; x++) {
			m.CellRef(x, y) = (m.CellVal(x, y) - mMin) * scl + min;
		}
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

ScaleDim::ScaleDim()
{}

ScaleDim::~ScaleDim()
{}

void ScaleDim::Gen_nearest(Mtx &mOut, const Mtx &mIn) const
{
	Vect2D<unsigned> dimIn  = mIn.GetDim();
	Vect2D<unsigned> dimOut = mOut.GetDim();

	DATA sX = (DATA)(dimIn.m_x-1) / (dimOut.m_x-1);
	DATA sY = (DATA)(dimIn.m_y-1) / (dimOut.m_y-1);

	for(unsigned y=0; y<dimOut.m_y; y++)
	{
		unsigned yInv = (unsigned)(y*sY + 0.5F);
		for(unsigned x=0; x<dimOut.m_x; x++)
		{
			unsigned xInv = (unsigned)(x*sX + 0.5F);
			mOut.CellRef(x, y) = mIn.CellVal(xInv, yInv);
		}
	}
}

void ScaleDim::Gen_bilinear(Mtx &mOut, const Mtx &mIn) const
{
	Vect2D<unsigned> dimIn  = mIn.GetDim();
	Vect2D<unsigned> dimOut = mOut.GetDim();

	DATA sX = (DATA)(dimIn.m_x-1) / (dimOut.m_x-1);
	DATA sY = (DATA)(dimIn.m_y-1) / (dimOut.m_y-1);

	for(unsigned y=0; y<dimOut.m_y; y++)
	{
		DATA yInv = y * sY;
		unsigned yB = (unsigned)yInv;
		unsigned yT = yB + 1; //(unsigned)(yInv + 0.5F);
		if(yT >= dimIn.m_y)
		{
			yT = yB;
		}

		DATA yBRatio = yInv - yB;
		for(unsigned x=0; x<dimOut.m_x; x++)
		{
			DATA xInv = x * sX;
			unsigned xL = (unsigned)xInv;
			unsigned xR = xL + 1; //(unsigned)(xInv + 0.5F);
			if(xR >= dimIn.m_x)
			{
				xR = xL;
			}

			DATA xLRatio = xInv - xL;

			DATA lb = mIn.CellVal(xL, yB);
			DATA rb = mIn.CellVal(xR, yB);
			DATA lt = mIn.CellVal(xL, yT);
			DATA rt = mIn.CellVal(xR, yT);
			mOut.CellRef(x, y) = myMath.Interpolate_linear(lb, rb, lt, rt, xLRatio, yBRatio);
		}
	}
}

void ScaleDim::Gen_blockAvg(Mtx &mOut, const Mtx &mIn, Mtx *pmBuf) const
{
	Vect2D<unsigned> dimIn  = mIn.GetDim();
	Vect2D<unsigned> dimOut = mOut.GetDim();
	if(dimOut.m_x>=dimIn.m_x || dimOut.m_y>dimIn.m_y)
	{
		Gen_bilinear(mOut, mIn);
		return;
	}
	else {}

	mtxOp.sum.Gen(*pmBuf, mIn);

	DATA sX = (DATA)(dimIn.m_x-1) / (dimOut.m_x-1);
	DATA sY = (DATA)(dimIn.m_y-1) / (dimOut.m_y-1);
	
	unsigned xLen = (unsigned)(sX + 0.5F);
	int xL = -((int)xLen / 2); 
	unsigned xR = xL + xLen - 1;

	unsigned yLen = (unsigned)(sY + 0.5F);
	int yB = -((int)yLen / 2);
	unsigned yT = yB + yLen - 1;

	for(unsigned y=0; y<dimOut.m_y; y++)
	{
		unsigned yCnt = (unsigned)(y*sY + 0.5F);
		int yLocB = (int)yCnt + yB;
		if(yLocB < 0)
		{
			yLocB = 0;
		}
		else {}

		unsigned yLocT = yCnt + yT;
		if(yLocT >= dimIn.m_y)
		{
			yLocT = dimIn.m_y - 1;
		}
		else {}

		for(unsigned x=0; x<dimOut.m_x; x++)
		{
			unsigned xCnt = (unsigned)(x*sX + 0.5F);
			int xLocL = xCnt + xL;
			if(xLocL < 0)
			{
				xLocL = 0;
			}
			else {}

			unsigned xLocR = xCnt + xR;
			if(xLocR >= dimIn.m_x)
			{
				xLocR = dimIn.m_x - 1;
			}
			else {}

			unsigned num = (yLocT - yLocB + 1) * (xLocR - xLocL + 1);
			DATA sum = pmBuf->CellVal(xLocR, yLocT) - pmBuf->CellVal(xLocR, yLocB) - pmBuf->CellVal(xLocL, yLocT) + pmBuf->CellVal(xLocL, yLocB);
			mOut.CellRef(x, y) = sum / num;
		}
	}
}

void ScaleDim::Gen(Mtx &mOut, const Mtx &mIn, METHOD m, Mtx *pmBuf) const
{
	switch(m)
	{
	case NEAREST:
		Gen_nearest(mOut, mIn);
		break;

	case BILINEAR:
		Gen_bilinear(mOut, mIn);
		break;

	case BLOCK_AVG:
		Gen_blockAvg(mOut, mIn, pmBuf);
		break;

	default:
		assert(0);
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

Add::Add()
{}

Add::~Add()
{}

void Add::Gen(Mtx &mtxBase, const Mtx &adder) const
{
	Vect2D<unsigned> dim = mtxBase.GetDim();

	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			mtxBase.CellRef(x, y) += adder.CellVal(x, y);
		}
	}
}

//*************************************************************************************************

Sub::Sub()
{}

Sub::~Sub()
{}

void Sub::Gen(Mtx &mtxBase, const Mtx &suber) const
{
	Vect2D<unsigned> dimBase = mtxBase.GetDim();
	Vect2D<unsigned> dimSub  = suber.GetDim();
	MyAssert(dimBase.m_x == dimSub.m_x &&
			 dimBase.m_y == dimSub.m_y);

	for(unsigned y = 0; y < dimBase.m_y; y++) {
		for(unsigned x = 0; x < dimBase.m_x; x++) {
			mtxBase.CellRef(x, y) -= suber.CellVal(x, y);
		}
	}
}

void Sub::Gen(DATA cNum, Mtx &mtxSuber) const
{
	Vect2D<unsigned> dim = mtxSuber.GetDim();

	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			mtxSuber.CellRef(x, y) = cNum - mtxSuber.CellVal(x, y);
		}
	}
}

//*************************************************************************************************

Mul::Mul()
{}

Mul::~Mul()
{}

void Mul::Gen(Mtx &mtxR, const Mtx &mtxA, const Mtx &mtxB) const
{
	Vect2D<unsigned> dimR = mtxR.GetDim();
	Vect2D<unsigned> dimA = mtxA.GetDim();
	Vect2D<unsigned> dimB = mtxB.GetDim();
	MyAssert(dimR.m_x == dimB.m_x && 
			 dimR.m_y == dimA.m_y && 
			 dimA.m_x == dimB.m_y);

	for (unsigned y = 0; y < dimR.m_y; y++) {
		for (unsigned x = 0; x < dimR.m_x; x++) {
			DATA r = 0;
			for (unsigned c = 0; c < dimA.m_x; c++) {
				r += mtxA.CellVal(c, y) * mtxB.CellVal(x, c);
			} // c
			mtxR.CellRef(x, y) = r;
		} // x
	} // y
}

void Mul::Gen(Mtx &mtx, DATA scl) const
{
	Vect2D<unsigned> dim = mtx.GetDim();
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			mtx.CellRef(x, y) *= scl;
		}
	}
}

void Mul::GenTA(Mtx &mtxR, const Mtx &mtxA) const
{
	Vect2D<unsigned> dimA = mtxA.GetDim();
	Vect2D<unsigned> dimR = mtxR.GetDim();
	MyAssert(dimR.m_x == dimR.m_y);
	MyAssert(dimR.m_y == dimA.m_x);

	for (unsigned y = 0; y < dimR.m_y; y++) {
		for (unsigned x = 0; x < dimR.m_x; x++) {
			DATA r = 0;
			for (unsigned c = 0; c < dimA.m_y; c++) {
				r += mtxA.CellVal(y, c) * mtxA.CellVal(x, c);
			} // c
			mtxR.CellRef(x, y) = r;
		} // x
	} // y
}

void Mul::GenTA(Mtx &mtxR, const Mtx &mtxA, const Mtx &mtxB) const
{
	Vect2D<unsigned> dimR = mtxR.GetDim();
	Vect2D<unsigned> dimA = mtxA.GetDim();
	Vect2D<unsigned> dimB = mtxB.GetDim();
	MyAssert(dimR.m_x == dimB.m_x && 
			 dimR.m_y == dimA.m_x && 
			 dimA.m_y == dimB.m_y);

	for(unsigned y = 0; y < dimR.m_y; y++) {
		for(unsigned x = 0; x < dimR.m_x; x++) {
			DATA r = 0;
			for(unsigned c = 0; c < dimA.m_y; c++) {
				r += mtxA.CellVal(y, c) * mtxB.CellVal(x, c);
			} // c
			mtxR.CellRef(x, y) = r;
		} // x
	} // y
}

void Mul::GenTB(Mtx &mtxR, const Mtx &mtxA, const Mtx &mtxB) const
{
	Vect2D<unsigned> dimR = mtxR.GetDim();
	Vect2D<unsigned> dimA = mtxA.GetDim();
	Vect2D<unsigned> dimB = mtxB.GetDim();
	MyAssert(dimR.m_x == dimB.m_y && 
			 dimR.m_y == dimA.m_y && 
			 dimA.m_x == dimB.m_x);

	for(unsigned y = 0; y < dimR.m_y; y++) {
		for(unsigned x = 0; x < dimR.m_x; x++) {
			DATA r = 0;
			for(unsigned c = 0; c < dimA.m_x; c++) {
				r += mtxA.CellVal(c, y) * mtxB.CellVal(c, x);
			} // c
			mtxR.CellRef(x, y) = r;
		} // x
	} // y
}

void Mul::Test()
{
	cout << "start Mul testing" << endl;

	const unsigned ax = 2;	const unsigned ay = 3;
	const unsigned bx = 3;  const unsigned by = ax;
	DATA aaA[ay][ax] = {{1.F, 2.F}, {3.F, 4.F}, {5.F, 6.F}};
	DATA aaB[by][bx] = {{7.F, 8.F, 9.F}, {10.F, 11.F, 12.F}};
	DATA aaR[ay][bx] = {{27.F, 30.F, 33.F}, {61.F, 68.F, 75.F}, {95.F, 106.F, 117.F}};

	Mtx mA(ax, ay);
	for(unsigned y=0; y<ay; y++)
	{
		for(unsigned x=0; x<ax; x++)
		{
			mA.CellRef(x, y) = aaA[y][x];
		}
	}

	Mtx mB(bx, by);
	for(unsigned y=0; y<by; y++)
	{
		for(unsigned x=0; x<bx; x++)
		{
			mB.CellRef(x, y) = aaB[y][x];
		}
	}

	Mtx mR(bx, ay);
	mtxOp.mul.Gen(mR, mA, mB);
	for(unsigned y=0; y<ay; y++)
	{
		for(unsigned x=0; x<bx; x++)
		{
			assert(myMath.IsEqual(mR.CellVal(x, y), aaR[y][x]));
		}
	}

	cout << "Mul test ok" << endl;
}

//*************************************************************************************************
//
//*************************************************************************************************

Hessian::Hessian()
{}

Hessian::~Hessian()
{}

void Hessian::Gen(Mtx2D &mtxH, const Mtx &mtx, unsigned x, unsigned y) const
{	
	Vect2D<unsigned> dim = mtx.GetDim();

	DATA xxDer, yyDer, xyDer;
	if(x>=1 && x<=dim.m_x-2 &&
	   y>=1 && y<=dim.m_y-2)
	{
		unsigned xR = x + 1;	unsigned xL = x - 1;
		unsigned yT = y + 1;	unsigned yB = y - 1;

		xxDer = mtx.CellVal(xL,  y) + mtx.CellVal(xR,  y) - 2 * mtx.CellVal(x, y);
		yyDer = mtx.CellVal( x, yB) + mtx.CellVal( x, yT) - 2 * mtx.CellVal(x, y);

		xyDer = (mtx.CellVal(xR, yT) - mtx.CellVal(xR, yB) - mtx.CellVal(xL, yT) + mtx.CellVal(xL, yB)) * 0.25F;
	}
	else
	{
		Mtx subMtx(3, 3);
		mtxOp.zero.Gen(subMtx);

		for(unsigned yy=0; yy<3; yy++)
		{
			int yLoc = yy - 1 + y;
			if(yLoc<0 || yLoc>=(int)dim.m_y)	continue;
			else {}

			for(unsigned xx=0; xx<3; xx++)
			{
				int xLoc = xx - 1 + x;
				if(xLoc<0 || xLoc>=(int)dim.m_x)	continue;
				else {}

				subMtx.CellRef(xx, yy) = mtx.CellVal(xLoc, yLoc);
			}
		}

		xxDer = subMtx.CellVal(0, 1) + subMtx.CellVal(2, 1) - 2 * subMtx.CellVal(1, 1);
		yyDer = subMtx.CellVal(1, 0) + subMtx.CellVal(1, 2) - 2 * subMtx.CellVal(1, 1);

		xyDer = (subMtx.CellVal(2, 2) - subMtx.CellVal(2, 0) - subMtx.CellVal(0, 2) + subMtx.CellVal(0, 0)) * 0.25F;
	}

	mtxH.CellRef(0, 0) = xxDer;		mtxH.CellRef(1, 0) = xyDer;	
	mtxH.CellRef(0, 1) = xyDer;		mtxH.CellRef(1, 1) = yyDer;	
}

//*************************************************************************************************
//
//*************************************************************************************************

Derivate::Derivate()
{}

Derivate::~Derivate()
{}

void Derivate::Gen(Mtx &mtxDer, const Mtx &mtx, unsigned x, unsigned y, bool bNormal) const
{	
	Vect2D<unsigned> dim = mtx.GetDim();

	if(x>=1 && x<=dim.m_x-2 &&
	   y>=1 && y<=dim.m_y-2)
	{
		unsigned xR = x + 1;	unsigned xL = x - 1;
		unsigned yT = y + 1;	unsigned yB = y - 1;

		mtxDer.CellRef(0, 0) = mtx.CellVal(xR,  y) - mtx.CellVal(xL,  y);
		mtxDer.CellRef(0, 1) = mtx.CellVal( x, yT) - mtx.CellVal( x, yB);
	}
	else
	{
		Mtx subMtx(3, 3);
		for(unsigned yy=0; yy<3; yy++)
		{
			for(unsigned xx=0; xx<3; xx++)
			{
				subMtx.CellRef(xx, yy) = mtx.CellVal(x, y);
			}
		}
		for(unsigned yy=0; yy<3; yy++)
		{
			int yLoc = yy - 1 + y;
			if(yLoc<0 || yLoc>=(int)dim.m_y)	continue;
			else {}

			for(unsigned xx=0; xx<3; xx++)
			{
				int xLoc = xx - 1 + x;
				if(xLoc<0 || xLoc>=(int)dim.m_x)	continue;
				else {}

				subMtx.CellRef(xx, yy) = mtx.CellVal(xLoc, yLoc);
			}
		}

		mtxDer.CellRef(0, 0) = subMtx.CellVal(2, 1) - subMtx.CellVal(0, 1);
		mtxDer.CellRef(0, 1) = subMtx.CellVal(1, 2) - subMtx.CellVal(1, 0);
	}

	if(bNormal)
	{
		mtxDer.CellRef(0, 0) *= 0.5F;
		mtxDer.CellRef(0, 1) *= 0.5F;
	}
	else {}
}

//*************************************************************************************************
//
//*************************************************************************************************

Median::Median()
{}

Median::~Median()
{}

int Compare(const void *pA, const void *pB)
{
	if(myMath.IsLess(*(DATA*)pA, *(DATA*)pB))
	{
		return -1;
	}
	else if(myMath.IsGreat(*(DATA*)pA, *(DATA*)pB))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

DATA Median::Gen(DATA *pData, unsigned size)
{
	qsort(pData, size, sizeof(DATA), Compare); 
	return pData[size/2];
}

MedFlt::MedFlt()
{}

MedFlt::~MedFlt()
{}

void MedFlt::Gen(Mtx &mtxOut, Mtx &mtxIn, unsigned winLen)
{
	Vect2D<unsigned> dimIn = mtxIn.GetDim();

	DATA *pNei = new DATA[winLen*winLen];
	for(unsigned y=0; y<dimIn.m_y; y++)
	{
		for(unsigned x=0; x<dimIn.m_x; x++)
		{
			unsigned neiNum = 0;
			for(unsigned yy=0; yy<winLen; yy++)
			{
				int yLoc = (int)y + yy - winLen/2;
				if(yLoc<0 || yLoc>=(int)dimIn.m_y)
				{
					continue;
				}
				else {}

				for(unsigned xx=0; xx<winLen; xx++)
				{
					int xLoc = (int)x + xx - winLen/2;
					if(xLoc<0 || xLoc>=(int)dimIn.m_x)
					{
						continue;
					}
					else {}

					pNei[neiNum++] = mtxIn.CellVal(xLoc, yLoc);
				} // xx
			} // yy

			mtxOut.CellRef(x, y) = mtxOp.median.Gen(pNei, neiNum);
		} // x
	} // y
	delete []pNei;
}

//*************************************************************************************************
//
//*************************************************************************************************

His::His()
{
	m_maxV = 1.F;
	m_minV = 0;
	m_vLen = m_maxV - m_minV;
	m_invVLen = 1.F / m_vLen;
}

His::~His()
{}

//***********************************************

void His::SetBound(DATA maxV, DATA minV)
{
	m_maxV = maxV;
	m_minV = minV;
	m_vLen = m_maxV - m_minV;
	m_invVLen = 1.F / m_vLen;
}

void His::Add(Mtx &mtxHis, DATA val) const
{
	Vect2D<unsigned> dim = mtxHis.GetDim();
	unsigned hisMax = dim.m_y - 1;

	int loc = (int)((val - m_minV) * m_invVLen * hisMax);
	if(loc < 0)
	{
		loc = 0;
	}
	else if(loc > (int)hisMax)
	{
		loc = hisMax;
	}
	else {}
	mtxHis.CellRef(0, loc) = mtxHis.CellVal(0, loc) + 1.F;
} 

//***********************************************

void His::Gen(Mtx &mtxHis, Mtx &mtxIn, DATA maxV, DATA minV)
{
	Vect2D<unsigned> dimHis = mtxHis.GetDim();
	assert(dimHis.m_x == 1);
	for(unsigned y=0; y<dimHis.m_y; y++)
	{
		mtxHis.CellRef(0, y) = 0;
	}

	assert(maxV > minV);
	DATA vLen = maxV - minV;
	DATA vLenRatio = 1.F / vLen;
	Vect2D<unsigned> dimIn = mtxIn.GetDim();
	for(unsigned y=0; y<dimIn.m_y; y++)
	{
		for(unsigned x=0; x<dimIn.m_x; x++)
		{
			int hisLoc = (int)((mtxIn.CellVal(x, y) - minV) * vLenRatio * (dimHis.m_y - 1));
			if(hisLoc<0 || hisLoc>=(int)dimHis.m_y)
			{
				continue;
			}
			else {}
			/*
			if(hisLoc < 0)
			{
				hisLoc = 0;
			}
			else if(hisLoc >= (int)dimHis.m_y)
			{
				hisLoc = dimHis.m_y - 1;
			}
			else {}
			*/
			mtxHis.CellRef(0, hisLoc) = mtxHis.CellVal(0, hisLoc) + 1.F;
		}
	}
}

//*************************************************************************************************

HisEqual::HisEqual()
	:m_levelNum(0), m_memSize(m_levelNum)
{
	m_pCount = 0;  //new unsigned[m_memSize];
}

HisEqual::~HisEqual()
{
	delete []m_pCount;
}

void HisEqual::Reval(Mtx &mtx, Mtx *pMSkip)
{
	unsigned sNum = m_pCount[m_levelNum-1];

	unsigned cMin = 0;
	for(unsigned i=0; i<m_levelNum; i++) {
		if(m_pCount[i] > 0) {
			cMin = m_pCount[i];
			break;
		} else {}
	}

	unsigned mtxSize = sNum; 
	DATA d_count = mtxSize - cMin;
	Vect2D<unsigned> dim = mtx.GetDim();
	for(unsigned y=0; y<dim.m_y; y++) {
		for(unsigned x=0; x<dim.m_x; x++) {
			if(pMSkip && pMSkip->CellVal(x, y) == 0) {
				continue;
			} else {}

			DATA val = mtx.CellVal(x, y);
			unsigned level = (unsigned)(val * (m_levelNum-1));

			if(m_pCount[level] == 0) {
				mtx.CellRef(x, y) = 0;
			} else {
				mtx.CellRef(x, y) = (m_pCount[level] - cMin) / d_count * (m_levelNum - 1);
			}
		}
	}
}

void HisEqual::Gen(Mtx &mtx, unsigned levelNum, Mtx *pMSkip)
{
	if(levelNum != m_levelNum) {
		if(m_memSize < levelNum) {
			delete []m_pCount;
			m_memSize = levelNum;
			m_pCount = new unsigned[m_memSize];
		} else {}
		m_levelNum = levelNum;
	} else {}
	for(unsigned i = 0; i < m_levelNum; i++) {
		m_pCount[i] = 0;
	}

	Vect2D<unsigned> dim = mtx.GetDim();

	DATA min = 1e10;
	DATA max = -1e10;
	mtxOp.rng.Gen(min, max, mtx);
	MyAssert(max > min);
	DATA d_m = max - min;
	for (unsigned y = 0; y < dim.m_y; y++) {
		for (unsigned x = 0; x < dim.m_x; x++) {
			mtx.CellRef(x, y) = (mtx.CellVal(x, y) - min) / d_m;
		} // x
	} // y

	for(unsigned y = 0; y < dim.m_y; y++) {
		for(unsigned x = 0; x < dim.m_x; x++) {
			if(pMSkip && pMSkip->CellVal(x, y) == 0) {
				continue;
			} else {}

			unsigned level = (unsigned)(mtx.CellVal(x, y) * (m_levelNum-1));
			m_pCount[level]++; 
		}
	}

	for(unsigned i = 1; i < m_levelNum; i++) {
		m_pCount[i] += m_pCount[i-1];
	}

	/*
	for (unsigned i = 0; i < m_levelNum; i++) {
		cout << i << ":" << m_pCount[i] << " ";
		PrintLine(i, 5);
	}
	cout << endl;
	*/

	Reval(mtx , pMSkip);
}

void HisEqual::Gen(unsigned char *pIn, unsigned size, unsigned char *pSkip)
{
	unsigned levelNum = 256;
	if(levelNum != m_levelNum) {
		if(m_memSize < levelNum) {
			delete []m_pCount;
			m_memSize = levelNum;
			m_pCount = new unsigned[m_memSize];
		} else {}
		m_levelNum = levelNum;
	} else {}
	for(unsigned i = 0; i < m_levelNum; i++) {
		m_pCount[i] = 0;
	}

	for (unsigned s = 0; s < size; s++) {
		if (pSkip && pSkip[s] == 0) {
			continue;
		} else {}

		m_pCount[pIn[s]]++;
	}

	for(unsigned i = 1; i < m_levelNum; i++) {
		m_pCount[i] += m_pCount[i - 1];
	}

	//*******************************************
	// reval
	//*******************************************

	unsigned sNum = m_pCount[m_levelNum - 1];

	unsigned cMin = 0;
	for(unsigned i=0; i<m_levelNum; i++) {
		if(m_pCount[i] > 0) {
			cMin = m_pCount[i];
			break;
		} else {}
	}

	unsigned inSize = sNum; 
	unsigned d_count = inSize - cMin;
	for(unsigned s = 0; s < size; s++) {
		if(pSkip && pSkip[s] == 0) {
			continue;
		} else {}

		unsigned level = pIn[s];

		if(m_pCount[level] == 0) {
			pIn[s] = 0;
		} else {
			pIn[s] = (m_levelNum - 1) * (m_pCount[level] - cMin) / d_count;
		}
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

Distance::Distance()
{}

Distance::~Distance()
{}

DATA Distance::Gen(Mtx &mtxA, Mtx &mtxB, METHOD m)
{
	Vect2D<unsigned> dimA = mtxA.GetDim();
	Vect2D<unsigned> dimB = mtxB.GetDim();
	assert(dimA.m_x==dimB.m_x && dimA.m_y==dimB.m_y);

	if(m == SSD)
	{
		DATA difSum = 0;
		for(unsigned y=0; y<dimA.m_y; y++)
		{
			for(unsigned x=0; x<dimA.m_x; x++)
			{
				DATA diff = mtxA.CellVal(x, y) - mtxB.CellVal(x, y);
				difSum += diff*diff;
			}
		}
		return difSum;
	}
	else if(m == SAD)
	{
		DATA difSum = 0;
		for(unsigned y=0; y<dimA.m_y; y++)
		{
			for(unsigned x=0; x<dimA.m_x; x++)
			{
				DATA diff = mtxA.CellVal(x, y) - mtxB.CellVal(x, y);
				difSum += fabs(diff);
			}
		}
		return difSum;
	}
	else
	{
		assert(0);
		return 0;
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

BinThrd::BinThrd()
{}

BinThrd::~BinThrd()
{}

void BinThrd::Gen(Mtx &mtx, DATA thrd, bool bRmLow)
{
	Vect2D<unsigned> dim = mtx.GetDim();
	for(unsigned y = 0; y < dim.m_y; y++) {
		for(unsigned x = 0; x < dim.m_x; x++) {
			if (bRmLow) {				
				mtx.CellRef(x, y) = (mtx.CellVal(x, y) < thrd) ?
					0 : 1.F;
			} else {
				mtx.CellRef(x, y) = (mtx.CellVal(x, y) > thrd) ?
					0 : 1.F;
			} // bRmLow
		} // x
	} // y
}

Clamp::Clamp()
{}

Clamp::~Clamp()
{}

void Clamp::Gen(Mtx &mtx, DATA thrd, bool bRmLow, DATA repV)
{
	Vect2D<unsigned> dim = mtx.GetDim();
	for(unsigned y = 0; y < dim.m_y; y++) {
		for(unsigned x = 0; x < dim.m_x; x++) {
			if (bRmLow) {				
				if (mtx.CellVal(x, y) < thrd) {
					mtx.CellRef(x, y) = repV;
				} else {}
			} else {
				if (mtx.CellVal(x, y) > thrd) {
					mtx.CellRef(x, y) = repV;
				} else {}
			} // bRmLow
		} // x
	} // y
}

//*************************************************************************************************

Sobel::Sobel() 
{}

Sobel::~Sobel()
{}

void Sobel::Gen(Layer &lyrOut, const Mtx &mtxIn, bool bMag, bool bDebug)
{
	if (bDebug) {
		cout << "Sobel" << endl;
	} else {}

	Vect2D<unsigned> dimIn = mtxIn.GetDim();
	Vect3D<unsigned> dimOut = lyrOut.GetDim();
	MyAssert(dimIn.m_x == dimOut.m_x &&
			 dimIn.m_y == dimOut.m_y);
	if (bMag) {
		MyAssert(dimOut.m_z >= 3);
	} else {
		MyAssert(dimOut.m_z >= 2);
	}

	// Gx = [-1 0 +1]
	//      [-2 0 +2]
	//      [-1 0 +1]
	//
	// Gy = [+1 +2 +1]
	//      [ 0  0  0]
	//      [-1 -2 -1]
	static Mtx mtxBound(3, 3);
	for (unsigned y = 0; y < dimOut.m_y; y++) {
		for (unsigned x = 0; x < dimOut.m_x; x++) {
			if (bDebug) {
				cout << "(" << x << "," << y << ") "; 
			} else {}

			Mtx *pMtxR = 0;
			if (x != 0 && x != dimOut.m_x - 1 &&
				y != 0 && y != dimOut.m_y - 1) {
				pMtxR = new Mtx(mtxIn, x - 1, y - 1, 3, 3);
			} else {
				for (unsigned yy = 0; yy < 3; yy++) {
					int yLoc = (int)y + yy - 1;
					if (yLoc < 0) {
						yLoc = 0;
					} else if (yLoc >= (int)dimOut.m_y) {
						yLoc = dimOut.m_y - 1;
					} else {}

					for (unsigned xx = 0; xx < 3; xx++) {
						int xLoc = (int)x + xx - 1;
						if (xLoc < 0) {
							xLoc = 0;
						} else if (xLoc >= (int)dimOut.m_x) {
							xLoc = dimOut.m_x - 1;
						} else {}

						mtxBound.CellRef(xx, yy) = mtxIn.CellVal(xLoc, yLoc);
					} // xx
				} // yy
				pMtxR = new Mtx(mtxBound);
				if (bDebug) {
					cout << "b ";
				} else {}
			} // if bound

			MyAssert(pMtxR != 0);
			DATA xD = 
				- pMtxR->CellVal(0, 0) - 2.F * pMtxR->CellVal(0, 1) - pMtxR->CellVal(0, 2)
				+ pMtxR->CellVal(2, 0) + 2.F * pMtxR->CellVal(2, 1) + pMtxR->CellVal(2, 2);

			DATA yD = 
				- pMtxR->CellVal(0, 0) - 2.F * pMtxR->CellVal(1, 0) - pMtxR->CellVal(2, 0)
				+ pMtxR->CellVal(0, 2) + 2.F * pMtxR->CellVal(1, 2) + pMtxR->CellVal(2, 2);

			xD /= 4.F;
			yD /= 4.F;
			lyrOut.CellRef(x, y, 0) = xD;
			lyrOut.CellRef(x, y, 1) = yD;

			delete pMtxR;

			if (bMag) {
				//MyAssert(0);
				lyrOut.CellRef(x, y, 2) = sqrt(xD * xD + yD * yD);
			} else {}
		} // x
	} // y

	if (bDebug) {
		cout << "Sobel ok" << endl;
	} else {}
}

SIFT_desc::SIFT_desc()
	: m_patchLen(m_blockLen * m_blockNum)
	, m_superPL((unsigned)(m_patchLen * 1.42F + 0.5F) + 1)
	//, m_mtxG(m_patchLen, m_patchLen)
	, m_mtxG(m_superPL, m_superPL)
{
	DATA delta = 1.8F; //1.5F;
	mtxOp.Gauss.Gen(m_mtxG, delta, delta, true);
}

SIFT_desc::~SIFT_desc()
{}

unsigned SIFT_desc::GetDscNum()
{
	return m_blockNum * m_blockNum * m_subDNum;
}

unsigned SIFT_desc::GetSuperPL()
{
	return m_superPL;
}

void AccountBlock(DATA aBin[], unsigned binNum, Mtx &mtxAng, Mtx &mtxMag, bool bDebug)
{
	if (bDebug) {
		cout << "AccountBlock" << endl;
	} else {}

	Vect2D<unsigned> dimAng = mtxAng.GetDim();
	Vect2D<unsigned> dimMag = mtxMag.GetDim();
	MyAssert(dimAng.m_x == dimMag.m_x &&
			 dimAng.m_y == dimMag.m_y);

	unsigned subAng = 360 / binNum;
	MyAssert(360 % binNum == 0);

	for (unsigned i = 0; i < binNum; i++) {
		aBin[i] = 0;
	}

	for (unsigned y = 0; y < dimAng.m_y; y++) {
		for (unsigned x = 0; x < dimAng.m_x; x++) {
			MyAssert(mtxAng.CellVal(x, y) >= 0);
			unsigned idx = (unsigned)mtxAng.CellVal(x, y) / subAng;
			idx = idx % (binNum);
			if (bDebug) {
				if (idx >= binNum) {
					cout << "idx: " << idx << " "
				 	 	 << "ang: " << mtxAng.CellVal(x, y) << endl;
				} else {}
			} else {}
			MyAssert(idx < binNum);
			aBin[idx] += mtxMag.CellVal(x, y);
		}
	}

	if (bDebug) {
		cout << "AccountBlock ok" << endl;
	} else {}
}

//void SIFT_desc::GetOrient(vector<unsigned> &ort, vector<DATA> &mag, Mtx &mtxAng, Mtx &mtxMag, bool bDebug)
void SIFT_desc::GetOrient(unsigned &ort, DATA &mag, Mtx &mtxAng, Mtx &mtxMag, bool bDebug)
{
	static DATA aOrtBin[m_ortNum];
	AccountBlock(aOrtBin, m_ortNum, mtxAng, mtxMag, bDebug);

	DATA angFst = 0;	unsigned afIdx = 0;
	DATA angSnd = 0;	unsigned asIdx = 0;
	for (unsigned i = 0; i < 36; i++) {
		DATA ac = aOrtBin[i];
		if (ac > angFst) {
			angSnd = angFst;	asIdx = afIdx;
			angFst = ac;		afIdx = i;
		} else if (ac > angSnd) {
			angSnd = ac;		asIdx = i;
		} else {}
	}	

	ort = afIdx * 360 / m_ortNum;	mag = angFst;
	/*
	ort[0] = afIdx * 360 / m_ortNum;	mag[0] = angFst;
	ort[1] = 0;							mag[1] = 0;
	if (angSnd > 0.8F * angFst) {
		ort[1] = asIdx * 360 / m_ortNum;
		mag[1] = angSnd;
	} else {}
	*/
}

void SIFT_desc::Gen(unsigned &rOrt, DATA &rMag, DATA dsc[], 
		Mtx &mtxIn, Layer &lyrGrd, Layer &lyrRot, bool bDebug)
{	
	if (bDebug) {
		cout << "sift_desc" << endl;
	} else {}

	Vect2D<unsigned> dimIn = mtxIn.GetDim();
	Vect3D<unsigned> dimGrd = lyrGrd.GetDim();
	Vect3D<unsigned> dimRot = lyrRot.GetDim();
	//MyAssert(dimIn.m_x == m_patchLen &&
	//		 dimIn.m_y == m_patchLen);
	MyAssert(dimIn.m_x == m_superPL &&
			 dimIn.m_y == m_superPL);
	MyAssert(dimRot.m_x == m_patchLen &&
			 dimRot.m_y == m_patchLen);
	MyAssert(dimIn.m_x == dimGrd.m_x &&
			 dimIn.m_y == dimGrd.m_y);
	MyAssert(dimGrd.m_z >= 3);
	MyAssert(dimRot.m_z >= 2);

	mtxOp.Sobel.Gen(lyrGrd, mtxIn, true);

	Mtx *pMtxMag = lyrGrd.GetMtx(2);
	mtxOp.cellX.Gen(*pMtxMag, m_mtxG);

	Mtx &mtxAng = mtxIn;
	for (unsigned y = 0; y < dimIn.m_y; y++) {
		for (unsigned x = 0; x < dimIn.m_x; x++) {
			DATA ang = atan2(lyrGrd.CellVal(x, y, 1), lyrGrd.CellVal(x, y, 0));
			ang = (ang + PI) * R2D;
			mtxAng.CellRef(x, y) = ang;
		} // x
	} // y

	//*******************************************
	// compute orientation of this patch
	//*******************************************
	unsigned ort;
	DATA mag;	
	GetOrient(ort, mag, mtxAng, *pMtxMag, false); //true);
	if (bDebug) {
		//cout << "orient: " << ort[0] << " " << ort[1] << endl;
		cout << "orient: " << ort << endl;
	} else {}
	rOrt = ort;
	rMag = mag;

	//*******************************************
	// rotate
	//*******************************************	
	if (bDebug) {
		cout << "rotate" << endl;
	} else {}
	for (unsigned y = 0; y < dimIn.m_y; y++) {
		for (unsigned x = 0; x < dimIn.m_x; x ++) {
			mtxAng.CellRef(x, y) -= rOrt;

			if (mtxAng.CellVal(x, y)  < 0) {
				mtxAng.CellRef(x, y) += 360.F;
			} else {}
		}
	}

	// x = cosA, 	y = sinA
	// xx = cosB,	yy = sinB
	// sin(A + B) = sinA * cosB + cosA * sinB = y * xx + x * yy
	// cos(A + B) = cosA * cosB - sinA * sinB = x * xx - y * yy
	DATA rA = rOrt * D2R;
	DATA xx = cos(rA);
	DATA yy = sin(rA);

	Mtx mtxAng2(*lyrRot.GetMtx(0));
	Mtx mtxMag2(*lyrRot.GetMtx(1));
	for (unsigned y = 0; y < dimRot.m_y; y++) {
		int yDis = (int)y - dimRot.m_y / 2;
		for (unsigned x = 0; x < dimRot.m_x; x++) {
			int xDis = (int)x - dimRot.m_x / 2;

			int xLoc = (int)(xDis * xx - yDis * yy + dimIn.m_x / 2.F + 0.5F);
			int yLoc = (int)(yDis * xx + xDis * yy + dimIn.m_y / 2.F + 0.5F);
			MyAssert(mtxAng.IsInside(xLoc, yLoc));

			mtxAng2.CellRef(x, y) = mtxAng.CellVal(xLoc, yLoc);
			mtxMag2.CellRef(x, y) = pMtxMag->CellVal(xLoc, yLoc);
		}
	}
	if (bDebug) {
		cout << "rotate ok" << endl;
	} else {}

	//*******************************************
	// account descriptor
	//*******************************************
	if (bDebug) {
		cout << "block no:" << endl;
	} else {}
	for (unsigned yB = 0; yB < m_blockNum; yB++) {
		for (unsigned xB = 0; xB < m_blockNum; xB++) {
			Mtx mtxAng_b(mtxAng2,  xB * m_blockLen, yB * m_blockLen, m_blockLen, m_blockLen);
			Mtx mtxMag_b(mtxMag2,  xB * m_blockLen, yB * m_blockLen, m_blockLen, m_blockLen);

			unsigned dscOff = (yB * m_blockNum + xB) * m_subDNum;
			AccountBlock(&dsc[dscOff], m_subDNum, mtxAng_b, mtxMag_b, bDebug);
		} // xB
	} // yB
	
	
	//*******************************************
	// normalization
	//*******************************************
	unsigned dscNum = m_patchLen * m_subDNum;
	MyAssert(dscNum == 128);
	DATA dscSum = 0;
	for (unsigned d = 0; d < dscNum; d++) {
		dscSum += dsc[d] * dsc[d];
	}

	if (dscSum < 1e-8) {
		DATA uu = 1.F / dscNum;
		uu = sqrt(uu);
		for (unsigned d = 0; d < dscNum; d++) {
			dsc[d] = uu;
		}
	} else {
		dscSum = sqrt(dscSum);

		bool bRe = false;
		DATA ds2 = 0;
		for (unsigned d = 0; d < dscNum; d++) {
			dsc[d] /= dscSum;
			if (dsc[d] > 0.2F) {
				dsc[d] = 0.2F;
				bRe = true;
			} else {}
			ds2 += dsc[d] * dsc[d];
		}

		if (bRe) {
			ds2 = sqrt(ds2);
			for (unsigned d = 0; d < dscNum; d++) {
				dsc[d] /= ds2;
			}
		} else {}
	}

	if (bDebug) {
		cout << "sift_desc ok" << endl;
	} else {}
	
}

Otsu::Otsu()
{}

Otsu::~Otsu()
{}

DATA Otsu::Gen(Mtx &mtx, DATA &thrd, unsigned lev, DATA valMax, DATA valMin)
{
	MyAssert(0);
	return 0;
	/*
	Vect2D<unsigned> dimIn = mtx.GetDim();
	//unsigned size = dimIn.m_x * dimIn.m_y;
	unsigned size = 0;
	for(unsigned y=0; y<dimIn.m_y; y++)
	{
		for(unsigned x=0; x<dimIn.m_x; x++)
		{
			DATA val = mtx.CellVal(x, y);
			if(val>=valMin && val<=valMax)
			{
				size++;
			}
		}
	}
	DATA invSize = 1.F / size;

	Vect2D<unsigned> dimHis(1, lev);
	Mtx mtxHis(dimHis);
	mtxOp.his.Gen(mtxHis, mtx, valMax, valMin);
	mtxOp.mul.Gen(mtxHis, invSize);

	Mtx mtxHVal(dimHis);
	DATA valLen = valMax - valMin;
	DATA valStep = valLen / (lev - 1);
	for(unsigned i=0; i<dimHis.m_y; i++)
	{
		mtxHVal.CellRef(0, i) = valStep * i + valMin;
		//cout << mtxHVal.CellVal(0, i) << " ";
	}

	Mtx mtxP(dimHis);
	mtxP.CellRef(0, 0) = mtxHis.CellVal(0, 0);
	for(unsigned i=1; i<dimHis.m_y; i++)
	{
		mtxP.CellRef(0, i) = mtxP.CellVal(0, i-1) + mtxHis.CellVal(0, i);
		//cout << mtxP.CellVal(0, i) << " ";
	}

	Mtx mtxM(dimHis);
	mtxM.CellRef(0, 0) = mtxHis.CellVal(0, 0) * mtxHVal.CellVal(0, 0);
	for(unsigned i=1; i<dimHis.m_y; i++)
	{
		mtxM.CellRef(0, i) = mtxM.CellVal(0, i-1) + mtxHis.CellVal(0, i) * mtxHVal.CellVal(0, i);
		//cout << mtxM.CellVal(0, i) << " ";
	}

	DATA mg = mtxM.CellVal(0, dimHis.m_y-1);
	DATA dg = 0;
	for(unsigned i=0; i<dimHis.m_y; i++)
	{
		DATA diff = mtxHVal.CellVal(0, i) - mg;
		dg += diff * diff * mtxHis.CellVal(0, i);
	}

	Mtx mtxDelta(dimHis);
	for(unsigned i=0; i<dimHis.m_y; i++)
	{
		DATA diff = mg*mtxP.CellVal(0, i) - mtxM.CellVal(0, i);
		DATA div = mtxP.CellVal(0, i) * (1.F - mtxP.CellVal(0, i));

		if(div <= 0)
		{
			mtxDelta.CellRef(0, i) = 0;
		}
		else 
		{
			mtxDelta.CellRef(0, i) = diff * diff / div;
		}
		//cout << mtxDelta.CellVal(0, i) << " ";
	}

	unsigned num = 1;
	unsigned sumI = 0; //mtxHVal.CellVal(0, 0);
	DATA maxD = mtxDelta.CellVal(0, 0);
	for(unsigned i=1; i<dimHis.m_y; i++)
	{
		if(mtxDelta.CellVal(0, i) == maxD)
		{
			num++;
			sumI += i; //mtxHVal.CellVal(0, i);
		}
		else if(mtxDelta.CellVal(0, i) > maxD)
		{
			num = 1;
			sumI = i; //mtxHVal.CellVal(0, i);
			maxD = mtxDelta.CellVal(0, i);
		}
		else {}
	}
	unsigned iLoc = (unsigned)((DATA)sumI / num + 0.5F);
	thrd = mtxHVal.CellVal(0, iLoc);

	mtxOp.thrd.Gen(mtx, thrd, true);

	DATA conf = mtxDelta.CellVal(0, iLoc) / dg;
	return conf;
	*/
}

//*************************************************************************************************
//
//*************************************************************************************************

Dilate::Dilate()
{}

Dilate::~Dilate()
{}

void Dilate::Gen(Mtx &mtxOut, Mtx &mtxIn, unsigned len)
{
	Vect2D<unsigned> dimIn = mtxIn.GetDim();

	mtxOp.zero.Gen(mtxOut);
	for(unsigned y=0; y<dimIn.m_y; y++)
	{
		for(unsigned x=0; x<dimIn.m_x; x++)
		{
			if(mtxIn.CellVal(x, y) > 0)
			{
				for(unsigned yy=0; yy<len; yy++)
				{
					int yLoc = (int)y + yy - len/2;
					if(yLoc<0 || yLoc>=(int)dimIn.m_y)
					{
						continue;
					}
					else {}

					for(unsigned xx=0; xx<len; xx++)
					{
						int xLoc = (int)x + xx - len/2;
						if(xLoc<0 || xLoc>=(int)dimIn.m_x)
						{
							continue;
						}
						else {}

						mtxOut.CellRef(xLoc, yLoc) = 1.F;
					}
				}
			}
			else {}
		}
	}
}

Erose::Erose()
{}

Erose::~Erose()
{}

void Erose::Gen(Mtx &mtxOut, Mtx &mtxIn, unsigned len)
{
	Vect2D<unsigned> dimIn = mtxIn.GetDim();

	mtxOp.zero.Gen(mtxOut);
	for(unsigned y=0; y<dimIn.m_y; y++)
	{
		for(unsigned x=0; x<dimIn.m_x; x++)
		{
			bool bAll = true;
			for(unsigned yy=0; yy<len; yy++)
			{
				int yLoc = (int)y + yy - len/2;
				if(yLoc<0 || yLoc>=(int)dimIn.m_y)
				{
					continue;
				}
				else {}

				for(unsigned xx=0; xx<len; xx++)
				{
					int xLoc = (int)x + xx - len/2;
					if(xLoc<0 || xLoc>=(int)dimIn.m_x)
					{
						continue;
					}
					else {}

					if(mtxIn.CellVal(xLoc, yLoc) <= 0)
					{
						bAll = false;
						break;
					} 
					else {}
				} // xx
				if(!bAll)
				{
					break;
				}
				else {}
			} // yy

			if(bAll)
			{
				mtxOut.CellRef(x, y) = 1.F;
				//cout << x << "," << y << "  ";
			}
			else {}
		} // x
	} // y
}

MorphGray::MorphGray()
{}

MorphGray::~MorphGray()
{}

void MorphGray::Gen(Mtx &mtxOut, Mtx &mtxIn, Mtx mtxKerl, bool bDilate) 
{
	Vect2D<unsigned> dimOut  = mtxOut.GetDim();
	Vect2D<unsigned> dimIn   = mtxIn.GetDim();
	Vect2D<unsigned> dimKerl = mtxKerl.GetDim();
	MyAssert(dimOut.m_x == dimIn.m_x &&
			 dimOut.m_y == dimIn.m_y);

	for (unsigned y = 0; y < dimIn.m_y; y++) {
		int yB = (int)y - dimKerl.m_y / 2;
		for (unsigned x = 0; x < dimIn.m_x; x++) {
			int xL = (int)x - dimKerl.m_x / 2;

			DATA maxV = -1e10;
			DATA minV = 1e10;
			for (unsigned yy = 0; yy < dimKerl.m_y; yy++) {
				int yLoc = yB + yy;
				if (!mtxIn.IsYInside(yLoc)) {
					continue;
				} else {}

				for (unsigned xx = 0; xx < dimKerl.m_x; xx++) {
					int xLoc = xL + xx;
					if (!mtxIn.IsXInside(xLoc)) {
						continue;
					} else {}

					if (mtxKerl.CellVal(xx, yy) != 0) {
						DATA v = mtxIn.CellVal((unsigned)xLoc, (unsigned)yLoc);
						if (v > maxV) {
							maxV = v;
						} else if (v < minV) {
							minV = v;
						} else {}
					} else {}
				} // xx
			} // yy

			if (maxV < 0) {
				maxV = 0;
			} else {}
			if (minV < 0) {
				minV = 0;
			} else {}
			mtxOut.CellRef(x, y) = (bDilate) ? maxV : minV;
		} // x
	} // y
}

//*************************************************************************************************
//
//*************************************************************************************************

DoG::DoG()
{}

DoG::~DoG()
{}

void DoG::Gen(Mtx &mtxIn, Mtx* apMtxG[], DATA s, DATA wSub, DATA sScl)
{
	mtxOp.Gauss.Gen(*apMtxG[1], s*sScl, s*sScl, true);	
	mtxOp.mul.Gen(*apMtxG[1], wSub);

	mtxOp.Gauss.Gen(*apMtxG[0], s, s, true);
	mtxOp.sub.Gen(*apMtxG[0], *apMtxG[1]);

	mtxOp.conv.Gen(mtxIn, *apMtxG[0]);
}

void DoG::GenNPR(Mtx &mtxIn, Mtx* apMtxG[], DATA s, DATA coeSharp)
{
	Vect2D<unsigned> dim = mtxIn.GetDim();

	mtxOp.mul.Gen(mtxIn, 255.F);	
	mtxOp.DoG.Gen(mtxIn, apMtxG, s);
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			if(mtxIn.CellVal(x, y) > 0)
			{
				mtxIn.CellRef(x, y) = 1.F;
			}
			else
			{
				DATA ex = exp(coeSharp * 2.F * mtxIn.CellVal(x, y));
				DATA rel = 1.F + (ex-1.F) / (ex+1.F);
				mtxIn.CellRef(x, y) = (rel>=0)? rel: 0;
			}
		}
	}	
}

//*************************************************************************************************
//
//*************************************************************************************************

CrossCorr::CrossCorr()
{}

CrossCorr::~CrossCorr()
{}

DATA CrossCorr::Gen(const Mtx &mtxIn, Mtx &mtxMsk, unsigned xOff, unsigned yOff)
{
	Vect2D<unsigned> dimMsk = mtxMsk.GetDim();
	
	Mtx mtxInTmp(mtxIn, xOff, yOff, dimMsk.m_x, dimMsk.m_y);
	DATA avgIn  = mtxOp.avg.Gen(mtxInTmp);
	DATA avgMsk = mtxOp.avg.Gen(mtxMsk);

	DATA devIn  = mtxOp.dev.Gen(mtxInTmp, false);
	DATA devMsk = mtxOp.dev.Gen(mtxMsk, false);

	DATA val = 0;
	for(unsigned y=0; y<dimMsk.m_y; y++)
	{
		for(unsigned x=0; x<dimMsk.m_x; x++)
		{
			val += (mtxInTmp.CellVal(x, y)-avgIn) * (mtxMsk.CellVal(x, y)-avgMsk);
		}
	}
	val /= devIn*devMsk;
	return val;
}

//*************************************************************************************************
//
//*************************************************************************************************

Inverse_Newton::Inverse_Newton()
{}

Inverse_Newton::~Inverse_Newton()
{}

void Inverse_Newton::Gen(Mtx &mtxOut, const Mtx &mtxIn, 
	Mtx &mtxV, Mtx &mtxTmp, unsigned itNum)
{
	Vect2D<unsigned> dimOut = mtxOut.GetDim();
	Vect2D<unsigned> dimIn  = mtxIn.GetDim();
	Vect2D<unsigned> dimTmp = mtxTmp.GetDim();
	MyAssert(dimOut.m_x == dimOut.m_y);
	MyAssert(dimOut.m_x == dimIn.m_x && 
			 dimOut.m_y == dimIn.m_y);
	MyAssert(dimTmp.m_x == dimIn.m_x &&
			 dimTmp.m_y == dimIn.m_y);

	//*******************************************
	// initial matrix
	//*******************************************
//cout << "******************************" << endl;
//cout << mtxIn.CellVal(50, 0) << " " << mtxIn.CellVal(0, 50) << endl;

	DATA maxX = -1e10;
	for (unsigned x = 0; x < dimIn.m_x; x++) {
		DATA sumX = 0;
		for (unsigned y = 0; y < dimIn.m_y; y++) {
			sumX += fabs(mtxIn.CellVal(x, y));
		}
		if (sumX > maxX) {
			maxX = sumX;
		} else {}
	}
	MyAssert(maxX > 0);

	DATA maxY = -1e10;
	for (unsigned y = 0; y < dimIn.m_y; y++) {
		DATA sumY = 0;
		for (unsigned x = 0; x < dimIn.m_x; x++) {
			sumY += fabs(mtxIn.CellVal(x, y));
		}
		if (sumY > maxY) {
			maxY = sumY;
		} else {}
	}
	MyAssert(maxY > 0);

	for (unsigned y = 0; y < dimOut.m_y; y++) {
		for (unsigned x = 0; x < dimOut.m_x; x++) {
			mtxOut.CellRef(x, y) = mtxIn.CellVal(y, x);
		}
	}
	mtxOp.mul.Gen(mtxOut, 1.F / (maxX * maxY));
	//cout << maxX << " " << maxY << endl;
	//getchar();

	//*******************************************
	// interation
	//*******************************************
	cout << "iteration" << endl;
	for (unsigned it = 0; it < itNum; it++) {
		cout << it << " ";
		if (it == itNum - 1) {
			cout << endl;
		}
		
		// V
		mtxV.CopyFrom(mtxOut);

		mtxOp.mul.Gen(mtxOut, mtxIn, mtxV);
		for (unsigned y = 0; y < dimOut.m_y; y++) {
			for (unsigned x = 0; x < dimOut.m_x; x++) {
				mtxOut.CellRef(x, y) = (x == y) ?
					2.F - mtxOut.CellVal(x, y) :
					-mtxOut.CellVal(x, y);
			}
		}

		mtxOp.mul.Gen(mtxTmp, mtxV, mtxOut);
		mtxOut.CopyFrom(mtxTmp);

		DATA vMax, vMin;
		mtxOp.rng.Gen(vMin, vMax, mtxOut);
		cout << vMin << " " << vMax << endl;

		// test if I
		mtxOp.mul.Gen(mtxTmp, mtxIn, mtxOut);
		if (it == itNum - 1) {
			//mtxOp.out << mtxTmp;
		} else {}
		bool bI = true;
		for (unsigned y = 0; y < dimTmp.m_y; y++) {
			for (unsigned x = 0; x < dimTmp.m_x; x++) {
				bI = (x == y) ?
					myMath.IsEqual(mtxTmp.CellVal(x, y), 1.F, 1e-3) :
					myMath.IsEqual(mtxTmp.CellVal(x, y), 0,   1e-3);
				if (!bI) {
					break;	
				} else {}
			} // x
			if (!bI) {
				break;
			} else {}
		} // y
		if (bI) {
			break;
		} else {}
	} // it
	cout << "iteration ok" << endl;
}

Solv_GElim::Solv_GElim()
{}

Solv_GElim::~Solv_GElim()
{}

// AX = B
int Solv_GElim::Gen(Mtx &mtxX, Mtx &mtxA, Mtx &mtxB, unsigned aIdxR[], bool bDebug)
{
	Vect2D<unsigned> dimA = mtxA.GetDim();
	Vect2D<unsigned> dimB = mtxB.GetDim();
	Vect2D<unsigned> dimX = mtxX.GetDim();
	MyAssert(dimA.m_x == dimA.m_y);
	MyAssert(dimB.m_x == 1 && dimX.m_x == 1);
	MyAssert(dimB.m_y == dimX.m_y);
	MyAssert(dimA.m_x == dimX.m_y);

	//unsigned lenIdx = sizeof(aIdxR) / sizeof(unsigned);
	//cout << lenIdx << endl;
	//MyAssert(dimX.m_y <= lenIdx);

	if (bDebug) {
		cout << "debug of Solv_GLime" << endl;
	} else {}

	for (unsigned i = 0; i < dimA.m_y; i++) {
		aIdxR[i] = i;
	}

	for (unsigned r = 0; r <= dimA.m_y - 2; r++) {
		//***************************************
		// exchange rows
		//***************************************
		DATA maxV = mtxA.CellVal(r, aIdxR[r]);
		unsigned maxIdxR = r;
		for (unsigned rr = r + 1; rr <= dimA.m_y - 1; rr++) {
			DATA v = mtxA.CellVal(r, aIdxR[rr]);
			if (v > maxV) {
				maxV = v;
				maxIdxR = rr;
			} else {}
		} // rr

		if (maxIdxR != r) {
			unsigned rTmp = aIdxR[maxIdxR];
			aIdxR[maxIdxR] = aIdxR[r];
			aIdxR[r] = rTmp;
		} else {}
		if (bDebug) {
			cout << "aIdxR: " << endl;
			for (unsigned rr = 0; rr < dimA.m_y; rr++ ) {
				cout << aIdxR[rr] << " ";
			}	
			cout << endl;
		} else {}

		bool bZPiv = false;
		DATA pivVal = mtxA.CellVal(r, aIdxR[r]);
		if (myMath.IsEqual(pivVal, 0)) {
			bZPiv = true;
		} else {}

		//***************************************
		// elimination
		//***************************************
		for (unsigned rElim = r + 1; rElim <= dimA.m_y - 1; rElim++) {
			DATA ratio = (bZPiv == false) ?
				mtxA.CellVal(r, aIdxR[rElim]) / pivVal :
				0; 

			for (unsigned col = r; col <= dimA.m_x - 1; col++) {
				mtxA.CellRef(col, aIdxR[rElim]) = 
					mtxA.CellVal(col, aIdxR[rElim]) - mtxA.CellVal(col, aIdxR[r]) * ratio;
			} // c

			mtxB.CellRef(0, aIdxR[rElim]) = 
				mtxB.CellVal(0, aIdxR[rElim]) - mtxB.CellVal(0, aIdxR[r]) * ratio;
		} // rr	
	} // r
	if (bDebug) {
		cout << "mtxA: " << endl;
		mtxOp.out << mtxA;

		cout << "mtxB: " << endl;
		mtxOp.out << mtxB;
	} else {}

	//*******************************************
	// reverse
	//*******************************************
	for (int r = dimA.m_y - 1; r >= 0; r--) {
		/*
		if (myMath.IsEqual(mtxA.CellVal(r, aIdxR[r]), 0)) {
			return -1;
		} else {}

		mtxX.CellRef(0, aIdxR[r]) = mtxB.CellVal(0, aIdxR[r]) / mtxA.CellVal(r, aIdxR[r]);
		*/
		bool bZPiv = false;
		if (!myMath.IsEqual(mtxA.CellVal(r, aIdxR[r]), 0, 1e-6)) {
			mtxX.CellRef(0, aIdxR[r]) = mtxB.CellVal(0, aIdxR[r]) / mtxA.CellVal(r, aIdxR[r]);
		} else {
			mtxX.CellRef(0, aIdxR[r]) = 0; //1.F;
			bZPiv = true;
		}

		if (!bZPiv) {
			for (int rr = 0; rr <= r - 1; rr++) {
				mtxB.CellRef(0, aIdxR[rr]) -= mtxX.CellVal(0, aIdxR[r]) * mtxA.CellVal(r, aIdxR[rr]);
			} // rr
		} else {
			for (int rr = 0; rr <= r - 1; rr++) {
				mtxB.CellRef(0, aIdxR[rr]) -= mtxB.CellVal(0, aIdxR[r]);
			} // rr
		}
	} // r 
	if (bDebug) {
		cout << "result: " << endl;
		mtxOp.out << mtxX;
	} else {}
	return 0;
}

QR_symmetric::QR_symmetric()
{}

QR_symmetric::~QR_symmetric()
{}

DATA RunDot(Mtx const &mA, Mtx const &mB)
{
	Vect2D<unsigned> dimA = mA.GetDim();
	Vect2D<unsigned> dimB = mB.GetDim();
	MyAssert(dimA.m_x == 1 && dimB.m_x == 1);
	MyAssert(dimA.m_y == dimB.m_y);

	DATA r = 0;
	for (unsigned y = 0; y < dimA.m_y; y++) {
		r += mA.CellVal(0, y) * mB.CellVal(0, y);
	}
	return r;
}

DATA RunNorm(Mtx &mA)
{
	Vect2D<unsigned> dimA = mA.GetDim();
	MyAssert(dimA.m_x == 1);

	DATA r = 0;
	for (unsigned y = 0; y < dimA.m_y; y++) {
		r += mA.CellVal(0, y) * mA.CellVal(0, y);
	}
	if (myMath.IsEqual(r, 0)) {
		return 0;
	} else {
		r = sqrt(r);
	}

	for (unsigned y = 0; y < dimA.m_y; y++) {
		mA.CellRef(0, y) /= r;
	}
	return r;
}

// Gram-Schmidt algorithm
int QR_symmetric::Gen(Mtx &mtxQ, Mtx &mtxR, Mtx const &mtxIn, bool bDebug)
{
	Vect2D<unsigned> dimIn = mtxIn.GetDim();
	Vect2D<unsigned> dimQ = mtxQ.GetDim();
	Vect2D<unsigned> dimR = mtxR.GetDim();
	MyAssert(dimQ.m_x == dimIn.m_x && dimQ.m_y == dimIn.m_y);
	MyAssert(dimR.m_x == dimIn.m_x && dimR.m_y == dimIn.m_y);

	DATA err = 1e-4;
	MyAssert(dimIn.m_x == dimIn.m_y);
	for (unsigned y = 0; y < dimIn.m_y; y++) {
		for (unsigned x = 0; x < dimIn.m_x; x++) {
			DATA fabIn = fabs(mtxIn.CellVal(x, y));
			if (fabIn > 1.F) {
				err *= fabIn;
			} else {}

			bool bZ = myMath.IsEqual(mtxIn.CellVal(x, y), mtxIn.CellVal(y, x),err);
			if (!bZ) {
				cout << mtxIn.CellVal(x, y) << " " << mtxIn.CellVal(y, x) << endl;
				MyAssert(0);
			} else {}
		} // x
	} // y

	if (bDebug) {
		cout << "qr" << endl;
	} else {}

	// compute Q
	for (unsigned x = 0; x < dimIn.m_x; x++) {
		// subtract mapping
		for (unsigned y = 0; y < dimIn.m_y; y++) {
			mtxQ.CellRef(x, y) = mtxIn.CellVal(x, y);
		}
		for (int xx = (int)x - 1; xx >= 0; xx--) {
			DATA dot = RunDot(Mtx(mtxIn, x,  0, 1, dimIn.m_y), 
							  Mtx(mtxQ,  xx, 0, 1, dimIn.m_y));

			for (unsigned y = 0; y < dimIn.m_y; y++) {
				mtxQ.CellRef(x, y) -= dot * mtxQ.CellVal(xx, y);
			}
		} // xx

		// normalize
		DATA len = RunNorm(Mtx(mtxQ, x, 0, 1, dimIn.m_y));
		if (myMath.IsEqual(len, 0)) {
			return -1;
		} else {}
	} // x

	// compute R
	mtxOp.zero.Gen(mtxR);
	for (unsigned x = 0; x < dimIn.m_x; x++) {
		for (unsigned y = 0; y <=x; y++) {
				mtxR.CellRef(x, y) = RunDot(Mtx(mtxIn, x, 0, 1, dimIn.m_y),
											Mtx(mtxQ,  y, 0, 1, dimIn.m_y));
		} // y
	} // x

	if (bDebug) {
		cout << "qr end" << endl;
	} else {}
	return 0;
}

Eigen_symmetric::Eigen_symmetric()
{}

Eigen_symmetric::~Eigen_symmetric()
{}

int Eigen_symmetric::Gen(Mtx &mtxEVal, Mtx &mtxQ, Mtx &mtxR, Mtx &mtxTmp, unsigned maxLoop, bool bDebug) 
{
	Vect2D<unsigned> dimEVal = mtxEVal.GetDim();

	//*******************************************
	// eigenvalues
	//*******************************************
	DATA zErr = 1e-4;
	mtxTmp.CopyFrom(mtxEVal);
	for (unsigned i = 0; i < maxLoop; i++) {
		int bErr = mtxOp.qr_sym.Gen(mtxQ, mtxR, mtxEVal, false);
		if (bErr != 0) {
			return bErr;
		} else {}

		mtxOp.mul.GenTA(mtxEVal, mtxQ, mtxTmp);
		mtxOp.mul.Gen(mtxTmp, mtxEVal, mtxQ);
		mtxEVal.CopyFrom(mtxTmp);
		if (bDebug) {
			cout << "mtxQ:" << endl;
			mtxOp.out << mtxQ;
			cout << "mtxR:" << endl;
			mtxOp.out << mtxR;
			cout << "mtxTmp:" << endl;
			mtxOp.out << mtxTmp;
		} else {}

		bool bUp = true;
		for (unsigned x = 0; x < dimEVal.m_x; x++) {
			for (unsigned y = x + 1; y < dimEVal.m_y ; y++) {
				if (!myMath.IsEqual(mtxEVal.CellVal(x, y), 0, zErr)) {
					bUp = false;
					//cout << bUp << " ";
					break;
				} else {}
			} // y
			if (!bUp) {
				break;
			} else {}
		} // x

		if (bUp) {
			break;
		} else {
			//mtxTmp.CopyFrom(mtxEVal);
		}
	} // i
	return 0;
}

LeastSquare::LeastSquare()
{}

LeastSquare::~LeastSquare()
{}

/*
void LeastSquare::Gen(Mtx &mtxX, Mtx &mtxA, Mtx &mtxB)
{
	Vect2D<unsigned> dimA = mtxA.GetDim();
	Vect2D<unsigned> dimB = mtxB.GetDim();
	Vect2D<unsigned> dimX = mtxX.GetDim();

	static Eigen::MatrixXf A(dimA.m_y, dimA.m_x);
	for (unsigned y = 0; y < dimA.m_y;  y++) {
		for (unsigned x = 0; x < dimA.m_x; x++) {
			A(y, x) = (int)mtxA.CellVal(x, y);
		}
	}
	//cout << A << endl;
	//getchar();

	static Eigen::VectorXf B(dimB.m_y);
	for (unsigned y = 0; y < dimB.m_y;  y++) {
		B(y) = (float)mtxB.CellVal(0, y);
	}
	//cout << B << endl;
	//getchar();

	static Eigen::VectorXf X; 
	//X = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B);
	X = (A.transpose() * A).ldlt().solve(A.transpose() * B);
	for (unsigned y = 0; y < dimX.m_y;  y++) {
		mtxX.CellRef(0, y) = X(y);
	}
}
*/

int LeastSquare::Gen(Mtx &mtxX, Mtx &mtxA, Mtx &mtxB, Mtx &mtxAA, Mtx &mtxAB, unsigned aIdx[], bool bDebug)
{
	Vect2D<unsigned> dimA = mtxA.GetDim();
	Vect2D<unsigned> dimB = mtxB.GetDim();
	Vect2D<unsigned> dimX = mtxX.GetDim();
	//MyAssert(dimA.m_x == dimA.m_y);
	MyAssert(dimB.m_x == 1 && dimX.m_x == 1);
	MyAssert(dimB.m_y == dimA.m_y);
	MyAssert(dimX.m_y == dimA.m_x);

	Vect2D<unsigned> dimAA = mtxAA.GetDim();
	Vect2D<unsigned> dimAB = mtxAB.GetDim();
	MyAssert(dimAA.m_x == dimAA.m_y);
	MyAssert(dimAA.m_y == dimA.m_x);
	MyAssert(dimAB.m_x == 1 && dimAB.m_y == dimA.m_x);

	if (bDebug) {
		cout << "debug of LeastSquare" << endl;
		cout << "mtxA:" << endl;
		mtxOp.out << mtxA;
		cout << "mtxB:" << endl;
		mtxOp.out << mtxB;
	} else {}

	mtxOp.mul.GenTA(mtxAA, mtxA);
	if (bDebug) {
		cout << "mtxAA:" << endl;
		mtxOp.out << mtxAA;
	} else {}

	mtxOp.mul.GenTA(mtxAB, mtxA, mtxB);
	if (bDebug) {
		cout << "mtxAB:" << endl;
		mtxOp.out << mtxAB;
	} else {}

	// replace the inversion
	int errNo = mtxOp.solv_GElim.Gen(mtxX, mtxAA, mtxAB, aIdx);
	for (unsigned i = 0; i < dimX.m_y; i++) {
		mtxAA.CellRef(i, 0) = mtxX.CellVal(0, aIdx[i]);
	}
	for (unsigned i = 0; i < dimX.m_y; i++) {
		mtxX.CellRef(0, i) = mtxAA.CellVal(i, 0);
	}

	/*
	Mtx mtxAAt(dimAA.m_x, dimAA.m_y);
	Mtx mtxT0(dimAA.m_x, dimAA.m_y);
	Mtx mtxT1(dimAA.m_x, dimAA.m_y);
	Mtx mtxT2(dimAA.m_x, dimAA.m_y);
	mtxOp.inverse_Newton.Gen(mtxAAt, mtxAA, mtxT0, mtxT1);
	cout << "inverse ok" << endl;
	mtxOp.mul.Gen(mtxX, mtxAAt, mtxAB);
	int errNo = 0;
	*/

	if (bDebug) {
		if (errNo != 0) {
			cout << "result:" << endl;
			mtxOp.out << mtxX;
		} else {
			cout << "no result" << endl;
		}
	} else {}
	return errNo;
}

RegionLabel::RegionLabel()
{}

RegionLabel::~RegionLabel()
{}

DATA VCell(Mtx &mtx, int x, int y) {
	if (!mtx.IsInside(x, y)) {
		return 0;
	} else {}

	return mtx.CellVal(x, y);
}
unsigned FindRoot(vector<unsigned> arr, unsigned i)
{
	unsigned idx = i;
	unsigned v = arr[idx];
	while (v != idx) {
		idx = v;
		v = arr[idx];
	}
	return idx;
}
void RegionLabel::Gen(Mtx &mtxLab, Mtx &mtxIn, vector<unsigned> &idx)
{
	Vect2D<unsigned> dimIn  = mtxIn.GetDim();
	Vect2D<unsigned> dimLab = mtxLab.GetDim();
	MyAssert(dimIn.m_x == dimLab.m_x &&
		  	 dimIn.m_y == dimLab.m_y);

	//*******************************************
	// first pass
	//*******************************************
	idx.clear();
	idx.push_back(0);
	mtxOp.zero.Gen(mtxLab);
	unsigned maxNo = 0;
	for (unsigned y = 0; y < dimIn.m_y; y++) {
		for (unsigned x = 0; x < dimIn.m_x; x++) {
			if (mtxIn.CellVal(x, y) == 0) {
				continue;
			} else {}

			int aL[] = {(int)x - 1, y};
			int aT[] = {x, (int)y - 1};
			DATA vL = VCell(mtxLab, aL[0], aL[1]);
			DATA vT = VCell(mtxLab, aT[0], aT[1]);
			if (vL == 0 && vT == 0) {
				maxNo++;
				mtxLab.CellRef(x, y) = maxNo;
				idx.push_back(maxNo);
			} else if (vL == 0) {
				mtxLab.CellRef(x, y) = (unsigned)vT;
			} else if (vT == 0) {
				mtxLab.CellRef(x, y) = (unsigned)vL;
			} else {				
				unsigned vMin = (vT < vL) ? (unsigned)vT : (unsigned)vL;
				unsigned vMax = (vT > vL) ? (unsigned)vT : (unsigned)vL;

				mtxLab.CellRef(x, y) = vMin;
				if (vMin != vMax) {
					idx[vMax] = vMin;
				} else {}
			}
			MyAssert(idx.size() == maxNo + 1);
		} // x
	} // y

	//*******************************************
	// second pass
	//*******************************************
	for (unsigned y = 0; y < dimIn.m_y; y++) {
		for (unsigned x = 0; x < dimIn.m_x; x++) {
			unsigned i = (unsigned)mtxLab.CellVal(x, y);
			i = FindRoot(idx, i);
			mtxLab.CellRef(x, y) = i;
		}
	}

	//for (unsigned i = 0; i < idx.size(); i++) {
	//	cout << i << ":" << idx[i] << " ";
	//}
	//cout << endl;
}

//*************************************************************************************************
//
//*************************************************************************************************

void Rotate2D(DATA aOut[], const DATA aIn[], DATA ang, const DATA aCnt[])
{
	DATA aTmp[2];
	aTmp[0] = aIn[0] - aCnt[0];
	aTmp[1] = aIn[1] - aCnt[1];

	aOut[0] =  cos(ang)*aTmp[0] + sin(ang)*aTmp[1];
	aOut[1] = -sin(ang)*aTmp[0] + cos(ang)*aTmp[1];

	aOut[0] += aCnt[0];
	aOut[1] += aCnt[1];
}

DATA Length2D(DATA x, DATA y)
{
	return sqrt( (DATA)(x*x + y*y) );
}

//*************************************************************************************************