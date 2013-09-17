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
	Vect2D<unsigned> dim = mtxBase.GetDim();

	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
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
	assert(dimR.m_x==dimB.m_x && dimR.m_y==dimA.m_y && dimA.m_x==dimB.m_y);

	for(unsigned y=0; y<dimR.m_y; y++)
	{
		for(unsigned x=0; x<dimR.m_x; x++)
		{
			DATA r = 0;
			for(unsigned c=0; c<dimA.m_x; c++)
			{
				r += mtxA.CellVal(c, y) * mtxB.CellVal(x, c);
			}
			mtxR.CellRef(x, y) = r;
		}
	}
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
	for(unsigned i=0; i<m_levelNum; i++)
	{
		if(m_pCount[i] > 0)
		{
			cMin = m_pCount[i];
			break;
		}
		else {}
	}

	unsigned mtxSize = sNum; 
	DATA invDiv = 1.F / (mtxSize - cMin);
	Vect2D<unsigned> dim = mtx.GetDim();
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			if(pMSkip && pMSkip->CellVal(x, y)==0)
			{
				continue;
			}
			else {}

			DATA val = mtx.CellVal(x, y);
			unsigned level = (unsigned)(val * (m_levelNum-1));

			if(m_pCount[level] == 0)
			{
				mtx.CellRef(x, y) = 0;
			}
			else 
			{
				mtx.CellRef(x, y) = (m_pCount[level] - cMin) * invDiv;
			}
		}
	}
}

void HisEqual::Gen(Mtx &mtx, unsigned levelNum, bool bReval, Mtx *pMSkip)
{
	if(levelNum != m_levelNum)
	{
		if(m_memSize < levelNum)
		{
			delete []m_pCount;
			m_memSize = levelNum;
			m_pCount = new unsigned[m_memSize];
		}
		else {}
		m_levelNum = levelNum;
	}
	else {}

	for(unsigned i=0; i<m_levelNum; i++)
	{
		m_pCount[i] = 0;
	}
	Vect2D<unsigned> dim = mtx.GetDim();
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			if(pMSkip && pMSkip->CellVal(x, y)==0)
			{
				continue;
			}
			else {}

			unsigned level = (unsigned)(mtx.CellVal(x, y) * (m_levelNum-1));
			m_pCount[level]++; 
		}
	}

	unsigned sNum = m_pCount[0];
	for(unsigned i=1; i<m_levelNum; i++)
	{
		sNum += m_pCount[i];
		m_pCount[i] += m_pCount[i-1];
	}

	if(bReval)
	{
		Reval(mtx , pMSkip);
	}
	else {}
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

Thrd::Thrd()
{}

Thrd::~Thrd()
{}

void Thrd::Gen(Mtx &mtx, DATA thrd)
{
	Vect2D<unsigned> dim = mtx.GetDim();
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			if(mtx.CellVal(x, y) >= thrd)
			{
				mtx.CellRef(x, y) = 1.F;
			}
			else
			{
				mtx.CellRef(x, y) = 0;
			}
		}
	}
}

//*************************************************************************************************

Otsu::Otsu()
{}

Otsu::~Otsu()
{}

DATA Otsu::Gen(Mtx &mtx, DATA &thrd, unsigned lev, DATA valMax, DATA valMin)
{
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

	mtxOp.thrd.Gen(mtx, thrd);

	DATA conf = mtxDelta.CellVal(0, iLoc) / dg;
	return conf;
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
			}
			else {}
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