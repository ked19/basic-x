#include "sift.h"

//*************************************************************************************************

Octave::Octave(Mtx &mtxInit, unsigned sclPerOct, unsigned octNum, DATA dev)
	:m_sclPerOct(sclPerOct), m_octNum(octNum)
	,m_GaussNum(m_sclPerOct+3), m_GaussLen(41)
{
	cout << "initialize octave.." << endl;
	cout << "number of oct: " << m_octNum << "\n"
		 << "scale per oct: " << m_sclPerOct << "\n"
		 << "dev: " << dev << endl;
	
	//*******************************************

	m_pDev = new DATA[m_GaussNum];
	DATA k = pow(2.F, 1.F/sclPerOct); 
	
	m_pDev[0] = dev;
	for(unsigned i=1; i<m_GaussNum; i++)
	{
		m_pDev[i] = m_pDev[i-1] * k;
	}
	
	/*
	DATA devPre = dev / k;
	m_pDev[0] = devPre;
	DATA devG = dev;
	for(unsigned i=1; i<m_GaussNum; i++)
	{
		m_pDev[i] = devG;
		devG *= k;
	}
	*/
	//*******************************************

	m_ppGauss = new Mtx*[m_GaussNum];	
	for(unsigned i=0; i<m_GaussNum; i++)
	{
		DATA d = m_pDev[i];
		Mtx *pG = new Mtx(m_GaussLen, m_GaussLen);
		mtxOp.Gauss.Gen(*pG, d, d, true);
		m_ppGauss[i] = pG;
	}

	//*******************************************

	Vect2D<unsigned> dimInit = mtxInit.GetDim();
	Vect2D<unsigned> dimNow = dimInit;
	m_pppBlur = new Mtx**[m_octNum];
	for(unsigned i=0; i<m_octNum; i++)
	{
		m_pppBlur[i] = new Mtx*[m_GaussNum];
		for(unsigned j=0; j<m_GaussNum; j++)
		{
			m_pppBlur[i][j] = new Mtx(dimNow);
		}
		dimNow.m_x /= 2;		
		dimNow.m_y /= 2;
	}

	for(unsigned i=0; i<m_octNum; i++)
	{
		mtxOp.scaleDim.Gen(*m_pppBlur[i][0], mtxInit);
		for(unsigned j=1; j<m_GaussNum; j++)
		{
			m_pppBlur[i][j]->CopyFrom(*m_pppBlur[i][0]);
		}

		for(unsigned j=0; j<m_GaussNum; j++)
		{
			mtxOp.conv.Gen(*m_pppBlur[i][j], *m_ppGauss[j]);
		}
	}

	//*******************************************

	m_pppOrnt = new Layer**[m_octNum];
	for(unsigned i=0; i<m_octNum; i++)
	{
		m_pppOrnt[i] = new Layer*[m_GaussNum];
		for(unsigned j=0; j<m_GaussNum; j++)
		{
			Vect2D<unsigned> dimBlur = m_pppBlur[i][j]->GetDim();
			m_pppOrnt[i][j] = new Layer(dimBlur.m_x, dimBlur.m_y, 2);
		}
	}

	for(unsigned i=0; i<m_octNum; i++)
	{
		for(unsigned j=0; j<m_GaussNum; j++)
		{
			lyrOp.ornt.Gen(*m_pppOrnt[i][j], *m_pppBlur[i][j]);
		}
	}

	//*******************************************

	unsigned DoGNum = m_GaussNum - 1;
	m_pppDoG = new Mtx**[m_octNum];
	for(unsigned i=0; i<m_octNum; i++)
	{
		m_pppDoG[i] = new Mtx*[DoGNum];
		for(unsigned j=0; j<DoGNum; j++)
		{
			Mtx *pBlurNow = m_pppBlur[i][j+1];
			Mtx *pBlurPre = m_pppBlur[i][j];
			
			Mtx *pDoG = new Mtx(pBlurNow->GetDim());
			pDoG->CopyFrom(*pBlurNow);
			mtxOp.sub.Gen(*pDoG, *pBlurPre);

			m_pppDoG[i][j] = pDoG;
		}
	}

	cout << "octave ok" << endl << endl;
}

Octave::~Octave()
{
	for(unsigned i=0; i<m_GaussNum; i++)
	{
		delete m_ppGauss[i];
	}
	delete []m_ppGauss;

	//*******************************************

	delete []m_pDev;

	//*******************************************

	for(unsigned i=0; i<m_octNum; i++)
	{
		for(unsigned j=0; j<m_GaussNum; j++)
		{
			delete m_pppBlur[i][j];
		}
		delete []m_pppBlur[i];
	}
	delete []m_pppBlur;

	//*******************************************

	for(unsigned i=0; i<m_octNum; i++)
	{
		for(unsigned j=0; j<m_GaussNum; j++)
		{
			delete m_pppOrnt[i][j];
		}
		delete []m_pppOrnt[i];
	}
	delete []m_pppOrnt;

	//*******************************************

	unsigned DoGNum = m_GaussNum - 1;
	for(unsigned i=0; i<m_octNum; i++)
	{
		for(unsigned j=0; j<DoGNum; j++)
		{
			delete m_pppDoG[i][j];
		}
		delete []m_pppDoG[i];
	}
	delete []m_pppDoG;
}

Mtx* Octave::GetBlur(unsigned octNo, unsigned GaussNo)
{
	return m_pppBlur[octNo][GaussNo];
}

Mtx* Octave::GetDoG(unsigned octNo, unsigned DoGNo)
{
	return m_pppDoG[octNo][DoGNo];
}

Layer* Octave::GetOrnt(unsigned octNo, unsigned GaussNo)
{
	return m_pppOrnt[octNo][GaussNo];
}

unsigned Octave::GetGaussNum()
{
	return m_GaussNum;
}

//*************************************************************************************************
//
//*************************************************************************************************

SmpWin::SmpWin(Octave &oct)
	:m_GaussNum(oct.m_GaussNum), m_GaussLen(45)
{
	m_pDev = new DATA[m_GaussNum];
	for(unsigned i=0; i<m_GaussNum; i++)
	{
		m_pDev[i] = oct.m_pDev[i] * 1.5F;
	}

	m_ppGauss = new Mtx*[m_GaussNum];
	for(unsigned i=0; i<m_GaussNum; i++)
	{
		DATA d = m_pDev[i];
		Mtx *pG = new Mtx(m_GaussLen, m_GaussLen);
		mtxOp.Gauss.Gen(*pG, d, d, true);
		m_ppGauss[i] = pG;
	}
}

SmpWin::~SmpWin()
{
	delete []m_pDev;

	for(unsigned i=0; i<m_GaussNum; i++)
	{
		delete m_ppGauss[i];
	}
	delete []m_ppGauss;
}

Mtx* SmpWin::GetGauss(unsigned GaussNo)
{
	return m_ppGauss[GaussNo];
}

//*************************************************************************************************
//
//*************************************************************************************************

Keypoint::Keypoint()
{}

Keypoint::~Keypoint()
{}

//*************************************************************************************************
//
//*************************************************************************************************

Sift::Sift()
	:m_pMtxIn(0), m_pOct(0), m_pOrntSwin(0)
	,SCALE_PER_OCTAVE(3), DEV(1.6F), CONTRAST_THRHLD(0.03F), EDGE_RATIO(12.1F) 
	,DIR_RATIO(0.8F), ORNT_BIN(36)
	,DSPT_TILE(4), DSPT_ORNT(8)
	,m_invOrnt(3, 3)
	,MIN_X(32), MIN_Y(32)
{
	DATA aaInv[3][3] = {{0.5F, -1.F, 0.5F}, {-0.5F, 0, 0.5F}, {0, 1.F, 0}};
	for(unsigned y=0; y<3; y++)
	{
		for(unsigned x=0; x<3; x++)
		{
			m_invOrnt.CellRef(x, y) = aaInv[y][x];
		}
	}

	m_orntRange = 360.F / ORNT_BIN;
	m_pDir = new DATA[ORNT_BIN];
}

Sift::~Sift()
{
	delete m_pMtxIn;
	delete m_pOct;
	delete m_pOrntSwin;
	delete []m_pDir;
}

void Sift::NormalizeCell()
{
	Vect2D<unsigned> dim = m_pMtxIn->GetDim();

	DATA max = -1e10;
	DATA min =  1e10;
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			DATA val =  m_pMtxIn->CellVal(x, y);		
			if(max < val) max = val;		else {}
			if(min > val) min = val;		else {}
		}
	}
	if(myMath.IsEqual(max, min)) return;		else {}
	
	DATA rangeInv = 1.F / (max - min);
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			m_pMtxIn->CellRef(x, y) = (m_pMtxIn->CellVal(x, y) - min) * rangeInv;
		}
	}
	
	mtxOp.hisEqual.Gen(*m_pMtxIn, 256);
}

unsigned Sift::ComputeOctNum(const Mtx &mtx) const
{
	Vect2D<unsigned> dim = mtx.GetDim();
	
	unsigned octNum = 1;
	unsigned xD = MIN_X;
	unsigned yD = MIN_Y;
	while(1)
	{
		if(xD>dim.m_x || yD>dim.m_y) break;
		else {}

		octNum++;
		xD *= 2;
		yD *= 2;
	}
	return octNum - 1;
}

vector<Keypoint> Sift::Gen(Mtx &mtxIn)
{
	cout << "sift...." << endl;

	unsigned octNum = ComputeOctNum(mtxIn);
	cout << "octNum = " << octNum << endl;

	delete m_pMtxIn;
	m_pMtxIn = new Mtx(mtxIn);
	NormalizeCell();

	delete m_pOct;
	m_pOct = new Octave(*m_pMtxIn, SCALE_PER_OCTAVE, octNum, DEV);
	unsigned GaussNum = m_pOct->GetGaussNum();
	
	delete m_pOrntSwin;
	m_pOrntSwin = new SmpWin(*m_pOct);

	m_vKeypoint.clear();
	for(m_octNo=0; m_octNo<octNum; m_octNo++)
	{
		for(m_DoGNo=1; m_DoGNo<GaussNum-2; m_DoGNo++)
		{
			vector<Vect2D<unsigned>> vExtre;
			GetExtrema(vExtre);

			vector<Vect3D<DATA>> vKeyLoc;
			RefineKeypoint(vKeyLoc, vExtre);

			vector<Keypoint> vKp;
			AssignOrient(vKp, vKeyLoc);

			DescriptNbr(vKp);
		}
	}

	ofstream outF("keypoint.txt");
	for(unsigned i=0; i<m_vKeypoint.size(); i++)
	{
		Keypoint *pKp = &m_vKeypoint[i];
		outF << pKp->m_x << ", " << pKp->m_y << ", " << pKp->m_ornt;
		outF << "\n\t";
		for(unsigned y=0; y<DSPT_TILE; y++)
		{
			for(unsigned x=0; x<DSPT_TILE; x++)
			{
				for(unsigned d=0; d<DSPT_ORNT; d++)
				{
					outF << pKp->m_aaaDspt[y][x][d] << ", ";
				}
			}
		}
		outF << "\n" << endl;
	}
	outF.close();

	return m_vKeypoint;
}

//*************************************************************************************************

void Sift::GetExtrema(vector<Vect2D<unsigned>> &vExtre)
{
	cout << "get extrema.. ";

	vExtre.clear();

	Mtx *pDoGPre = m_pOct->GetDoG(m_octNo, m_DoGNo-1);
	Mtx *pDoGNow = m_pOct->GetDoG(m_octNo, m_DoGNo);
	Mtx *pDoGNxt = m_pOct->GetDoG(m_octNo, m_DoGNo+1);

	Vect2D<unsigned> dim = pDoGNow->GetDim();
	for(unsigned y=0; y<dim.m_y; y++)
	{
		int yB = y - 1;			if(yB<0) yB = 0;	else {}
		unsigned yT = y + 1;	if(yT>=dim.m_y) yT = dim.m_y - 1;	else{}
		for(unsigned x=0; x<dim.m_x; x++)
		{
			int xL = x - 1;			if(xL<0) xL = 0;	else {}
			unsigned xR = x + 1;	if(xR>=dim.m_x)	xR = dim.m_x - 1;	else{}

			DATA valNow = pDoGNow->CellVal(x, y);

			bool bMax = true;
			bool bMin = true;
			for(unsigned yy=(unsigned)yB; yy<=yT; yy++)
			{
				for(unsigned xx=(unsigned)xL; xx<=xR; xx++)
				{
					if(myMath.IsELess (valNow, pDoGPre->CellVal(xx, yy))) bMax = false;		else {}
					if(myMath.IsELess (valNow, pDoGNxt->CellVal(xx, yy))) bMax = false;		else {}
					if(myMath.IsEGreat(valNow, pDoGPre->CellVal(xx, yy))) bMin = false;		else {}
					if(myMath.IsEGreat(valNow, pDoGNxt->CellVal(xx, yy))) bMin = false;		else {}

					if(yy==y && xx==x) continue;	else {}
					if(myMath.IsELess (valNow, pDoGNow->CellVal(xx, yy))) bMax = false;		else {}
					if(myMath.IsEGreat(valNow, pDoGNow->CellVal(xx, yy))) bMin = false;		else {}

				}
			}

			
			if(bMax || bMin)
			{
				vExtre.push_back(Vect2D<unsigned>(x, y));
			}
			else {}
		}
	}

	cout << "ok" << endl;
}

//*************************************************************************************************

bool Sift::IsLowContrast(DATA smpVal, const Mtx &mtxDer, const Mtx &mtxOpt)
{
	DATA dVal = smpVal + (mtxDer.CellVal(0, 0) * mtxOpt.CellVal(0, 0) + 
						  mtxDer.CellVal(0, 1) * mtxOpt.CellVal(0, 1) +
						  mtxDer.CellVal(0, 2) * mtxOpt.CellVal(0, 2)) * 0.5F;


	if(myMath.IsLess(fabs(dVal), CONTRAST_THRHLD))	return true;
	else											return false;
}

bool Sift::IsEdge(Mtx &mtx, unsigned x, unsigned y)
{
	Mtx mtxH(2, 2);
	mtxOp.Hess.Gen(mtxH, mtx, x, y);

	DATA dx2  = mtxH.CellVal(0, 0);
	DATA dy2  = mtxH.CellVal(1, 1);
	DATA dxdy = mtxH.CellVal(1, 0);
	
	DATA tr = dx2 + dy2;
	DATA det = dx2*dy2 - dxdy*dxdy;
	
	if(myMath.IsEqual(det, 0))
	{
		return true;
	}
	else {}

	DATA r = tr*tr / det;
	if(myMath.IsELess(r, EDGE_RATIO))	return false;
	else								return true;
}

void Sift::RefineKeypoint(vector<Vect3D<DATA>> &vKeyLoc , vector<Vect2D<unsigned>> &vExtre)
{
	cout << "refine keypoint.. ";

	vKeyLoc.clear();
	for(unsigned i=0; i<vExtre.size(); i++)
	{
		unsigned x = vExtre[i].m_x;
		unsigned y = vExtre[i].m_y;

		Mtx *pDoGPre = m_pOct->GetDoG(m_octNo, m_DoGNo-1);
		Mtx *pDoGNow = m_pOct->GetDoG(m_octNo, m_DoGNo);
		Mtx *pDoGNxt = m_pOct->GetDoG(m_octNo, m_DoGNo+1);
		Mtx aM[] = {Mtx(*pDoGPre), Mtx(*pDoGNow), Mtx(*pDoGNxt)};
		Layer lyrDoG(aM, 3);

		Vect2D<unsigned> dim = pDoGNow->GetDim();
		if(x<=0 || x>=dim.m_x-1 ||
		   y<=0 || y>=dim.m_y-1)
		{
			continue;
		}
		else {}
		
		Mtx mtxHess(3, 3);			lyrOp.Hess.Gen(mtxHess, lyrDoG, x, y, 1);
		Mtx mtxHInv(3, 3);			gslMtxOp.inv.Gen(mtxHInv, mtxHess);
		Mtx mtxDer(1, 3);			lyrOp.der.Gen(mtxDer, lyrDoG, x, y, 1);
		Mtx mtxOffset(1, 3);		mtxOp.mul.Gen(mtxOffset, mtxHInv, mtxDer);
		mtxOp.mul.Gen(mtxOffset, -1.F);

		if(mtxOffset.CellVal(0, 0)<-1.F || mtxOffset.CellVal(0, 0)>1.F ||
		   mtxOffset.CellVal(0, 1)<-1.F || mtxOffset.CellVal(0, 1)>1.F ||
		   mtxOffset.CellVal(0, 2)<-1.F || mtxOffset.CellVal(0, 2)>1.F)
		{
			mtxOp.zero.Gen(mtxOffset);
			//continue;
		}
		else {}

		unsigned ss = 1 << m_octNo;
		unsigned ssOffset = ss - 1;
		DATA xx = x*ss + ssOffset;
		DATA yy = y*ss + ssOffset;

		if(IsLowContrast(lyrDoG.CellVal(x, y, 1), mtxDer, mtxOffset))
		{
			continue;
		}
		else {}

		DATA xRefine = x + mtxOffset.CellVal(0, 0);
		DATA yRefine = y + mtxOffset.CellVal(0, 1);
		DATA cRefine = mtxOffset.CellVal(0, 2); 

		unsigned xPix = (unsigned)(xRefine + 0.5F);
		unsigned yPix = (unsigned)(yRefine + 0.5F);
		unsigned cPix = (unsigned)(1 + cRefine + 0.5F);
		if(IsEdge(*lyrDoG.GetMtx(cPix), xPix, yPix)) 
		{
			continue;
		}
		else {}
		
		Vect3D<DATA> kLoc(xRefine, yRefine, cRefine);
		vKeyLoc.push_back(kLoc);
	}
	vExtre.clear();

	cout << "ok" << endl;
}

//*************************************************************************************************

DATA Sift::LDir(DATA *pDir, unsigned now)
{
	unsigned l = (now + 1) % ORNT_BIN;
	return pDir[l];
}
DATA Sift::RDir(DATA *pDir, unsigned now)
{
	unsigned r = (now - 1 + ORNT_BIN) % ORNT_BIN;
	return pDir[r];
}

DATA Sift::SolveOrnt(DATA lVal, DATA nVal, DATA rVal, unsigned idx)
{
	Mtx mtxVal(1, 3);
	mtxVal.CellRef(0, 0) = rVal;
	mtxVal.CellRef(0, 1) = nVal;
	mtxVal.CellRef(0, 2) = lVal;

	Mtx mtxR(1, 3);
	mtxOp.mul.Gen(mtxR, m_invOrnt, mtxVal);	
	DATA offset = myMath.IsEqual(mtxR.CellVal(0, 0), 0)?
		0: -mtxR.CellVal(0, 1) / mtxR.CellVal(0, 0) * 0.5F;

	DATA ornt = (offset + idx) * m_orntRange * D2R;
	return ornt;
}

void Sift::AssignOrient(vector<Keypoint> &vKp, vector<Vect3D<DATA>> &vKeyLoc)
{	
	cout << "assign orient.. ";

	DATA d2rScl = m_orntRange * D2R;
	DATA invRange = 1.F / m_orntRange;

	unsigned newKp = 0;
	vKp.clear();
	for(unsigned i=0; i<vKeyLoc.size(); i++)
	{
		for(unsigned j=0; j<ORNT_BIN; j++)
		{
			m_pDir[j] = 0;
		}

		DATA cOffset = vKeyLoc[i].m_z;
		unsigned GaussNo = (cOffset>=0)?	m_DoGNo+1: m_DoGNo;
		const Layer *pOrnt = m_pOct->GetOrnt(m_octNo, GaussNo);
		Mtx *pOrntSwin = m_pOrntSwin->GetGauss(GaussNo);
	
		Vect2D<unsigned> dimSmpwin = pOrntSwin->GetDim();
		int xR = (dimSmpwin.m_x-1) / 2;		int xL = -xR;
		int yT = (dimSmpwin.m_y-1) / 2;		int yB = -yT;

		unsigned xCnt = (unsigned)(vKeyLoc[i].m_x + 0.5F);		
		unsigned yCnt = (unsigned)(vKeyLoc[i].m_y + 0.5F);
		Vect3D<unsigned> dimOrnt = pOrnt->GetDim();
		for(int y=yB, ySmpLoc=0; y<=yT; y++, ySmpLoc++)
		{
			int yLoc = yCnt + y;
			if(yLoc<0 || yLoc>=(int)dimOrnt.m_y) 
			{
				continue;
			}
			else {}

			for(int x=xL, xSmpLoc=0; x<=xR; x++, xSmpLoc++)
			{
				int xLoc = xCnt + x;
				if(xLoc<0 || xLoc>=(int)dimOrnt.m_x)
				{
					continue;
				}
				else {}

				DATA dir = pOrnt->CellVal(xLoc, yLoc, 1) * R2D + 360.F;
				unsigned dirIdx = (unsigned)(dir*invRange + 0.5F) % ORNT_BIN;

				DATA mag = pOrnt->CellVal(xLoc, yLoc, 0);
				m_pDir[dirIdx] += mag * pOrntSwin->CellVal(xSmpLoc, ySmpLoc);
			}
		}

		unsigned maxIdx = 0;
		for(unsigned j=1; j<ORNT_BIN; j++)
		{
			if(m_pDir[j] > m_pDir[maxIdx])
			{
				maxIdx = j;
			}
			else{}
		}
		DATA dirThrhld = m_pDir[maxIdx] * DIR_RATIO;
		
		for(unsigned j=0; j<ORNT_BIN; j++)
		{
			DATA nVal = m_pDir[j];
			DATA lVal = LDir(m_pDir, j);
			DATA rVal = RDir(m_pDir, j);
			if(lVal>=nVal || rVal>=nVal)
			{
				continue;
			}
			else{}

			if(m_pDir[j] > dirThrhld)
			{
				DATA ornt = SolveOrnt(lVal, nVal, rVal, j);

				Keypoint kp;
				kp.m_octNo		= m_octNo;
				kp.m_GaussNo	= GaussNo;
				kp.m_x			= vKeyLoc[i].m_x;
				kp.m_y			= vKeyLoc[i].m_y;
				kp.m_ornt		= ornt;
				vKp.push_back(kp);
				++newKp;
			}
			else {}
		}
	}

	DATA addPercent = (vKeyLoc.size())?
		(DATA)newKp / vKeyLoc.size(): 1.F;
	cout << "add " << addPercent << " times ok" << endl;

	vKeyLoc.clear();
}

//*************************************************************************************************

void Sift::ZeroDspt(Keypoint &kp)
{
	for(unsigned y=0; y<DSPT_TILE; y++)
	{
		for(unsigned x=0; x<DSPT_TILE; x++)
		{
			for(unsigned d=0; d<DSPT_ORNT; d++)
			{
				kp.m_aaaDspt[y][x][d] = 0;
			}
		}
	}
}

void Sift::DescriptNbr(vector<Keypoint> &vKp)
{	
	unsigned scl = 1 << m_octNo;
	unsigned sOffset = scl - 1;

	DATA orntWidth = 2.F*PI / DSPT_ORNT;
	for(unsigned i=0; i<vKp.size(); i++)
	{
		Keypoint *pKp = &vKp[i];
		Layer *pOrnt = m_pOct->GetOrnt(pKp->m_octNo, pKp->m_GaussNo);
		Vect3D<unsigned> dim = pOrnt->GetDim();

		unsigned xCnt = (unsigned)(pKp->m_x + 0.5F);
		unsigned yCnt = (unsigned)(pKp->m_y + 0.5F);
		int xL = xCnt - 2*DSPT_TILE;	// xL = xCnt - 4 / 2 * DSPT_TILE
		int yB = yCnt - 2*DSPT_TILE;

		ZeroDspt(*pKp);
		DATA orntKp = pKp->m_ornt;
		for(unsigned yTiIdx=0, yOffset=0; yTiIdx<DSPT_TILE; yTiIdx++, yOffset+=4)
		{
			for(unsigned xTiIdx=0, xOffset=0; xTiIdx<DSPT_TILE; xTiIdx++, xOffset+=4)
			{
				for(unsigned yy=0; yy<4; yy++)
				{
					int yCount = yB + yOffset + yy;
					int yTiNAdd = (yy<2)? -1: 1;

					for(unsigned xx=0; xx<4; xx++)
					{
						int xCount = xL + xOffset + xx;
						int xTiNAdd = (xx<2)? -1: 1;

						DATA aIn[] = {xCount, yCount};
						DATA aCnt[] = {xCnt, yCnt};
						DATA aOut[2];
						Rotate2D(aOut, aIn, -orntKp, aCnt);
						int xLoc = (aOut[0]>0)? (int)(aOut[0]+0.5F): (int)(aOut[0]-0.5F);
						int yLoc = (aOut[1]>0)? (int)(aOut[1]+0.5F): (int)(aOut[1]-0.5F);
						if(xLoc<0 || xLoc>=(int)dim.m_x ||
						   yLoc<0 || yLoc>=(int)dim.m_y)
						{
							continue;
						}
						else {}

						DATA magCell = pOrnt->CellVal(xLoc, yLoc, 0);
						//if(magCell > 0.2F) magCell = 0.2F;		else {}

						DATA orntCell = pOrnt->CellVal(xLoc, yLoc, 1) - orntKp;
						while(orntCell < 0)
						{
							orntCell += 2 * PI;
						}
						orntCell /= orntWidth;
						if(orntCell == DSPT_ORNT)
						{
							orntCell = DSPT_ORNT - 1;
						}
						else {}
						assert(orntCell>=0 && orntCell<DSPT_ORNT);

						unsigned oNow = (unsigned)(orntCell + 0.5F);
						int oNbr = ((DATA)oNow > orntCell)? 
							oNow-1: oNow+1;
						DATA ratioNow = 1.F - fabs(orntCell - oNow);
						DATA ratioNbr = 1.F - fabs(orntCell - oNbr);
						oNow = oNow % DSPT_ORNT;
						oNbr = (oNbr + DSPT_ORNT) % DSPT_ORNT;

						for(int yA=0, yN=0; yA<2; yA++, yN=yTiNAdd)
						{
							int yTiNIdx = yTiIdx + yN;
							if(yTiNIdx<0 || yTiNIdx>=(int)DSPT_TILE) 
							{
								continue;
							}
							else {}
							int yDiff = yN*4 + 2 - yy;

							for(int xA=0, xN=0; xA<2; xA++, xN=xTiNAdd)
							{
								int xTiNIdx = xTiIdx + xN;
								if(xTiNIdx<0 || xTiNIdx>=(int)DSPT_TILE)
								{
									continue;
								}
								else {}

								int xDiff = xN*4 + 2 - xx;
								DATA dis = Length2D(xDiff, yDiff);
								DATA ratioDis = 1.F - dis/6.F;

								pKp->m_aaaDspt[yTiNIdx][xTiNIdx][oNow] += ratioNow * ratioDis * magCell;
								pKp->m_aaaDspt[yTiNIdx][xTiNIdx][oNbr] += ratioNbr * ratioDis * magCell;
							}
						}
					}
				}
			}
		}

		bool bConti = true;
		while(bConti)
		{
			DATA dsptSum = 0;
			for(unsigned y=0; y<DSPT_TILE; y++)
			{
				for(unsigned x=0; x<DSPT_TILE; x++)
				{
					for(unsigned d=0; d<DSPT_ORNT; d++)
					{
						dsptSum += pKp->m_aaaDspt[y][x][d] * pKp->m_aaaDspt[y][x][d];
					}
				}
			}

			bConti = false;
			if(dsptSum)
			{
				DATA invS = 1.F / sqrt(dsptSum);
				for(unsigned y=0; y<DSPT_TILE; y++)
				{
					for(unsigned x=0; x<DSPT_TILE; x++)
					{
						for(unsigned d=0; d<DSPT_ORNT; d++)
						{
							DATA &dspt = pKp->m_aaaDspt[y][x][d];
							dspt *= invS;
							if(dspt > 0.2F)
							{
								dspt = 0.2F;
								bConti = true;
							}
							else {}
						}
					}
				}
			}
			else {}
		}

		pKp->m_x = pKp->m_x*scl + sOffset;
		pKp->m_y = pKp->m_y*scl + sOffset;
		m_vKeypoint.push_back(*pKp);		
	}
	vKp.clear();
}

//*************************************************************************************************

void Sift::Test()
{
}
