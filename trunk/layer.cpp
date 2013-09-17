#include "layer.h"

Layer2D::Layer2D(unsigned xDim, unsigned yDim, unsigned cDim)
	:m_xDim(xDim), m_yDim(yDim), m_cDim(cDim)
{
	m_vpMtx.clear();
	for(unsigned cc=0; cc<m_cDim; cc++)
	{
		m_vpMtx.push_back(new Mtx(m_xDim, m_yDim));
	}
}

Layer2D::Layer2D(const Layer2D &layOrg)
	:m_xDim(layOrg.m_xDim), m_yDim(layOrg.m_yDim), m_cDim(layOrg.m_cDim)
{
	m_vpMtx.clear();
	for(unsigned c=0; c<m_cDim; c++)
	{
		m_vpMtx.push_back(new Mtx2D(*layOrg.m_vpMtx[c]));
	}
}

Layer2D::Layer2D(const Layer2D &layOrg, unsigned xOffset, unsigned yOffset, unsigned cOffset, unsigned xDim, unsigned yDim, unsigned cDim)
	:m_xDim(xDim), m_yDim(yDim), m_cDim(cDim)
{
	m_vpMtx.clear();
	for(unsigned c=0; c<m_cDim; c++)
	{
		m_vpMtx.push_back(new Mtx2D(*layOrg.m_vpMtx[c+cOffset], xOffset, yOffset, xDim, yDim));
	}
}

Layer2D::Layer2D(const Mtx2D &mtx)
	:m_xDim(mtx.GetDim().m_x), m_yDim(mtx.GetDim().m_y), m_cDim(1)
{
	m_vpMtx.clear();
	m_vpMtx.push_back(new Mtx2D(mtx));
}

Layer2D::Layer2D(const Mtx2D *ptrMtx, unsigned num)
	:m_xDim(ptrMtx[0].GetDim().m_x), m_yDim(ptrMtx[0].GetDim().m_y), m_cDim(num) 
{
	m_vpMtx.clear();
	for(unsigned i=0; i<num; i++)
	{
		m_vpMtx.push_back(new Mtx2D(ptrMtx[i]));
	}
}

Layer2D::~Layer2D()
{
	for(unsigned c=0; c<m_cDim; c++)
	{
		delete m_vpMtx[c];
	}
	m_vpMtx.clear();
}

//*************************************************************************************************

bool Layer2D::IsXInside(int x) const
{
	if(x>=0 && x<(int)m_xDim)	return true;
	else						return false;
}
bool Layer2D::IsYInside(int y) const
{
	if(y>=0 && y<(int)m_yDim)	return true;
	else						return false;
}
bool Layer2D::IsCInside(int c) const
{
	if(c>=0 && c<(int)m_cDim)	return true;
	else						return false;
}
bool Layer2D::IsInside(int x, int y, int c) const
{
	if(x>=0 && x<(int)m_xDim && 
	   y>=0 && y<(int)m_yDim &&
	   c>=0 && c<(int)m_cDim)		return true;
	else							return false;
}

void Layer2D::ClipBound(int &xL, unsigned &xR, int &yB, unsigned &yT, int &zN, unsigned &zF)
{
	if(xL < 0)			xL = 0;				else {}
	if(xR >= m_xDim)	xR = m_xDim - 1;	else {}
	if(yB < 0)			yB = 0;				else {}
	if(yT >= m_yDim)	yT = m_yDim - 1;	else {}
	if(zN < 0)			zN = 0;				else {}
	if(zF >= m_cDim)	zF = m_cDim - 1;	else {}
}

Vect3D<unsigned> Layer2D::GetDim() const
{
	Vect3D<unsigned> dim(m_xDim, m_yDim, m_cDim);
	return dim;
}

DATA& Layer2D::CellRef(unsigned x, unsigned y, unsigned c)
{
	return m_vpMtx[c]->CellRef(x, y);
}
DATA Layer2D::CellVal(unsigned x, unsigned y, unsigned c) const
{
	return m_vpMtx[c]->CellVal(x, y);
}

void Layer2D::CopyFrom(const Layer2D &lyrFrom)
{
	assert(m_xDim==lyrFrom.m_xDim && m_yDim==lyrFrom.m_yDim && m_cDim==lyrFrom.m_cDim);

	for(unsigned c=0; c<m_cDim; c++)
	{
		m_vpMtx[c]->CopyFrom(*lyrFrom.m_vpMtx[c]);
	}
}
void Layer2D::CopyFrom_zFirst(float *pFrom)
{
	unsigned loc = 0;
	for(unsigned y=0; y<m_yDim; y++)
	{
		for(unsigned x=0; x<m_xDim; x++)
		{
			for(unsigned c=0; c<m_cDim; c++)
			{
				CellRef(x, y, c) = (DATA)pFrom[loc];
				loc++;
			}
		}
	}
}
void Layer2D::CopyTo_zFirst(float *pTo, bool bNormal) const
{
	Vect3D<unsigned> dim = GetDim();
	
	unsigned loc = 0;
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			for(unsigned c=0; c<dim.m_z; c++)
			{
				pTo[loc] = (float)CellVal(x, y, c);
				loc++;
			}
		}
	}

	if(bNormal)
	{
		unsigned size = dim.m_x * dim.m_y * dim.m_z;
		for(unsigned i=0; i<size; i++)
		{
			pTo[i] = (float)(pTo[i] * INV_255);
		}
	}
	else {}
}

void Layer2D::CopyTo_zLast(float *pTo) const
{
	Vect3D<unsigned> dim = GetDim();
	
	unsigned loc = 0;
	for(unsigned z=0; z<dim.m_z; z++)
	{
		for(unsigned y=0; y<dim.m_y; y++)
		{
			for(unsigned x=0; x<dim.m_x; x++)
			{
				pTo[loc] = (float)CellVal(x, y, z);
				loc++;
			}
		}
	}	
}

Mtx* Layer2D::GetMtx(unsigned int c) const
{
	Mtx *pMtx = new Mtx(*m_vpMtx[c]);
	return pMtx;
}

//*************************************************************************************************

/*
Layer2D* Layer2D::GetLayer2D() const
{
	//return this;
	
	Layer2D *pLayer = new Layer2D(*this);
	return pLayer;
}
*/

//*************************************************************************************************
//
//*************************************************************************************************

CplxMtx2D::CplxMtx2D(unsigned xDim, unsigned yDim)
	:m_xDim(xDim), m_yDim(yDim)
{
	m_pLay = new Layer2D(xDim, yDim, 2);
}

CplxMtx2D::CplxMtx2D(const CplxMtx2D &cMtxOrg, unsigned xOffset, unsigned yOffset, unsigned xDim, unsigned yDim)
	:m_xDim(xDim), m_yDim(yDim)
{
	m_pLay = new Layer2D(*cMtxOrg.m_pLay, xOffset, yOffset, 0, xDim, yDim, 2);	
}

CplxMtx2D::~CplxMtx2D()
{
	delete m_pLay;
}

bool CplxMtx2D::IsXInside(int x) const
{
	if(x>=0 && x<(int)m_xDim)	return true;
	else						return false;
}
bool CplxMtx2D::IsYInside(int y) const
{
	if(y>=0 && y<(int)m_yDim)	return true;
	else						return false;
}
bool CplxMtx2D::IsInside(int x, int y) const
{
	if(x>=0 && x<(int)m_xDim && y>=0 && y<(int)m_yDim)	return true;
	else												return false;
}

Vect2D<unsigned> CplxMtx2D::GetDim() const
{
	Vect2D<unsigned> dim(m_xDim, m_yDim);
	return dim;
}

Vect2D<DATA&> CplxMtx2D::CellRef(unsigned x, unsigned y)
{
	Vect2D<DATA&> cCell(m_pLay->CellRef(x, y, 0), m_pLay->CellRef(x, y, 1));
	return cCell;
}

Vect2D<DATA> CplxMtx2D::CellVal(unsigned x, unsigned y) const
{
	Vect2D<DATA> cCell(m_pLay->CellVal(x, y, 0), m_pLay->CellVal(x, y, 1));
	return cCell;
}

Mtx* CplxMtx2D::GetRealMtx() const
{
	Mtx *pMtx = m_pLay->GetMtx(0);
	return pMtx;
}
Mtx* CplxMtx2D::GetImgMtx() const
{
	Mtx *pMtx = m_pLay->GetMtx(1);
	return pMtx;
}

//*************************************************************************************************
//
//*************************************************************************************************

MyImg::MyImg(unsigned xDim, unsigned yDim, unsigned cDim)
	:Layer2D(xDim, yDim, cDim)
{}

MyImg::MyImg(const MyImg &imgOrg)
	:Layer2D(imgOrg)
{}

MyImg::MyImg(const Layer &lyr)
	:Layer2D(lyr)
{}

MyImg::MyImg(const MyImg &imgOrg, unsigned xOffset, unsigned yOffset, unsigned xDim, unsigned yDim)
	:Layer2D(imgOrg, xOffset, yOffset, 0, xDim, yDim, imgOrg.GetDim().m_z)
{}

MyImg::MyImg(const Mtx &mtx)
	:Layer2D(mtx)
{}

MyImg::~MyImg()
{}

//*************************************************************************************************

MyImg* MyImg::ConvertGray() const
{
	if(m_cDim == 1)
	{
		MyImg *pGray = new MyImg(*this);
		return pGray;
	}
	else if(m_cDim == 3 || m_cDim == 4)
	{
		MyImg *pGray = new MyImg(m_xDim, m_yDim, 1);
		for(unsigned y=0; y<m_yDim; y++)
		{
			for(unsigned x=0; x<m_xDim; x++)
			{
				DATA r = CellVal(x, y, 0);
				DATA g = CellVal(x, y, 1);
				DATA b = CellVal(x, y, 2);
				pGray->CellRef(x, y, 0) = 0.3F*r + 0.59F*g + 0.11F*b;
			}
		}
		return pGray;
	}
	else
	{
		assert(0);
	}
	return 0;
}

//*************************************************************************************************
//
//*************************************************************************************************

VolumeData::VolumeData(Layer &lyrData, Vect3D<DATA> step)
	:m_lyrData(lyrData), m_vStep(step)
{}

VolumeData::VolumeData(Vect3D<unsigned> dim, Vect3D<DATA> step)
	:m_lyrData(dim.m_x, dim.m_y, dim.m_z), m_vStep(step)
{}

VolumeData::~VolumeData()
{}

//*************************************************************************************************

Vect3D<unsigned> VolumeData::GetDim() const
{
	return m_lyrData.GetDim();
}

Vect3D<DATA> VolumeData::GetStep() const
{
	return m_vStep;
}

DATA& VolumeData::CellRef(unsigned x, unsigned y, unsigned c)
{
	return m_lyrData.CellRef(x, y, c);
}

DATA VolumeData::CellVal(unsigned x, unsigned y, unsigned c) const
{
	return m_lyrData.CellVal(x, y, c);
}

Layer& VolumeData::GetLyrRef()
{
	return m_lyrData;
}
const Layer& VolumeData::GetLyrVal() const
{
	return m_lyrData;
}