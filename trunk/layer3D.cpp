#include "layer3D.h"

Layer3D::Layer3D(unsigned xDim, unsigned yDim, unsigned zDim, unsigned cDim)
	:m_xDim(xDim), m_yDim(yDim), m_zDim(zDim), m_cDim(cDim)
{
	m_vpLyr.clear();
	for(unsigned cc=0; cc<m_cDim; cc++)
	{
		m_vpLyr.push_back(new Layer2D(m_xDim, m_yDim, m_zDim));
	}
}

Layer3D::Layer3D(const Layer3D &lay3Org)
	:m_xDim(lay3Org.m_xDim), m_yDim(lay3Org.m_yDim), m_zDim(lay3Org.m_zDim), m_cDim(lay3Org.m_cDim)
{
	m_vpLyr.clear();
	for(unsigned c=0; c<m_cDim; c++)
	{
		m_vpLyr.push_back(new Layer2D(*lay3Org.m_vpLyr[c]));
	}
}

Layer3D::Layer3D(const Layer3D &lay3Org, unsigned xOffset, unsigned yOffset, unsigned zOffset, unsigned cOffset, 
				 unsigned xDim, unsigned yDim, unsigned zDim, unsigned cDim)
	:m_xDim(xDim), m_yDim(yDim), m_zDim(zDim), m_cDim(cDim)
{
	m_vpLyr.clear();
	for(unsigned c=0; c<m_cDim; c++)
	{
		m_vpLyr.push_back(new Layer2D(*lay3Org.m_vpLyr[c+cOffset], xOffset, yOffset, zOffset, xDim, yDim, zDim));
	}
}

Layer3D::~Layer3D()
{
	for(unsigned c=0; c<m_cDim; c++)
	{
		delete m_vpLyr[c];
	}
	m_vpLyr.clear();
}

//*************************************************************************************************

bool Layer3D::IsXInside(int x) const
{
	if(x>=0 && x<(int)m_xDim)	return true;
	else						return false;
}
bool Layer3D::IsYInside(int y) const
{
	if(y>=0 && y<(int)m_yDim)	return true;
	else						return false;
}
bool Layer3D::IsZInside(int z) const
{
	if(z>=0 && z<(int)m_zDim)	return true;
	else						return false;
}
bool Layer3D::IsCInside(int c) const
{
	if(c>=0 && c<(int)m_cDim)	return true;
	else						return false;
}
bool Layer3D::IsInside(int x, int y, int z, int c) const
{
	if(x>=0 && x<(int)m_xDim && 
	   y>=0 && y<(int)m_yDim &&
	   z>=0 && z<(int)m_zDim &&
	   c>=0 && c<(int)m_cDim)		return true;
	else							return false;
}

Vect4D<unsigned> Layer3D::GetDim() const
{
	Vect4D<unsigned> dim(m_xDim, m_yDim, m_zDim, m_cDim);
	return dim;
}

DATA& Layer3D::CellRef(unsigned x, unsigned y, unsigned z, unsigned c)
{
	return m_vpLyr[c]->CellRef(x, y, z);
}
DATA Layer3D::CellVal(unsigned x, unsigned y, unsigned z, unsigned c) const
{
	return m_vpLyr[c]->CellVal(x, y, z);
}

void Layer3D::CopyFrom(const Layer3D &lyr3From)
{
	assert(m_xDim==lyr3From.m_xDim && m_yDim==lyr3From.m_yDim && 
		   m_zDim==lyr3From.m_zDim && m_cDim==lyr3From.m_cDim);

	for(unsigned c=0; c<m_cDim; c++)
	{
		m_vpLyr[c]->CopyFrom(*lyr3From.m_vpLyr[c]);
	}
}
/*
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
*/

Layer* Layer3D::GetLyr(unsigned int c) const
{
	Layer *pLyr = new Layer(*m_vpLyr[c]);
	return pLyr;
}

//*************************************************************************************************
//
//*************************************************************************************************

CplxLyr2D::CplxLyr2D(unsigned xDim, unsigned yDim, unsigned zDim)
	:m_xDim(xDim), m_yDim(yDim), m_zDim(zDim)
{
	m_pLyr3 = new Layer3D(xDim, yDim, zDim, 2);
}

CplxLyr2D::CplxLyr2D(const CplxLyr2D &cLyrOrg, unsigned xOffset, unsigned yOffset, unsigned zOffset, unsigned xDim, unsigned yDim, unsigned zDim)
	:m_xDim(xDim), m_yDim(yDim), m_zDim(zDim)
{
	m_pLyr3 = new Layer3D(*cLyrOrg.m_pLyr3, xOffset, yOffset, zOffset, 0, xDim, yDim, zDim, 2);	
}

CplxLyr2D::~CplxLyr2D()
{
	delete m_pLyr3;
}

bool CplxLyr2D::IsXInside(int x) const
{
	if(x>=0 && x<(int)m_xDim)	return true;
	else						return false;
}
bool CplxLyr2D::IsYInside(int y) const
{
	if(y>=0 && y<(int)m_yDim)	return true;
	else						return false;
}
bool CplxLyr2D::IsZInside(int z) const
{
	if(z>=0 && z<(int)m_zDim)	return true;
	else						return false;
}
bool CplxLyr2D::IsInside(int x, int y, int z) const
{
	if(x>=0 && x<(int)m_xDim && 
	   y>=0 && y<(int)m_yDim && 
	   z>=0 && z<(int)m_zDim)	return true;
	else						return false;
}

Vect3D<unsigned> CplxLyr2D::GetDim() const
{
	Vect3D<unsigned> dim(m_xDim, m_yDim, m_zDim);
	return dim;
}

Vect2D<DATA&> CplxLyr2D::CellRef(unsigned x, unsigned y, unsigned z)
{
	Vect2D<DATA&> cCell(m_pLyr3->CellRef(x, y, z, 0), m_pLyr3->CellRef(x, y, z, 1));
	return cCell;
}

Vect2D<DATA> CplxLyr2D::CellVal(unsigned x, unsigned y, unsigned z) const
{
	Vect2D<DATA> cCell(m_pLyr3->CellVal(x, y, z, 0), m_pLyr3->CellVal(x, y, z, 1));
	return cCell;
}

Layer* CplxLyr2D::GetRealLyr() const
{
	Layer *pLyr = m_pLyr3->GetLyr(0);
	return pLyr;
}
Layer* CplxLyr2D::GetImgLyr() const
{
	Layer *pLyr = m_pLyr3->GetLyr(1);
	return pLyr;
}

//*************************************************************************************************
//
//*************************************************************************************************

VectorData::VectorData(Layer3D &l3dData, Vect3D<DATA> step)
	:m_l3dData(l3dData), m_vStep(step)
{}

VectorData::VectorData(Vect3D<unsigned> dim, Vect3D<DATA> step)
	:m_l3dData(dim.m_x, dim.m_y, dim.m_z, 3), m_vStep(step)
{}

VectorData::~VectorData()
{}

Vect3D<unsigned> VectorData::GetDim() const
{
	Vect4D<unsigned> dim = m_l3dData.GetDim();
	Vect3D<unsigned> dimOut(dim.m_r, dim.m_g, dim.m_b);
	return dimOut;
}

Vect3D<DATA> VectorData::GetStep() const
{
	return m_vStep;
}

DATA& VectorData::CellRef(unsigned x, unsigned y, unsigned z, unsigned c)
{
	return m_l3dData.CellRef(x, y, z, c);
}

DATA VectorData::CellVal(unsigned x, unsigned y, unsigned z, unsigned c) const
{
	return m_l3dData.CellVal(x, y, z, c);
}

Layer3D& VectorData::GetL3dRef()
{
	return m_l3dData;
}

const Layer3D& VectorData::GetL3dVal() const
{
	return m_l3dData;
}