#include "volumeData.h"

/*
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
*/

//*************************************************************************************************
//
//*************************************************************************************************

VolumeFile::VolumeFile(string f)
	:m_vDim(0, 0, 0), m_vStep(0, 0, 0), m_vDimStep(0, 0, 0)
	,m_pVolInt(0), m_pVolScl(0)
{
	string line;
	stringstream ss;

	ifstream inF(f.c_str());

	getline(inF, line);
	ss.clear();
 	ss << line.substr( line.find(":")+1 );
	ss >> m_vDim.m_x >> m_vDim.m_y >> m_vDim.m_z;
	cout << "dimension:\t" << m_vDim.m_x << " " << m_vDim.m_y << " " << m_vDim.m_z << endl;

	//*******************************************

	getline(inF, line);
	ss.clear();
	ss << line.substr( line.find(":")+1 );
	ss >> m_vStep.m_x >> m_vStep.m_y >> m_vStep.m_z;
	cout << "step:\t" << m_vStep.m_x << " " << m_vStep.m_y << " " << m_vStep.m_z << endl;

	//*******************************************

	getline(inF, line);
	ss.clear();
	ss << line.substr( line.find(":")+1 );
	ss >> m_vDimStep.m_x >> m_vDimStep.m_y >> m_vDimStep.m_z;
	cout << "dimStep:\t" << m_vDimStep.m_x << " " << m_vDimStep.m_y << " " << m_vDimStep.m_z << endl;

	//*******************************************

	getline(inF, line);
	ss.clear();
	ss << line.substr( line.find(":")+1 );
	ss >> m_format;
	cout << "format:\t" << m_format << endl;

	//*******************************************

	getline(inF, line);
	ss.clear();
	ss << line.substr( line.find(":")+1 );
	ss >> m_unitSize;
	cout << "unitSize:\t" << m_unitSize << endl;

	//*******************************************

	getline(inF, line);
	ss.clear();
	ss << line.substr( line.find(":")+1 );
	ss >> m_fTplate;
	cout << "template:\t" << m_fTplate << endl;
	inF.close();

	//*******************************************

	LoadIntensity();
	GenSclVolume();
}

VolumeFile::~VolumeFile()
{
	delete m_pVolInt;
	delete m_pVolScl;
}

//***********************************************

short ChangeEndian(short a)
{
	short b = 0;
	b = b | (a&0x00ff);

	a = a >> 8;
	short c = 0;
	c = c | (a&0x00ff);

	b = b<<8;
	b = b + c;

	return b;
}

void VolumeFile::GetFType()
{
	size_t tLoc = m_fTplate.find_first_of(".");
	if(tLoc != m_fTplate.npos)
	{
		m_fType = m_fTplate.substr(tLoc+1);
	}
	else
	{
		m_fType = "stanford";
	}
}

unsigned VolumeFile::GetUSize() const
{
	return m_unitSize;
}
Vect3D<unsigned> VolumeFile::GetDimStep() const
{
	return m_vDimStep;
}

Vect3D<unsigned> VolumeFile::GetOrgDim() const
{
	return m_vDim;
}
Vect3D<unsigned> VolumeFile::GetSclDim() const
{
	return m_pVolScl->GetDim();
}

Vect3D<DATA> VolumeFile::GetOrgStep() const
{
	return m_vStep;
}
Vect3D<DATA> VolumeFile::GetSclStep() const
{
	return Vect3D<DATA>(
		m_vStep.m_x*m_vDimStep.m_x,
		m_vStep.m_y*m_vDimStep.m_y,
		m_vStep.m_z*m_vDimStep.m_z);
}

void VolumeFile::LoadIntensity()
{
	string line;
	stringstream ss;

	delete m_pVolInt;

	Layer lyrInt(m_vDim.m_x, m_vDim.m_y, m_vDim.m_z);
	m_pVolInt = new VolumeData(lyrInt, m_vStep);

	//*******************************************

	string sfName;
	string no;
	ifstream fSlice;

	GetFType();
	if(m_fType.compare("raw") == 0 ||
	   m_fType.compare("p") == 0 ||
	   m_fType.compare("lyr") == 0)
	{
		fSlice.open(m_fTplate.c_str(), ios::binary);
	}
	else {}

	if(m_fType.compare("lyr") == 0)
	{
		unsigned aDim[3];
		fSlice.read((char*)aDim, sizeof(unsigned)*3);
		//cout << aDim[0] << " " << aDim[1] << " " << aDim[2] << endl;
		//getchar();
		//cout << sizeof(DATA) << endl;
	}
	else {}

	for(unsigned z=1; z<=m_vDim.m_z; z++)
	{
		if(m_fType.compare("stanford") == 0)
		{
			ss.clear();
			ss << z;
			ss >> no;
			sfName = m_fTplate + "." + no;
			cout << sfName << endl;
			fSlice.open(sfName.c_str(), ios::binary);
		}
		else {}

		for(unsigned y=0; y<m_vDim.m_y; y++)
		{
			for(unsigned x=0; x<m_vDim.m_x; x++)
			{
				if(m_format.compare("integer") == 0)
				{
					if(m_unitSize == 2)
					{
						unsigned short a;
						fSlice.read((char*)&a, m_unitSize);

						unsigned short b;
						if(m_fType.compare("stanford") == 0 ||
						   m_fType.compare("p") == 0)
						{
							b = ChangeEndian(a);
						}
						else if(m_fType.compare("raw") == 0)
						{
							b = a;
						}
						else
						{
							assert(0);
						}
						lyrInt.CellRef(x, y, z-1) = (DATA)b;
					}
					else if(m_unitSize == 1)
					{
						unsigned char a;
						fSlice.read((char*)&a, m_unitSize);

						assert(m_fType.compare("raw") == 0);
						lyrInt.CellRef(x, y, z-1) = (DATA)a;
					}
					else {}
				}
				else if(m_format.compare("float") == 0)
				{
					if(m_unitSize == 8)
					{	
						double a;
						assert(sizeof(double) == m_unitSize);
						fSlice.read((char*)&a, sizeof(double));
						//cout << a << " ";

						assert(m_fType.compare("lyr") == 0);
						lyrInt.CellRef(x, y, z-1) = (DATA)a;
					}
					else {}
				}
				else
				{
					assert(0);
				}
			} // x
		} // y
		if(m_fType.compare("stanford") == 0)
		{
			fSlice.close();
			//cout << "no"
		}
		else {}
	} // z
	if(m_fType.compare("raw") == 0 ||
	   m_fType.compare("p") == 0)
	{
		fSlice.close();
	}
	else {}
	cout << "read intensity ok" << endl;
}

//***********************************************

const VolumeData& VolumeFile::GetOrgVolVal() const
{
	return *m_pVolInt;
}
const VolumeData& VolumeFile::GetSclVolVal() const
{
	return *m_pVolScl;
}

/*
void VolumeFile::GenSclVolume()
{
	Vect3D<unsigned> dimScl((m_vDim.m_x-1) / m_vDimStep.m_x + 1,
							(m_vDim.m_y-1) / m_vDimStep.m_y + 1,
							(m_vDim.m_z-1) / m_vDimStep.m_z + 1);
	Vect3D<DATA> stepScl(m_vStep.m_x * m_vDimStep.m_x,
						 m_vStep.m_y * m_vDimStep.m_y,
						 m_vStep.m_z * m_vDimStep.m_z);

	delete m_pVolScl;
	m_pVolScl = new VolumeData(dimScl, stepScl);

	Vect3D<unsigned> dimG(7, 7, 7);
	VolumeData volG(dimG, stepScl);
	lyrOp.Gauss3D.Gen(volG, 2.F);

	for(unsigned z=0; z<dimScl.m_z; z++)
	{
		unsigned zOrg = z * m_vDimStep.m_z;
		for(unsigned y=0; y<dimScl.m_y; y++)
		{
			unsigned yOrg = y * m_vDimStep.m_y;
			for(unsigned x=0; x<dimScl.m_x; x++)
			{
				unsigned xOrg = x * m_vDimStep.m_x;
				
				DATA sum = 0;
				for(int zz=0; zz<dimG.m_z; zz++)
				{
					int zLoc = (int)zOrg + zz - dimG.m_z/2;
					if(!m_pVolInt->GetLyrVal().IsCInside(zLoc))	continue;
					for(int yy=0; yy<dimG.m_y; yy++)
					{
						int yLoc = (int)yOrg + yy - dimG.m_y/2;
						if(!m_pVolInt->GetLyrVal().IsYInside(yLoc))	continue;
						for(int xx=0; xx<dimG.m_x; xx++)
						{
							int xLoc = (int)xOrg + xx - dimG.m_x/2;
							if(!m_pVolInt->GetLyrVal().IsXInside(xLoc))	continue;
							
							sum += m_pVolInt->CellVal(xLoc, yLoc, zLoc) * volG.CellVal(xx, yy, zz);
						}
					}
				}
				m_pVolScl->CellRef(x, y, z) = sum;
			}
		}
	}
	cout << "compute scale intensity ok" << endl;
}
*/

void VolumeFile::GenSclVolume()
{
	Vect3D<unsigned> dimScl((m_vDim.m_x-1) / m_vDimStep.m_x + 1,
							(m_vDim.m_y-1) / m_vDimStep.m_y + 1,
							(m_vDim.m_z-1) / m_vDimStep.m_z + 1);
	Vect3D<DATA> stepScl(m_vStep.m_x * m_vDimStep.m_x,
						 m_vStep.m_y * m_vDimStep.m_y,
						 m_vStep.m_z * m_vDimStep.m_z);
	
	delete m_pVolScl;
	m_pVolScl = new VolumeData(dimScl, stepScl);

	Layer lyrScl = m_pVolScl->GetLyrRef();
	for(unsigned z=0; z<dimScl.m_z; z++)
	{
		unsigned zOrg = z * m_vDimStep.m_z;
		for(unsigned y=0; y<dimScl.m_y; y++)
		{
			unsigned yOrg = y * m_vDimStep.m_y;
			for(unsigned x=0; x<dimScl.m_x; x++)
			{
				unsigned xOrg = x * m_vDimStep.m_x;
					
				lyrScl.CellRef(x, y, z) = m_pVolInt->CellVal(xOrg, yOrg, zOrg);
			}
		}
	}
	cout << "compute scale intensity ok" << endl;
}