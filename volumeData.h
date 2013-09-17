#ifndef _VOLUME_DATA_H
#define _VOLUME_DATA_H

#include "layer.h"
#include "layerOperation.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/*
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
*/

class VolumeFile
{
public:
	VolumeFile(string f);
	~VolumeFile();

	const VolumeData& GetSclVolVal() const;
	const VolumeData& GetOrgVolVal() const;
	
	unsigned GetUSize() const;
	Vect3D<unsigned> GetDimStep() const;

	Vect3D<unsigned> GetOrgDim() const;
	Vect3D<unsigned> GetSclDim() const;
	Vect3D<DATA> GetOrgStep() const;
	Vect3D<DATA> GetSclStep() const;

private:
	void LoadIntensity();
	void GenSclVolume();
	void GetFType();

	Vect3D<unsigned> m_vDim;
	Vect3D<DATA> m_vStep;
	Vect3D<unsigned> m_vDimStep;

	string m_format;
	unsigned m_unitSize;
	string m_fTplate;
	string m_fType;

	VolumeData *m_pVolInt;
	VolumeData *m_pVolScl;
};

#endif