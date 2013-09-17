#ifndef _SIFT_H
#define _SIFT_H

#include "denseMatrix.h"
#include "matrixOperation.h"
#include "gslMtxOp.h"
#include "layerOperation.h"

#include <vector>
#include <cstdio>
#include <fstream>

using namespace std;

class Octave
{
friend class SmpWin;

public:
	Octave(Mtx &mtxInit, unsigned GaussNum, unsigned octNum, DATA dev);
	~Octave();

	Mtx*	GetBlur(unsigned octNo, unsigned GaussNo);
	Mtx*	GetDoG(unsigned octNo, unsigned DoGNo);
	Layer*	GetOrnt(unsigned octNo, unsigned GaussNo);
	unsigned GetGaussNum();

private:
	const unsigned m_octNum;
	const unsigned m_sclPerOct;
	const unsigned m_GaussNum;
	const unsigned m_GaussLen;

	DATA	*m_pDev;
	Mtx		**m_ppGauss;
	Mtx		***m_pppBlur;
	Mtx		***m_pppDoG;
	Layer	***m_pppOrnt;
};

class SmpWin
{
public:
	SmpWin(Octave &oct);
	~SmpWin();

	Mtx* GetGauss(unsigned GaussNo);

private:
	const unsigned m_GaussNum;
	const unsigned m_GaussLen;

	DATA *m_pDev;
	Mtx **m_ppGauss;
};

class Keypoint
{
public:
	Keypoint();
	~Keypoint();

	unsigned m_octNo;
	unsigned m_GaussNo;

	DATA m_x;		
	DATA m_y;
	DATA m_ornt;
	DATA m_aaaDspt[4][4][8];

private:
};

//*************************************************************************************

class Sift
{
public:
	Sift();
	~Sift();

	vector<Keypoint> Gen(Mtx &mtx);

	static void Test();

private:
	unsigned ComputeOctNum(const Mtx &mtx) const;

	void GetExtrema(vector<Vect2D<unsigned>> &extre);
	void NormalizeCell();

	void RefineKeypoint(vector<Vect3D<DATA>> &vKeyLoc, vector<Vect2D<unsigned>> &vExtre);
	bool IsLowContrast(DATA smpVal, const Mtx &mtxDer, const Mtx &mtxOpt);
	bool IsEdge(Mtx &mtx, unsigned x, unsigned y);

	DATA LDir(DATA *pDir, unsigned now);
	DATA RDir(DATA *pDir, unsigned now);
	void AssignOrient(vector<Keypoint> &vKp, vector<Vect3D<DATA>> &vKeyLoc);
	DATA SolveOrnt(DATA lVal, DATA nVal, DATA rVal, unsigned idx);

	void DescriptNbr(vector<Keypoint> &vKp);
	void ZeroDspt(Keypoint &kp);

	const unsigned SCALE_PER_OCTAVE;
	const DATA DEV;
	const DATA CONTRAST_THRHLD;
	const DATA EDGE_RATIO;
	const DATA DIR_RATIO;
	const unsigned ORNT_BIN;
	const unsigned DSPT_TILE;
	const unsigned DSPT_ORNT;
	const unsigned MIN_X;
	const unsigned MIN_Y;

	Mtx *m_pMtxIn;
	Octave *m_pOct;
	SmpWin *m_pOrntSwin;
	Mtx m_invOrnt;

	DATA m_contrastThrhld;
	unsigned m_octNum;
	unsigned m_octNo;
	unsigned m_DoGNo;
	DATA m_orntRange;
	DATA m_dspORange;
	DATA *m_pDir;
	vector<Keypoint> m_vKeypoint;
};

#endif