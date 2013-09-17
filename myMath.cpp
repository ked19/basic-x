#include "myMath.h"
#include <iostream>

MyMath myMath;

MyMath::MyMath()
{}

MyMath::~MyMath()
{}

//***********************************************

void MyMath::CopyNV(DATA aTo[], const DATA aFrom[], unsigned n) const
{
	for(unsigned i=0; i<n; i++)
	{
		aTo[i] = aFrom[i];
	}
}
void MyMath::Copy3V(DATA aTo[], const DATA aFrom[]) const
{
	CopyNV(aTo, aFrom, 3);
}
void MyMath::Copy4V(DATA aTo[], const DATA aFrom[]) const
{
	CopyNV(aTo, aFrom, 4);
}

#define CopyNM(aaTo, aaFrom, n)							\
{														\
    for(unsigned j=0; j<(n); j++)						\
    {													\
        for(unsigned i=0; i<(n); i++)					\
        {												\
            (aaTo)[j][i] = (aaFrom)[j][i];				\
        }												\
    }													\
}												
void MyMath::Copy4M(DATA aaTo[][4], const DATA aaFrom[][4]) const
{
    CopyNM(aaTo, aaFrom, 4);
}

//***********************************************

void MyMath::AddNV(DATA aA[], const DATA aB[], unsigned n) const
{
	for(unsigned i=0; i<n; i++)
	{
		aA[i] += aB[i];
	}
}
void MyMath::Add3V(DATA aA[], const DATA aB[]) const
{
	AddNV(aA, aB, 3);
}

void MyMath::SubNV(DATA aA[], const DATA aB[], unsigned n) const
{
	for(unsigned i=0; i<n; i++)
	{
		aA[i] -= aB[i];
	}
}
void MyMath::Sub3V(DATA aA[], const DATA aB[]) const
{
	SubNV(aA, aB, 3);
}

//***********************************************

void MyMath::MulNV(DATA aA[], DATA m, unsigned n) const
{
	for(unsigned i=0; i<n; i++)
	{
		aA[i] *= m;
	}
}
void MyMath::Mul3V(DATA aA[], DATA m) const
{
	MulNV(aA, m, 3);
}

#define MulNMV(aR, aaMtx, aV, n)						\
{														\
    for(unsigned j=0; j<(n); j++)						\
    {													\
        (aR)[j] = 0;									\
        for(unsigned i=0; i<(n); i++)					\
        {												\
            (aR)[j] += (aaMtx)[j][i] * (aV)[i];			\
        }												\
    }													\
}
void MyMath::Mul3MV(DATA aR[], const DATA aaMtx[][3], const DATA aV[]) const
{
	MulNMV(aR, aaMtx, aV, 3);
}
void MyMath::Mul4MV(DATA aR[], const DATA aaMtx[][4], const DATA aV[]) const
{
	MulNMV(aR, aaMtx, aV, 4);
}

#define MulMM(aaMR, aaMA, aaMB, n)									\
{																	\
	for(unsigned y=0; y<(n); y++)									\
	{																\
		for(unsigned x=0; x<(n); x++)								\
		{															\
			(aaMR)[y][x] = 0;										\
			for(unsigned v=0; v<(n); v++)							\
			{														\
				(aaMR)[y][x] += (aaMA)[y][v] * (aaMB)[v][x];		\
			}														\
		}															\
	}																\
}
void MyMath::Mul4MM(DATA aaR[][4], const DATA aaMA[][4], const DATA aaMB[][4]) const
{
	MulMM(aaR, aaMA, aaMB, 4);
}
void MyMath::Mul3MM(DATA aaR[][3], const DATA aaMA[][3], const DATA aaMB[][3]) const
{
	MulMM(aaR, aaMA, aaMB, 3);
}

//*************************************************************************************************

void MyMath::GetRotXMtx(DATA aaMtx[][4], DATA ang) const
{
    DATA cosA = cos(ang);
    DATA sinA = sin(ang);
    DATA aaRM[][4] = {{1.f,    0,     0,   0},
                      {  0, cosA, -sinA,   0},
                      {  0, sinA,  cosA,   0},
                      {  0,    0,     0, 1.f}};
    Copy4M(aaMtx, aaRM);
}
void MyMath::RotateX(DATA aOut[], const DATA aIn[], DATA ang) const
{
    DATA aaMtx[4][4];
    GetRotXMtx(aaMtx, ang);
    Mul4MV(aOut, aaMtx, aIn);
}

void MyMath::GetRotYMtx(DATA aaMtx[][4], DATA ang) const
{
    DATA cosA = cos(ang);
    DATA sinA = sin(ang);
    DATA aaRM[][4] = {{ cosA,   0, sinA,   0},
                      {    0, 1.f,    0,   0},
                      {-sinA,   0, cosA,   0},
                      {    0,   0,    0, 1.f}};
    Copy4M(aaMtx, aaRM);
}
void MyMath::RotateY(DATA aOut[], const DATA aIn[], DATA ang) const
{
    DATA aaMtx[4][4];
    GetRotYMtx(aaMtx, ang);
    Mul4MV(aOut, aaMtx, aIn);
}

void MyMath::GetRotZMtx(DATA aaMtx[][4], DATA ang) const
{
    DATA cosA = cos(ang);
    DATA sinA = sin(ang);
    DATA aaRM[][4] = {{ cosA, -sinA,   0,   0},
                      { sinA,  cosA,   0,   0},
                      {    0,     0, 1.f,   0},
                      {    0,     0,   0, 1.f}};
    Copy4M(aaMtx, aaRM);
}
void MyMath::RotateZ(DATA aOut[], const DATA aIn[], DATA ang) const
{
    DATA aaMtx[4][4];
    GetRotZMtx(aaMtx, ang);
    Mul4MV(aOut, aaMtx, aIn);
}

void MyMath::GetRotVecMtx(DATA aaMtx[][4], DATA ang, const DATA aVec[], bool bNorm) const
{
    DATA aNormV[3];
	Copy3V(aNormV, aVec);
	if(bNorm)
	{
		Normal3V(aNormV);
	}
	else {}

    DATA cosA = cos(ang);
    DATA sinA = sin(ang);
    DATA aaTable[3][3] =
       {{            cosA,  aNormV[2]*sinA,  -aNormV[1]*sinA},
         {-aNormV[2]*sinA,            cosA,   aNormV[0]*sinA},
         { aNormV[1]*sinA, -aNormV[0]*sinA,            cosA}};
		/*
        {{           cosA, -aNormV[2]*sinA,  aNormV[1]*sinA},
         {-aNormV[2]*sinA,            cosA, -aNormV[0]*sinA},
         {-aNormV[1]*sinA,  aNormV[0]*sinA,            cosA}};
		 */
    
	DATA aaRM[4][4] = {0};
    for(unsigned j=0; j<3; j++)
    {
        for(unsigned i=0; i<3; i++)
        {
            aaRM[j][i] = aNormV[j]*aNormV[i]*(1.F-cosA) + aaTable[j][i];
        }
    }
    aaRM[3][3] = 1.F;

    Copy4M(aaMtx, aaRM);
}
void MyMath::RotateVec(DATA aOut[], const DATA aIn[], DATA ang, const DATA aVec[], bool bNorm) const
{
    DATA aaMtx[4][4];
    GetRotVecMtx(aaMtx, ang, aVec, bNorm);
    Mul4MV(aOut, aaMtx, aIn);
}

//*************************************************************************************************

void MyMath::RGB2HSV(DATA aHsv[], const DATA aRgb[], bool bNormH) const
{
	DATA cMax, cMin;
	if(aRgb[0] >= aRgb[1])
	{
		cMax = aRgb[0];
		cMin = aRgb[1];
	}
	else
	{
		cMax = aRgb[1];
		cMin = aRgb[0];
	}
	if(aRgb[2] > cMax)
	{
		cMax = aRgb[2];
	}
	else if(aRgb[2] < cMin)
	{
		cMin = aRgb[2];
	}

        
	DATA dist = cMax - cMin;
	if(IsEqual(cMax, cMin))
	{
		aHsv[0] = 0;
	}
	else
	{
		DATA invDist = 1.F / dist;
		if(cMax == aRgb[0])
		{
			if(aRgb[1] >= aRgb[2])
			{
				aHsv[0] = 60.F * (aRgb[1]-aRgb[2]) * invDist;
			}
			else
			{
				aHsv[0] = 60.F * (aRgb[1]-aRgb[2]) * invDist + 360.F;
			}
		}
		else if(cMax == aRgb[1])
		{
			aHsv[0] = 60.F * (aRgb[2]-aRgb[0]) * invDist + 120.F;
		}
		else if(cMax == aRgb[2])
		{
			aHsv[0] = 60.F * (aRgb[0]-aRgb[1]) * invDist + 240.F;
		}
		else
		{
			assert(0);
		}
	}
	if(bNormH)
	{
		aHsv[0] *= INV_360;
	}
	else {}

	aHsv[2] = cMax;
	aHsv[1] = (IsEqual(cMax, 0))? 0: dist/cMax;
}

void MyMath::HSV2RGB(DATA aRgb[], const DATA aHsv[], bool bNormH) const
{
	if(aHsv[0]==0 && aHsv[1]==0 && aHsv[2]==0)
	{
		aRgb[0] = 0;	aRgb[1] = 0;	aRgb[2] = 0;
		return;
	}
	else {}

	DATA c = aHsv[2] * aHsv[1];
	
	DATA h = (bNormH)? aHsv[0]*360.F: aHsv[0];
	DATA hp = h / 60.F;
	
	DATA hp2 = hp;
	while(hp2 > 2.F)
	{
		hp2 -= 2.F;
	}
	DATA x = c * (1.F - fabs(hp2 - 1.F));

	if(hp < 1.F)
	{
		aRgb[0] = c;	aRgb[1] = x;	aRgb[2] = 0;
		//cout << aRgb[0] << endl;
	}
	else if(hp < 2)
	{
		aRgb[0] = x;	aRgb[1] = c;	aRgb[2] = 0;
	}
	else if(hp < 3)
	{
		aRgb[0] = 0;	aRgb[1] = c;	aRgb[2] = x;
	}
	else if(hp < 4)
	{
		aRgb[0] = 0;	aRgb[1] = x;	aRgb[2] = c;
	}
	else if(hp < 5)
	{
		aRgb[0] = x;	aRgb[1] = 0;	aRgb[2] = c;
	}
	else
	{
		aRgb[0] = c;	aRgb[1] = 0;	aRgb[2] = x;
	}

	DATA m = aHsv[2] - c;
	aRgb[0] += m;
	aRgb[1] += m;
	aRgb[2] += m;
}

DATA MyMath::RGB2Gray(const DATA aRgb[]) const
{
	return 0.3F*aRgb[0] + 0.59F*aRgb[1] + 0.11F*aRgb[2];
}

//*************************************************************************************************

DATA MyMath::Dot3D(const DATA aA[], const DATA aB[]) const
{
	DATA r = aA[0]*aB[0] + aA[1]*aB[1] + aA[2]*aB[2];
	return r;
}

void MyMath::Cross3V(DATA aR[], const DATA aA[], const DATA aB[]) const
{
	aR[0] = aA[1]*aB[2] - aA[2]*aB[1];
	aR[1] = aA[2]*aB[0] - aA[0]*aB[2];
	aR[2] = aA[0]*aB[1] - aA[1]*aB[0];
}

DATA MyMath::LenNV(const DATA aV[], unsigned n) const
{
	DATA len = 0;
	for(unsigned i=0; i<n; i++)
	{
		len += aV[i]*aV[i];
	}
	return sqrt(len);
}
DATA MyMath::Len2V(const DATA aV[]) const
{
	return LenNV(aV, 2);
}
DATA MyMath::Len3V(const DATA aV[]) const
{
	return LenNV(aV, 3);
}

void MyMath::NormalNV(DATA aV[], unsigned n) const
{
	DATA len = LenNV(aV, n);
	if(IsGreat(len, 0))
	{
		for(unsigned i=0; i<n; i++)
		{
			aV[i] /= len;
		}
	}
	else {}
}

void MyMath::Normal2V(DATA aV[]) const
{
	NormalNV(aV, 2);
}

void MyMath::Normal3V(DATA aV[]) const
{
	NormalNV(aV, 3);
}

bool MyMath::IsEqual(const DATA &a, const DATA & b, const DATA &delta) const
{
	DATA diff = fabs(a) - fabs(b);
	if(diff <=  delta &&
	   diff >= -delta)		return true;
	else					return false;
}
bool MyMath::IsEGreat(const DATA &a, const DATA &b, const DATA &delta) const
{
	if(a+delta >= b)		return true;
	else					return false;
}
bool MyMath::IsELess(const DATA &a, const DATA &b, const DATA &delta) const
{
	if(a-delta <= b)		return true;
	else					return false;
}
bool MyMath::IsGreat(const DATA &a, const DATA &b, const DATA &delta) const
{
	if(a-delta > b)			return true;
	else					return false;
}
bool MyMath::IsLess(const DATA &a, const DATA &b, const DATA &delta) const
{
	if(a+delta < b)			return true;
	else					return false;
}

DATA MyMath::Interpolate_linear(DATA &a, DATA &b, DATA ratioA) const
{
	DATA ratioB = 1.F - ratioA;
	DATA r = a*ratioB + b*ratioA;
	return r;
}

DATA MyMath::Interpolate_linear(DATA &lb, DATA &rb, DATA &lt, DATA &rt, DATA &ratioL, DATA &ratioB) const
{
	DATA b = Interpolate_linear(lb, rb, ratioL);
	DATA t = Interpolate_linear(lt, rt, ratioL);

	DATA val = Interpolate_linear(b, t, ratioB);
	return val;
}

DATA MyMath::Interpolate_linear(DATA &lbn, DATA &rbn, DATA &ltn, DATA &rtn, DATA &lbf, DATA &rbf, DATA &ltf, DATA &rtf, DATA &ratioL, DATA &ratioB, DATA &ratioN) const
{
	DATA n = Interpolate_linear(lbn, rbn, ltn, rtn, ratioL, ratioB);
	DATA f = Interpolate_linear(lbf, rbf, ltf, rtf, ratioL, ratioB);

	DATA val = Interpolate_linear(n, f, ratioN);
	return val;
}

bool MyMath::IntersectLinePlane3D(DATA aPInt[], const DATA aOri[], const DATA aV[], const DATA aPlane[], DATA tMin, DATA tMax) const
{
	DATA t = -(Dot3D(aPlane, aOri) + aPlane[3]) / Dot3D(aPlane, aV);
	if((tMin!=DATA_NAN && t<tMin) ||
	   (tMax!=DATA_NAN && t>tMax))
	{
		return false;
	}
	else {}

	for(unsigned i=0; i<3; i++)
	{
		aPInt[i] = aOri[i] + t*aV[i];
	}
	return true;
}
bool MyMath::IntersectLinePlane3D(DATA aPInt[], const DATA aPA[], const DATA aPB[], const DATA aPlane[]) const
{
	DATA aV[3];
	for(unsigned i=0; i<3; i++)
	{
		aV[i] = aPB[i] - aPA[i];
	}
	return IntersectLinePlane3D(aPInt, aPA, aV, aPlane, 0, 1.F);
}

DATA MyMath::ComputeLenRatio(DATA val, DATA min, DATA len)
{
	return (val - min) / len;
}

void MyMath::Range(DATA &min, DATA &max, DATA aV[], unsigned size)
{
	min = 1e10;
	max = -1e10;
	for(unsigned i=0; i<size; i++)
	{
		if(min > aV[i])	min = aV[i];	else {}
		if(max < aV[i])	max = aV[i];	else {}
	}
}
void MyMath::Range(float &min, float &max, float aV[], unsigned size)
{
	min = 1e10;
	max = -1e10;
	for(unsigned i=0; i<size; i++)
	{
		if(min > aV[i])	min = aV[i];	else {}
		if(max < aV[i])	max = aV[i];	else {}
	}
}