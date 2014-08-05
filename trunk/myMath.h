#ifndef _MY_MATH_H
#define _MY_MATH_H

#include "define.h"
#include <cassert>
#include <cstdlib>
#include <ctime>

class MyMath
{
public:
    MyMath();
    ~MyMath();

    bool IsEqual (const DATA &a, const DATA &b, const DATA &delta = (DATA)1e-10) const;
    bool IsEGreat(const DATA &a, const DATA &b, const DATA &delta = (DATA)1e-10) const;
    bool IsELess (const DATA &a, const DATA &b, const DATA &delta = (DATA)1e-10) const;
    bool IsGreat (const DATA &a, const DATA &b, const DATA &delta = (DATA)1e-10) const;
    bool IsLess  (const DATA &a, const DATA &b, const DATA &delta = (DATA)1e-10) const;

    void CopyNV(DATA aTo[], const DATA aFrom[], unsigned n) const;
    void Copy3V(DATA aTo[], const DATA aFrom[]) const;
    void Copy4V(DATA aTo[], const DATA aFrom[]) const;
    void Copy4M(DATA aaTo[][4], const DATA aaFrom[][4]) const;

    void AddNV(DATA aA[], const DATA aB[], unsigned n) const;
    void Add3V(DATA aA[], const DATA aB[]) const;
    void SubNV(DATA aA[], const DATA aB[], unsigned n) const;
    void Sub3V(DATA aA[], const DATA aB[]) const;

    void MulNV(DATA aA[], DATA m, unsigned n) const;
    void Mul3V(DATA aA[], DATA m) const;
    void Mul3MV(DATA aR[], const DATA aaMtx[][3], const DATA aV[]) const;
    void Mul4MV(DATA aR[], const DATA aaMtx[][4], const DATA aV[]) const;
    void Mul3MM(DATA aaMR[][3], const DATA aaMA[][3], const DATA aaMB[][3]) const;
    void Mul4MM(DATA aaMR[][4], const DATA aaMA[][4], const DATA aaMB[][4]) const;

    DATA Dot3D(const DATA aA[], const DATA aB[]) const;
    void Cross3V(DATA aR[], const DATA aA[], const DATA aB[]) const;
    DATA LenNV(const DATA aV[], unsigned n) const;
    DATA Len2V(const DATA aV[]) const;
    DATA Len3V(const DATA aV[]) const;
    void NormalNV(DATA aV[], unsigned n) const;
    void Normal2V(DATA aV[]) const;
    void Normal3V(DATA aV[]) const;

    void GetRotXMtx(DATA aaMtx[][4], DATA ang) const;
    void GetRotYMtx(DATA aaMtx[][4], DATA ang) const;
    void GetRotZMtx(DATA aaMtx[][4], DATA ang) const;
    void GetRotVecMtx(DATA faaMtx[][4], DATA angle, const DATA aVect[], bool bNorm = false) const;

    void RotateX(DATA aOut[], const DATA aIn[], DATA ang) const;
    void RotateY(DATA aOut[], const DATA aIn[], DATA ang) const;
    void RotateZ(DATA aOut[], const DATA aIn[], DATA ang) const;
    void RotateVec(DATA aOut[], const DATA aIn[], DATA ang, const DATA aVec[], bool bNorm = false) const;

    void RGB2HSV(DATA aHsv[], const DATA aRgb[], bool bNormH) const;
    void HSV2RGB(DATA aRgb[], const DATA aHsv[], bool bNormH) const;
    DATA RGB2Gray(const DATA aRgb[]) const;
    DATA RGB2Gray2(const DATA aRgb[]) const;

    DATA Interpolate_linear(DATA a, DATA b, DATA ratioA) const;
    DATA Interpolate_linear(DATA lb, DATA rb, DATA lt, DATA rt, DATA ratioL, DATA ratioB) const;
    DATA Interpolate_linear(DATA lbn, DATA rbn, DATA ltn, DATA rtn, DATA lbf, DATA rbf, DATA ltf, DATA rtf, DATA ratioL, DATA ratioB, DATA ratioN) const;

    bool IntersectLinePlane3D(DATA aPInt[], const DATA aOri[], const DATA aV[], const DATA aPlane[], DATA tMin, DATA tMax) const;
    bool IntersectLinePlane3D(DATA aPInt[], const DATA aPA[], const DATA aPB[], const DATA aPlane[]) const;

    DATA ComputeLenRatio(DATA val, DATA min, DATA len);

    void Range(DATA &min, DATA &max, DATA aV[], unsigned size);
    void Range(float &min, float &max, float aV[], unsigned size);

	//DATA Max3(DATA v0, DATA v1, DATA v2);
	//DATA MIn3(DATA v0, DATA v1, DATA v2);

    DATA Rnd();

private:
};

extern MyMath myMath;

#endif