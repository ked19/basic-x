#include "myDraw.h"

MyDraw myDraw;

LineDraw::LineDraw()
{}

LineDraw::~LineDraw()
{}

void LineDraw::Gen(Layer &lyrD, Vect2D<DATA> &pFrom, Vect2D<DATA> &pTo, 
							Vect4D<DATA> &cFrom, Vect4D<DATA> &cTo, bool bDebug)
{
	if (bDebug) {
		cout << "line" << endl;
	} else {}

	MyAssert(cFrom.m_r <= 1.F &&
			 cFrom.m_g <= 1.F &&
			 cFrom.m_b <= 1.F &&
			 cFrom.m_a <= 1.F);
	MyAssert(cTo.m_r <= 1.F &&
			 cTo.m_g <= 1.F &&
			 cTo.m_b <= 1.F &&
			 cTo.m_a <= 1.F);

	Vect3D<unsigned> dimD = lyrD.GetDim();
	MyAssert(dimD.m_z == 4);

	DATA dX = pTo.m_x - pFrom.m_x;
	DATA dY = pTo.m_y - pFrom.m_y;

	bool bXMajor = (fabs(dX) >= fabs(dY)) ? true : false;
	if (bDebug) {
		cout << "bXMajor: " << bXMajor << endl;
	} else {}

	int fFrom = (bXMajor) ? 
		(int)(pFrom.m_x + 0.5F) : (int)(pFrom.m_y + 0.5F);
	int fTo = (bXMajor) ?
		(int)(pTo.m_x + 0.5F) : (int)(pTo.m_y + 0.5F);
	int fStp = (fFrom <= fTo) ?
		1 : -1;

	DATA sFrom = (bXMajor) ?
		pFrom.m_y : pFrom.m_x;
	DATA sDis = (bXMajor) ?
		pTo.m_y - pFrom.m_y : pTo.m_x - pFrom.m_x;
	Vect4D<DATA> cDis(cTo.m_r - cFrom.m_r,
					  cTo.m_g - cFrom.m_g,
					  cTo.m_b - cFrom.m_b,
					  cTo.m_a - cFrom.m_a);
	for (int f = fFrom; f != fTo; f += fStp) {
		if (bXMajor && !lyrD.IsXInside(f)) {
			continue;
		} else if (!bXMajor && !lyrD.IsYInside(f)) {
			continue;
		} else {}

		DATA fRatio = (DATA)(f - fFrom) / (fTo - fFrom);
		int s = (int)(sDis * fRatio + sFrom + 0.5F); 
		if (bXMajor && !lyrD.IsYInside(s)) {
			continue;
		} else if (!bXMajor && !lyrD.IsXInside(s)) {
			continue;
		} else {}

		Vect4D<DATA> c(cDis.m_r * fRatio + cFrom.m_r,
					   cDis.m_g * fRatio + cFrom.m_g,
					   cDis.m_b * fRatio + cFrom.m_b,
					   cDis.m_a * fRatio + cFrom.m_a);

		unsigned x = (bXMajor) ? f : s;
		unsigned y = (bXMajor) ? s : f;
		for (unsigned i = 0; i < dimD.m_z; i++) {
			lyrD.CellRef(x, y, i) *= (1.F - c.m_a);
		}
		lyrD.CellRef(x, y, 0) += c.m_r * c.m_a;
		lyrD.CellRef(x, y, 1) += c.m_g * c.m_a;
		lyrD.CellRef(x, y, 2) += c.m_b * c.m_a;
		lyrD.CellRef(x, y, 3) += c.m_a;
	} // f

	if (bDebug) {
		cout << "line ok" << endl;
	} else {}
}

RectangleDraw::RectangleDraw()
{}

RectangleDraw::~RectangleDraw()
{}

void RectangleDraw::Gen(Layer &lyrD, Vect2D<DATA> &p0, Vect2D<DATA> &p1, Vect4D<DATA> &c, bool bFull)
{
	Vect3D<unsigned> dimD = lyrD.GetDim();
	MyAssert(dimD.m_z == 4);

	DATA xL = (p0.m_x < p1.m_x) ? p0.m_x : p1.m_x;
	DATA xR = (p0.m_x < p1.m_x) ? p1.m_x : p0.m_x;
	DATA yB = (p0.m_y < p1.m_y) ? p0.m_y : p1.m_y;
	DATA yT = (p0.m_y < p1.m_y) ? p1.m_y : p0.m_y;

	if (!bFull) {
		Vect2D<DATA> plb(xL, yB);
		Vect2D<DATA> prb(xR, yB);
		Vect2D<DATA> plt(xL, yT);
		Vect2D<DATA> prt(xR, yT);
		myDraw.line.Gen(lyrD, plb, prb, c, c);
		myDraw.line.Gen(lyrD, prb, prt, c, c);
		myDraw.line.Gen(lyrD, prt, plt, c, c);
		myDraw.line.Gen(lyrD, plt, plb, c, c);
	} else {
		int xxL = (int)(xL + 0.5F);
		int xxR = (int)(xR + 0.5F);
		int yyB = (int)(yB + 0.5F);
		int yyT = (int)(yT + 0.5F);
		for (int y = yyB; y <= yyT; y++) {
			if (!lyrD.IsYInside(y)) {
				continue;
			} else {}

			for (int x = xxL; x <= xxR; x++) {
				if (!lyrD.IsXInside(x)) {
					continue;
				} else {}

				for (unsigned i = 0; i < dimD.m_z; i++) {
					lyrD.CellRef(x, y, i) *= (1.F - c.m_a);
				}
				lyrD.CellRef(x, y, 0) += c.m_r * c.m_a;
				lyrD.CellRef(x, y, 1) += c.m_g * c.m_a;
				lyrD.CellRef(x, y, 2) += c.m_b * c.m_a;
				lyrD.CellRef(x, y, 3) += c.m_a;
			}
		}
	}
}

ArrowDraw::ArrowDraw()
{}

ArrowDraw::~ArrowDraw()
{}

void ArrowDraw::Gen(Layer &lyrD, Vect2D<DATA> &pFrom, Vect2D<DATA> &pTo, 
				Vect4D<DATA> &cA, DATA aScl, DATA aAng, bool bDebug)
{
	DATA dX = pTo.m_x - pFrom.m_x;
	DATA dY = pTo.m_y - pFrom.m_y;
	DATA len = sqrt(dX * dX + dY * dY);
	DATA ang = atan2(-dY, -dX);
	myDraw.line.Gen(lyrD, pFrom, pTo, cA, cA);

	DATA ang2 = ang + aAng * D2R;
	DATA len2 = len * aScl;
	Vect2D<DATA> pTo2(pTo.m_x + len2 * cos(ang2),
					  pTo.m_y + len2 * sin(ang2));
	myDraw.line.Gen(lyrD, pTo, pTo2, cA, cA);

	DATA ang3 = ang - aAng * D2R;
	DATA len3 = len * aScl;
	Vect2D<DATA> pTo3(pTo.m_x + len3 * cos(ang3),
					  pTo.m_y + len3 * sin(ang3));
	myDraw.line.Gen(lyrD, pTo, pTo3, cA, cA);
}

PointDraw::PointDraw()
{}

PointDraw::~PointDraw()
{}

void PointDraw::Gen(Layer &lyrD, Vect2D<DATA> &pLoc, DATA r, Vect4D<DATA> &cA, bool bDebug)
{
	Vect3D<unsigned> dimD = lyrD.GetDim();
	MyAssert(dimD.m_z == 4 || dimD.m_z == 3);

	int xL = (int)(pLoc.m_x - r + 0.5F);
	int xR = (int)(pLoc.m_x + r + 0.5F);
	int yB = (int)(pLoc.m_y - r + 0.5F);
	int yT = (int)(pLoc.m_y + r + 0.5F);

	for (int y = yB; y <= yT; y++)  {
		if (!lyrD.IsYInside(y)) {
			continue;
		} else {}
		DATA yDis = y - pLoc.m_y;
		DATA yyD = yDis * yDis;

		for (int x = xL; x <= xR; x++) {
			if (!lyrD.IsXInside(x)) {
				continue;
			} else {}
			DATA xDis = x - pLoc.m_x;
			DATA xxD = xDis * xDis;

			if (xxD + yyD > r * r) {
				continue;
			} else {}

			for (unsigned i = 0; i < dimD.m_z; i++) {
				lyrD.CellRef(x, y, i) *= (1.F - cA.m_a);
			}
			lyrD.CellRef(x, y, 0) += cA.m_r * cA.m_a;
			lyrD.CellRef(x, y, 1) += cA.m_g * cA.m_a;
			lyrD.CellRef(x, y, 2) += cA.m_b * cA.m_a;
			if (dimD.m_z == 4) {
				lyrD.CellRef(x, y, 3) += cA.m_a;
			} else {}
		}
	}
}

BlendDraw::BlendDraw()
{}

BlendDraw::~BlendDraw()
{}

void BlendDraw::Gen(Layer &lyrDraw, const Mtx &mtxIn)
{
	Vect3D<unsigned> dimD = lyrDraw.GetDim();
	Vect2D<unsigned> dimIn = mtxIn.GetDim();
	MyAssert(dimD.m_x == dimIn.m_x &&
			 dimD.m_y == dimIn.m_y);
	MyAssert(dimD.m_z == 4);

	for (unsigned y = 0; y < dimD.m_y; y++) {
		for (unsigned x = 0; x < dimD.m_x; x++) {
			for (unsigned c = 0; c < 3; c++) {
				lyrDraw.CellRef(x, y, c) += mtxIn.CellVal(x, y) * (1.F - lyrDraw.CellVal(x, y, 3));
			}
		}
	}
}