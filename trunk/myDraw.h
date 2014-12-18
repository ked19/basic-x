#ifndef _MY_DRAW_H
#define _MY_DRAW_H

#include "layer.h"
#include "define.h"

#include <iostream>
#include <cmath>

using namespace std;

class LineDraw
{
public:
	LineDraw();
	~LineDraw();

	void Gen(Layer &lyrD, Vect2D<DATA> &pFrom, Vect2D<DATA> &pTo, 
						  Vect4D<DATA> &cFrom, Vect4D<DATA> &cTo, bool bDebug = false);
};

class RectangleDraw
{
public:
	RectangleDraw();
	~RectangleDraw();

	void Gen(Layer &lyrD, Vect2D<DATA> &p0, Vect2D<DATA> &p1, Vect4D<DATA> &c, bool bFull = false);
};

class ArrowDraw
{
public:
	ArrowDraw();
	~ArrowDraw();

	void Gen(Layer &lyrD, Vect2D<DATA> &pFrom, Vect2D<DATA> &pTo, 
			 Vect4D<DATA> &cA, DATA aScl = 0.7F, DATA aAng = 30.F, bool bDebug = false);
};

class PointDraw
{
public:
	PointDraw();
	~PointDraw();

	void Gen(Layer &lyrD, Vect2D<DATA> &pLoc, DATA r, Vect4D<DATA> &cA, bool bDebug = false);

};

class BlendDraw
{
public:
	BlendDraw();
	~BlendDraw();

	void Gen(Layer &lyrDraw, const Mtx &mtxIn);
};

class MyDraw
{
public:
	LineDraw				line;
	RectangleDraw 			rect;
	ArrowDraw 				arrow;
	PointDraw 				point;
	BlendDraw 				blend;

private:
};

extern MyDraw myDraw;

#endif