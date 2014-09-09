#ifndef _MY_DRAW_H
#define _MY_DRAW_H

#include "layer.h"
#include "define.h"

#include <iostream>
#include <cmath>

using namespace std;

class Line
{
public:
	Line();
	~Line();

	void Gen(Layer &lyrD, Vect2D<DATA> &pFrom, Vect2D<DATA> &pTo, 
						  Vect4D<DATA> &cFrom, Vect4D<DATA> &cTo, bool bDebug = false);
};

class Rectangle
{
public:
	Rectangle();
	~Rectangle();

	void Gen(Layer &lyrD, Vect2D<DATA> &p0, Vect2D<DATA> &p1, Vect4D<DATA> &c, bool bFull = false);
};

class Arrow
{
public:
	Arrow();
	~Arrow();

	void Gen(Layer &lyrD, Vect2D<DATA> &pFrom, Vect2D<DATA> &pTo, 
			 Vect4D<DATA> &cA, DATA aScl = 0.7F, DATA aAng = 30.F, bool bDebug = false);
};

class Point
{
public:
	Point();
	~Point();

	void Gen(Layer &lyrD, Vect2D<DATA> &pLoc, DATA r, Vect4D<DATA> &cA, bool bDebug = false);

};

class Blend
{
public:
	Blend();
	~Blend();

	void Gen(Layer &lyrDraw, const Mtx &mtxIn);
};

class MyDraw
{
public:
	Line				line;
	Rectangle 			rect;
	Arrow 				arrow;
	Point 				point;
	Blend 				blend;

private:
};

extern MyDraw myDraw;

#endif