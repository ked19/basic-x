#ifndef _IMAGE_IO_H
#define _IMAGE_IO_H

#include "layer.h"
//#include "matrixOperation.h"

#include "FreeImage.h"
#include <iostream>
#include <string>
#include <cassert>

using namespace std;

class ImgIO
{
public:
	ImgIO();
	~ImgIO();

	MyImg* Read(const string &fName) const;
	void Read(MyImg &img, const string &fName) const;

	void Write(const string &fName, MyImg &img) const;	

	bool IsBmpFormat(MyImg &img) const;

	static void Test();

private:
	FIBITMAP* Load2Bmp(const string &fName) const;	
	FREE_IMAGE_FORMAT DecideExtFormat(const string &fName) const;
};

extern ImgIO imgIO;

#endif