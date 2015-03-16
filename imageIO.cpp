#include "imageIO.h"

ImgIO imgIO;

ImgIO::ImgIO()
{
}

ImgIO::~ImgIO()
{
}

//*************************************************************************************************

FIBITMAP* ImgIO::Load2Bmp(const string &fName) const
{
	FREE_IMAGE_FORMAT aFormat[] = {FIF_BMP, FIF_JPEG, FIF_PPM, FIF_PNG, FIF_TIFF, FIF_PGMRAW, FIF_PNG};

	unsigned fNum = sizeof(aFormat) / sizeof(FREE_IMAGE_FORMAT);
	for(unsigned i = 0; i < fNum; i++) {
		FIBITMAP *bmp = FreeImage_Load(aFormat[i], fName.c_str());
		if(bmp) {
			return bmp;
		} else {}
	}
	MyAssert(0);
	return 0;
}

void GetPalette1V(MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	BYTE pIdx;
	FreeImage_GetPixelIndex(&bmp, x, y, &pIdx);
	img.CellRef(x, y, 0) = (float)pIdx;
}

void GetPalette3V(MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	BYTE pIdx;
	FreeImage_GetPixelIndex(&bmp, x, y, &pIdx);
	RGBQUAD color = palette[pIdx];
	img.CellRef(x, y, 0) = (float)color.rgbRed;
	img.CellRef(x, y, 1) = (float)color.rgbGreen;
	img.CellRef(x, y, 2) = (float)color.rgbBlue;
}

void GetNonpalette3V(MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	RGBQUAD color;
	FreeImage_GetPixelColor(&bmp, x, y, &color);
	img.CellRef(x, y, 0) = (float)color.rgbRed;
	img.CellRef(x, y, 1) = (float)color.rgbGreen;
	img.CellRef(x, y, 2) = (float)color.rgbBlue;
}

void GetNonpalette4V(MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	RGBQUAD color;
	FreeImage_GetPixelColor(&bmp, x, y, &color);
	img.CellRef(x, y, 0) = (float)color.rgbRed;
	img.CellRef(x, y, 1) = (float)color.rgbGreen;
	img.CellRef(x, y, 2) = (float)color.rgbBlue;
	img.CellRef(x, y, 3) = (float)color.rgbReserved;
}

void GetUint16(MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	unsigned short *pBits = (unsigned short*)FreeImage_GetScanLine(&bmp, y);
	img.CellRef(x, y, 0) = (float)pBits[x];
}

void ImgIO::Read(MyImg **ppImg, MyImg &imgIn, const string &fName) const
{
	FIBITMAP *pBmp = Load2Bmp(fName);

	unsigned xDim = FreeImage_GetWidth(pBmp);
	unsigned yDim = FreeImage_GetHeight(pBmp);

	FREE_IMAGE_TYPE dataType = FreeImage_GetImageType(pBmp);
	assert(dataType==FIT_BITMAP // 1, 4, 8, 16, 32 bits
		|| dataType==FIT_UINT16);

	FREE_IMAGE_COLOR_TYPE colorType = FreeImage_GetColorType(pBmp);
	unsigned cDim;
	FreeImage_Unload(pBmp);	

	if(dataType == FIT_BITMAP) {
		if(colorType == FIC_MINISBLACK)			{cDim = 1;}
		else if(colorType == FIC_MINISWHITE)	{cDim = 1;}
		else if(colorType == FIC_PALETTE)		{cDim = 3;}
		else if(colorType == FIC_RGB)			{cDim = 3;}
		else if(colorType == FIC_RGBALPHA)		{cDim = 4;}
		else if(colorType == FIC_CMYK)			assert(0);
		else									assert(0);
	} else if(dataType == FIT_UINT16) {
		if(colorType == FIC_MINISBLACK)			{cDim = 1;}
		else									assert(0);
	} else {
		assert(0);
	}
	
	Vect3D<unsigned> dimImg = imgIn.GetDim();
	MyAssert(dimImg.m_x >= xDim &&
			 dimImg.m_y >= yDim &&
			 dimImg.m_z == cDim);
	*ppImg = new MyImg(imgIn, 0, 0, xDim, yDim);
	Read(**ppImg, fName);
}

void ImgIO::Read(MyImg &img, const string &fName) const
{
	cout << "read image file: " << fName << endl;

	FIBITMAP *pBmp = Load2Bmp(fName);

	unsigned xDim = FreeImage_GetWidth(pBmp);
	unsigned yDim = FreeImage_GetHeight(pBmp);
	unsigned zDim = 1;
	cout << "xDim:" << xDim << " " <<
		    "yDim:" << yDim << " " <<
			"zDim:" << zDim << endl;

	FREE_IMAGE_TYPE dataType = FreeImage_GetImageType(pBmp);
	assert(dataType==FIT_BITMAP // 1, 4, 8, 16, 32 bits
		|| dataType==FIT_UINT16);

	FREE_IMAGE_COLOR_TYPE colorType = FreeImage_GetColorType(pBmp);
	unsigned cDim;
	void (*GetValue)(MyImg&, FIBITMAP&, RGBQUAD*, unsigned, unsigned) = 0;
	if(dataType == FIT_BITMAP)
	{
		if(colorType == FIC_MINISBLACK)			{cDim = 1;		GetValue = GetPalette1V;}
		else if(colorType == FIC_MINISWHITE)	{cDim = 1;		GetValue = GetPalette1V;}
		else if(colorType == FIC_PALETTE)		{cDim = 3;		GetValue = GetPalette3V;}
		else if(colorType == FIC_RGB)			{cDim = 3;		GetValue = GetNonpalette3V;}
		else if(colorType == FIC_RGBALPHA)		{cDim = 4;		GetValue = GetNonpalette4V;}
		else if(colorType == FIC_CMYK)			assert(0);
		else									assert(0);
	}
	else if(dataType == FIT_UINT16)
	{
		if(colorType == FIC_MINISBLACK)			{cDim = 1;		GetValue = GetUint16;}
		else									assert(0);
	}
	else
	{
		assert(0);
	}
	cout << "cDim:" << cDim << endl;

	Vect3D<unsigned> dimImg = img.GetDim();
	assert(
		dimImg.m_x==xDim &&
		dimImg.m_y==yDim &&
		dimImg.m_z==cDim);

	RGBQUAD *palette = FreeImage_GetPalette(pBmp);
	if(palette)	cout<<"palette:yes" << endl;
	else		cout<<"palette:no"  << endl;

	cout << "read value.." << endl;
	for(unsigned y=0; y<yDim; y++)
	{
		for(unsigned x=0; x<xDim; x++)
		{
			GetValue(img, *pBmp, palette, x, y); 
		}
	}

	FreeImage_Unload(pBmp);	
	cout << "read ok" << endl << endl;
}

MyImg* ImgIO::Read(const std::string &fName) const
{
	cout << "read image file: " << fName << endl;

	FIBITMAP *pBmp = Load2Bmp(fName);

	unsigned xDim = FreeImage_GetWidth(pBmp);
	unsigned yDim = FreeImage_GetHeight(pBmp);
	unsigned zDim = 1;

	FREE_IMAGE_TYPE dataType = FreeImage_GetImageType(pBmp);
	assert(dataType==FIT_BITMAP // 1, 4, 8, 16, 32 bits
		|| dataType==FIT_UINT16);

	FREE_IMAGE_COLOR_TYPE colorType = FreeImage_GetColorType(pBmp);
	unsigned cDim;
	if(dataType == FIT_BITMAP)
	{
		if(colorType == FIC_MINISBLACK)			{cDim = 1;}
		else if(colorType == FIC_MINISWHITE)	{cDim = 1;}
		else if(colorType == FIC_PALETTE)		{cDim = 3;}
		else if(colorType == FIC_RGB)			{cDim = 3;}
		else if(colorType == FIC_RGBALPHA)		{cDim = 4;}
		else if(colorType == FIC_CMYK)			assert(0);
		else									assert(0);
	}
	else if(dataType == FIT_UINT16)
	{
		if(colorType == FIC_MINISBLACK)			{cDim = 1;}
		else									assert(0);
	}
	else
	{
		assert(0);
	}

	MyImg *pImg = new MyImg(xDim, yDim, cDim);
	Read(*pImg, fName);
	return pImg;
}

//*************************************************************************************************

FREE_IMAGE_FORMAT ImgIO::DecideExtFormat(const string &fName) const
{
	size_t loc = fName.rfind('.');
	MyAssert(loc != string::npos);
	string ext = fName.substr(loc+1);

	if (!ext.compare("bmp")) {
		return FIF_BMP;
	} else if (!ext.compare("tif") || !ext.compare("tiff")) {
		return FIF_TIFF;
	} else if (!ext.compare("jpg") || !ext.compare("jpeg")) {
		return FIF_JPEG;
	} else if (!ext.compare("pgm")) {
		return FIF_PGMRAW;
	} else if (!ext.compare("png")) {
		return FIF_PNG;
	} else {
		MyAssert(0);
	}
	return FIF_UNKNOWN;
}

void SetNonpalette1V(const Mtx &mtx, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	RGBQUAD color;
	color.rgbRed		= (BYTE)mtx.CellVal(x, y);
	color.rgbGreen		= color.rgbRed;
	color.rgbBlue		= color.rgbRed;
	color.rgbReserved	= 0;
	FreeImage_SetPixelColor(&bmp, x, y, &color);	
}

void SetNonpalette1V(const MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	RGBQUAD color;
	color.rgbRed		= (BYTE)img.CellVal(x, y, 0);
	color.rgbGreen		= color.rgbRed;
	color.rgbBlue		= color.rgbRed;
	color.rgbReserved	= 0;
	FreeImage_SetPixelColor(&bmp, x, y, &color);	
}

void SetNonpalette2V(const MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	RGBQUAD color;
	color.rgbRed		= (BYTE)img.CellVal(x, y, 0);
	color.rgbGreen		= (BYTE)img.CellVal(x, y, 1);
	color.rgbBlue		= 0;
	color.rgbReserved	= 0;
	FreeImage_SetPixelColor(&bmp, x, y, &color);	
}

void SetNonpalette3V(const MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	RGBQUAD color;
	color.rgbRed		= (BYTE)img.CellVal(x, y, 0);
	color.rgbGreen		= (BYTE)img.CellVal(x, y, 1);
	color.rgbBlue		= (BYTE)img.CellVal(x, y, 2);
	color.rgbReserved	= 0;
	FreeImage_SetPixelColor(&bmp, x, y, &color);	
}

void SetNonpalette4V(const MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	RGBQUAD color;
	color.rgbRed		= (BYTE)img.CellVal(x, y, 0);
	color.rgbGreen		= (BYTE)img.CellVal(x, y, 1);
	color.rgbBlue		= (BYTE)img.CellVal(x, y, 2);
	color.rgbReserved	= (BYTE)img.CellVal(x, y, 3);
	FreeImage_SetPixelColor(&bmp, x, y, &color);	
}

void SetUint8(const MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	unsigned char *pBits = (unsigned char*)FreeImage_GetScanLine(&bmp, y);
	pBits[x] = (unsigned char)img.CellVal(x, y, 0);
}

void SetUint16(const MyImg &img, FIBITMAP &bmp, RGBQUAD *palette, unsigned x, unsigned y)
{
	unsigned short *pBits = (unsigned short*)FreeImage_GetScanLine(&bmp, y);
	pBits[x] = (unsigned short)img.CellVal(x, y, 0);
}

void ImgIO::Write(const std::string &fName, MyImg &img) const
{
	cout << "write image file: " << fName.c_str() << endl;

	Vect3D<unsigned> dim = img.GetDim();
	cout << "Dim(x,y,c)= " << dim.m_x << ", " 
						   << dim.m_y << ", " 
						   << dim.m_z << endl;

	FREE_IMAGE_FORMAT dataFormat = DecideExtFormat(fName);
	FREE_IMAGE_TYPE dataType;
	unsigned bpp;
	int flag;
	if (dataFormat == FIF_BMP) {
		dataType = FIT_BITMAP;
		bpp = 32;
		flag = BMP_DEFAULT;
	} else if (dataFormat == FIF_TIFF) {
		dataType = FIT_UINT16;
		bpp = 16;
		flag = TIFF_NONE;
	} else if (dataFormat == FIF_JPEG) {
		dataType = FIT_BITMAP;
		bpp = 24;
		flag = JPEG_QUALITYSUPERB;
	} else if (dataFormat == FIF_PGMRAW) {
		dataType = FIT_BITMAP;
		bpp = 8;
		flag = PNM_DEFAULT;
	} else if (dataFormat == FIF_PNG) {
		dataType = FIT_BITMAP;
		bpp = 24;
		flag = PNG_Z_NO_COMPRESSION;
	} else {
		MyAssert(0);
	}

	void (*SetValue)(const MyImg&, FIBITMAP&, RGBQUAD*, unsigned, unsigned) = 0;
	if (dataFormat == FIF_BMP) {
		if (dim.m_z == 1)		SetValue = SetNonpalette1V;
		else if (dim.m_z == 2)	SetValue = SetNonpalette2V;
		else if (dim.m_z == 3)	SetValue = SetNonpalette3V;
		else if (dim.m_z == 4)	SetValue = SetNonpalette4V;
		else					MyAssert(0);
	} else if (dataFormat == FIF_TIFF) {
		SetValue = SetUint16;
	} else if (dataFormat == FIF_JPEG) {
		if (dim.m_z == 1)		SetValue = SetNonpalette1V;
		else if (dim.m_z == 3)	SetValue = SetNonpalette3V;
		else					MyAssert(0);
	} else if (dataFormat == FIF_PGMRAW) {
		if (dim.m_z == 1)		SetValue = SetUint8;
		else 					MyAssert(0);
	} else if (dataFormat == FIF_PNG) {
		if (dim.m_z == 1)		SetValue = SetNonpalette1V;
		else if (dim.m_z == 3)	SetValue = SetNonpalette3V;
		else					MyAssert(0);
	} else {
		MyAssert(0);
	}

	cout << "write value.." << endl;
	FIBITMAP *pBmp = FreeImage_AllocateT(dataType, dim.m_x, dim.m_y, (int)bpp);
	RGBQUAD *palette = 0;
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			SetValue(img, *pBmp, palette, x, y);
		}
	}
	FreeImage_Save(dataFormat, pBmp, fName.c_str(), flag);
	FreeImage_Unload(pBmp);

	cout << "write ok" << endl << endl;
}

//*************************************************************************************************

bool ImgIO::IsBmpFormat(MyImg &img) const
{
	cout << "test if matrix is bitmap format.." << endl;

	bool bBmp = true;
	Vect3D<unsigned> dim = img.GetDim();
	if(dim.m_z > 4)
	{
		bBmp = false;
		cout << "number of channel: " << dim.m_z << endl;
	}
	else {}

	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			bool bPass = true;
			for(unsigned c=0; c<dim.m_z; c++)
			{
				DATA val = img.CellVal(x, y, c);
				if(val == DATA_NAN)
				{
					bPass = false;
					break;	
				}
				else if(val<0 || val>255.F)
				{
					bPass = false;
					break;
				}
				else {}

				if(!bPass)
				{
					assert(0);
					/*
					cout << "img["<<x<<"]" << "["<<y<<"]" << "["<<c<<"] = ";
					Output nOut;
					nOut << val;
					cout << endl;
					*/
				}
				else {}
			} // c
		} // x
	} // y
	cout << "test end"  << endl;
	return bBmp;	
}

void ImgIO::Test()
{
	ImgIO imgIO;
	cout << "start ImageIO testing.." << endl;

	string inBmpName  = TEST_DATA_DIR	+ "\\bmpImg.bmp"; 
	string outBmpName = OUTPUT_DIR		+ "\\bmpOut.bmp";
	MyImg *pImgBmp = imgIO.Read(inBmpName);
	assert(imgIO.IsBmpFormat(*pImgBmp));
	imgIO.Write(outBmpName, *pImgBmp);
	delete pImgBmp;
	cout << "IO of bmp ok" << endl << endl;

	//************************************************************************************************
	
	string inJpgName  = TEST_DATA_DIR	+ "\\jpgImg.jpg";
	string outJpgName = OUTPUT_DIR		+ "\\jpgOut.jpg";
	MyImg *pImgJpg = imgIO.Read(inJpgName);
	assert(imgIO.IsBmpFormat(*pImgJpg));
	imgIO.Write(outJpgName, *pImgJpg);
	delete pImgJpg;
	cout << "IO of jpg ok" << endl << endl;

	cout << "ImageIO testing ok" << endl << endl;
}