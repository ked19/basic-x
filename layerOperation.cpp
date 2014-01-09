#include "layerOperation.h"

G_LayerOp lyrOp;

//*************************************************************************************************

SaveLyr::SaveLyr()
{}

SaveLyr::~SaveLyr()
{}

void SaveLyr::Gen(const Layer &lyr, string f, bool bB) const
{
	if(bB)
	{
		ofstream outF(f.c_str(), ios::binary);

		Vect3D<unsigned> dim = lyr.GetDim();
		outF.write((char*)&dim.m_x, sizeof(unsigned));
		outF.write((char*)&dim.m_y, sizeof(unsigned));
		outF.write((char*)&dim.m_z, sizeof(unsigned));

		for(unsigned z=0; z<dim.m_z; z++)
		{
			for(unsigned y=0; y<dim.m_y; y++)
			{
				for(unsigned x=0; x<dim.m_x; x++)
				{
					DATA v = lyr.CellVal(x, y, z);
					/*
					cout << v << " ";
					if(x==dim.m_x-1 && y==dim.m_y-1 && z==dim.m_z-1)
					{
						cout << endl;
					}
					else {}
					*/

					outF.write((char*)&v, sizeof(DATA));
				}
			}
		}

		outF.close();
	}
	else
	{
		MyAssert(0);
	}
}

SaveImgLyr::SaveImgLyr()
{}

SaveImgLyr::~SaveImgLyr()
{}

void SaveImgLyr::Gen(const Layer & lyr, string f, string ext) const
{
	stringstream ss;
	string sN;
	Vect3D<unsigned> dim = lyr.GetDim();
    for (unsigned i = 0; i < dim.m_z; i++) {
        ss.clear();
        ss << setfill('0') << setw(4) << i;
        ss >> sN;
        string sF = f + "_" + sN + ".bmp";
        imgIO.Write(sF, MyImg(*lyr.GetMtx(i)));
    } // i
}

//*************************************************************************************************

ZeroLyr::ZeroLyr()
{}

ZeroLyr::~ZeroLyr()
{}

void ZeroLyr::Gen(Layer &lyr) const
{
	Vect3D<unsigned> dim = lyr.GetDim();
	for(unsigned c=0; c<dim.m_z; c++)
	{
		Mtx *pMtx = lyr.GetMtx(c);
		m_zero.Gen(*pMtx);
	}
}

OneLyr::OneLyr()
{}

OneLyr::~OneLyr()
{}

void OneLyr::Gen(Layer &lyr) const
{
	Vect3D<unsigned> dim = lyr.GetDim();
	for(unsigned c=0; c<dim.m_z; c++)
	{
		Mtx *pMtx = lyr.GetMtx(c);
		mtxOp.one.Gen(*pMtx);
	}
}

RndLyr::RndLyr()
{}

RndLyr::~RndLyr()
{}

void RndLyr::Gen(Layer &lyr) const
{
	srand((unsigned)time(0));

	Vect3D<unsigned> dim = lyr.GetDim();
	for(unsigned z=0; z<dim.m_z; z++)
	{
		for(unsigned y=0; y<dim.m_y; y++)
		{
			for(unsigned x=0; x<dim.m_x; x++)
			{
				lyr.CellRef(x, y, z) = (DATA)rand() / RAND_MAX;
			}
		}
	}
}

Gauss3DLyr::Gauss3DLyr()
{}

Gauss3DLyr::~Gauss3DLyr()
{}

void Gauss3DLyr::Gen(Layer &lyr, DATA s, bool bNormal) const
{
	Vect3D<unsigned> vDim = lyr.GetDim();
	Vect3D<DATA> vCnt((vDim.m_x-1.F)/2, (vDim.m_y-1.F)/2, (vDim.m_z-1.F)/2);

	DATA sum = 0;
	for(unsigned z=0; z<vDim.m_z; z++)
	{
		DATA zDis = (z - vCnt.m_z);
		for(unsigned y=0; y<vDim.m_y; y++)
		{
			DATA yDis = (y - vCnt.m_y);
			for(unsigned x=0; x<vDim.m_x; x++)
			{
				DATA xDis = (x - vCnt.m_x);

				DATA ex = -(xDis*xDis + yDis*yDis + zDis*zDis) / (2*s*s);
				lyr.CellRef(x, y, z) = exp(ex);
				//cout << lyr.CellVal(x, y, z) << " ";
				sum += lyr.CellVal(x, y, z);
			}
		}
	}

	if(bNormal)
	{
		for(unsigned z=0; z<vDim.m_z; z++)
		{
			for(unsigned y=0; y<vDim.m_y; y++)
			{
				for(unsigned x=0; x<vDim.m_x; x++)
				{
					lyr.CellRef(x, y, z) = lyr.CellVal(x, y, z) / sum;
				}
			}
		}
	}
	else {}
}

void Gauss3DLyr::Gen(VolumeData &vol, DATA s, bool bNormal) const
{
	Vect3D<unsigned> vDim = vol.GetDim();
	Vect3D<DATA> vStp = vol.GetStep();
	Vect3D<DATA> vCnt((vDim.m_x-1.F)*vStp.m_x/2, 
					  (vDim.m_y-1.F)*vStp.m_y/2, 
					  (vDim.m_z-1.F)*vStp.m_z/2);

	DATA sum = 0;
	for(unsigned z=0; z<vDim.m_z; z++)
	{
		DATA zDis = (z*vStp.m_z - vCnt.m_z);
		for(unsigned y=0; y<vDim.m_y; y++)
		{
			DATA yDis = (y*vStp.m_y - vCnt.m_y);
			for(unsigned x=0; x<vDim.m_x; x++)
			{
				DATA xDis = (x*vStp.m_x - vCnt.m_x);

				DATA ex = -(xDis*xDis + yDis*yDis + zDis*zDis) / (2*s*s);
				vol.CellRef(x, y, z) = exp(ex);
				//cout << lyr.CellVal(x, y, z) << " ";
				sum += vol.CellVal(x, y, z);
			}
		}
	}

	if(bNormal)
	{
		for(unsigned z=0; z<vDim.m_z; z++)
		{
			for(unsigned y=0; y<vDim.m_y; y++)
			{
				for(unsigned x=0; x<vDim.m_x; x++)
				{
					vol.CellRef(x, y, z) = vol.CellVal(x, y, z) / sum;
				}
			}
		}
	}
	else {}
}

Norm3DLyr::Norm3DLyr()
{}

Norm3DLyr::~Norm3DLyr()
{}

DATA Norm3DLyr::Gen(Layer &lyr, bool bNor) const
{
	Vect3D<unsigned> dim = lyr.GetDim();

	DATA len = 0;
	for(unsigned z=0; z<dim.m_z; z++)
	{
		for(unsigned y=0; y<dim.m_y; y++)
		{
			for(unsigned x=0; x<dim.m_x; x++)
			{
				len += lyr.CellVal(x, y, z) * lyr.CellVal(x, y, z);
			}
		}
	}
	len = sqrt(len);

	if(bNor)
	{
		if(myMath.IsEqual(len, 0))
		{
			return len;
		}
		else
		{
			lyrOp.mul.Gen(lyr, 1.F/len);
			return len;
		}
	}
	else 
	{
		return len;
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

HessianLyr::HessianLyr()
{}

HessianLyr::~HessianLyr()
{}

void HessianLyr::Gen(Mtx2D &mtxH, const Layer &lyr, unsigned x, unsigned y, unsigned c) const
{	
	Vect3D<unsigned> dim = lyr.GetDim();

	DATA xxDer, yyDer, zzDer, xyDer, xzDer, yzDer;
	if(x>=1 && x<=dim.m_x-2 &&
	   y>=1 && y<=dim.m_y-2 &&
	   c>=1 && c<=dim.m_z-2)
	{
		unsigned xR = x + 1;	unsigned xL = x - 1;
		unsigned yT = y + 1;	unsigned yB = y - 1;
		unsigned cF = c + 1;	unsigned cN = c - 1;

		xxDer = lyr.CellVal(xL,  y,  c) + lyr.CellVal(xR,  y,  c) - 2 * lyr.CellVal(x, y, c);
		yyDer = lyr.CellVal( x, yB,  c) + lyr.CellVal( x, yT,  c) - 2 * lyr.CellVal(x, y, c);
		zzDer = lyr.CellVal( x,  y, cN) + lyr.CellVal( x,  y, cF) - 2 * lyr.CellVal(x, y, c);

		xyDer = (lyr.CellVal(xR, yT,  c) - lyr.CellVal(xR, yB,  c) - lyr.CellVal(xL, yT,  c) + lyr.CellVal(xL, yB,  c)) * 0.25F;
		xzDer = (lyr.CellVal(xR,  y, cF) - lyr.CellVal(xR,  y, cN) - lyr.CellVal(xL,  y, cF) + lyr.CellVal(xL,  y, cN)) * 0.25F;
		yzDer = (lyr.CellVal( x, yT, cF) - lyr.CellVal( x, yT, cN) - lyr.CellVal( x, yB, cF) + lyr.CellVal(x,  yB, cN)) * 0.25F;
	}
	else
	{
		Layer subLyr(3, 3, 3);
		lyrOp.zero.Gen(subLyr);

		for(unsigned cc=0; cc<3; cc++)
		{
			int cLoc = cc - 1 + c;
			if(cLoc<0 || cLoc>=(int)dim.m_z)	continue;
			else {}

			for(unsigned yy=0; yy<3; yy++)
			{
				int yLoc = yy - 1 + y;
				if(yLoc<0 || yLoc>=(int)dim.m_y)	continue;
				else {}

				for(unsigned xx=0; xx<3; xx++)
				{
					int xLoc = xx - 1 + x;
					if(xLoc<0 || xLoc>=(int)dim.m_x)	continue;
					else {}

					subLyr.CellRef(xx, yy, cc) = lyr.CellVal(xLoc, yLoc, cLoc);
				}
			}
		}

		xxDer = subLyr.CellVal(0, 1, 1) + subLyr.CellVal(2, 1, 1) - 2 * subLyr.CellVal(1, 1, 1);
		yyDer = subLyr.CellVal(1, 0, 1) + subLyr.CellVal(1, 2, 1) - 2 * subLyr.CellVal(1, 1, 1);
		zzDer = subLyr.CellVal(1, 1, 0) + subLyr.CellVal(1, 1, 2) - 2 * subLyr.CellVal(1, 1, 1);

		xyDer = (subLyr.CellVal(2, 2, 1) - subLyr.CellVal(2, 0, 1) - subLyr.CellVal(0, 2, 1) + subLyr.CellVal(0, 0, 1)) * 0.25F;
		xzDer = (subLyr.CellVal(2, 1, 2) - subLyr.CellVal(2, 1, 0) - subLyr.CellVal(0, 1, 2) + subLyr.CellVal(0, 1, 0)) * 0.25F;
		yzDer = (subLyr.CellVal(1, 2, 2) - subLyr.CellVal(1, 2, 0) - subLyr.CellVal(1, 0, 2) + subLyr.CellVal(1, 0, 0)) * 0.25F;
	}

	mtxH.CellRef(0, 0) = xxDer;		mtxH.CellRef(1, 0) = xyDer;		mtxH.CellRef(2, 0) = xzDer;
	mtxH.CellRef(0, 1) = xyDer;		mtxH.CellRef(1, 1) = yyDer;		mtxH.CellRef(2, 1) = yzDer;
	mtxH.CellRef(0, 2) = xzDer;		mtxH.CellRef(1, 2) = yzDer;		mtxH.CellRef(2, 2) = zzDer;
}

//*************************************************************************************************
//
//*************************************************************************************************

DerivateLyr::DerivateLyr()
{}

DerivateLyr::~DerivateLyr()
{}

void DerivateLyr::Gen(Mtx &mtxDer, const Layer &lyr, unsigned x, unsigned y, unsigned c) const
{	
	Vect3D<unsigned> dim = lyr.GetDim();

	if(x>=1 && x<=dim.m_x-2 &&
	   y>=1 && y<=dim.m_y-2 &&
	   c>=1 && c<=dim.m_z-2)
	{
		unsigned xR = x + 1;	unsigned xL = x - 1;
		unsigned yT = y + 1;	unsigned yB = y - 1;
		unsigned cF = c + 1;	unsigned cN = c - 1;

		mtxDer.CellRef(0, 0) = (lyr.CellVal(xR,  y,  c) - lyr.CellVal(xL,  y,  c)) * 0.5F;
		mtxDer.CellRef(0, 1) = (lyr.CellVal( x, yT,  c) - lyr.CellVal( x, yB,  c)) * 0.5F;
		mtxDer.CellRef(0, 2) = (lyr.CellVal( x,  y, cF) - lyr.CellVal( x,  y, cN)) * 0.5F;
	}
	else
	{
		Layer subLyr(3, 3, 3);
		for(unsigned cc=0; cc<3; cc++)
		{
			for(unsigned yy=0; yy<3; yy++)
			{
				for(unsigned xx=0; xx<3; xx++)
				{
					subLyr.CellRef(xx, yy, cc) = lyr.CellVal(x, y, c);
				}
			}
		}
		for(unsigned cc=0; cc<3; cc++)
		{
			int cLoc = cc - 1 + c;
			if(cLoc<0 || cLoc>=(int)dim.m_z)	continue;
			else {}

			for(unsigned yy=0; yy<3; yy++)
			{
				int yLoc = yy - 1 + y;
				if(yLoc<0 || yLoc>=(int)dim.m_y)	continue;
				else {}

				for(unsigned xx=0; xx<3; xx++)
				{
					int xLoc = xx - 1 + x;
					if(xLoc<0 || xLoc>=(int)dim.m_x)	continue;
					else {}

					subLyr.CellRef(xx, yy, cc) = lyr.CellVal(xLoc, yLoc, cLoc);
				}
			}
		}

		mtxDer.CellRef(0, 0) = (subLyr.CellVal(2, 1, 1) - subLyr.CellVal(0, 1, 1)) * 0.5F;
		mtxDer.CellRef(0, 1) = (subLyr.CellVal(1, 2, 1) - subLyr.CellVal(1, 0, 1)) * 0.5F;
		mtxDer.CellRef(0, 2) = (subLyr.CellVal(1, 1, 2) - subLyr.CellVal(1, 1, 0)) * 0.5F;
	}
}

void DerivateLyr::Gen(Layer &lyrDer, const Mtx &mtx, bool bNormal) const
{
	Vect2D<unsigned> dim = mtx.GetDim();
	
	Mtx mtxDer(1, 2);
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			mtxOp.der.Gen(mtxDer, mtx, x, y, bNormal);
			lyrDer.CellRef(x, y, 0) = mtxDer.CellVal(0, 0);
			lyrDer.CellRef(x, y, 1) = mtxDer.CellVal(0, 1);
		}
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

Orient::Orient()
{}

Orient::~Orient()
{}

void Orient::Gen(Layer &lyrOrnt, const Mtx &mtx)
{
	Vect2D<unsigned> dim = mtx.GetDim();

	Mtx mtxDer(1, 2);
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			mtxOp.der.Gen(mtxDer, mtx, x, y, false);
			DATA dx = mtxDer.CellVal(0, 0);
			DATA dy = mtxDer.CellVal(0, 1);

			lyrOrnt.CellRef(x, y, 0) = sqrt(dx*dx + dy*dy);
			lyrOrnt.CellRef(x, y, 1) = atan2(dy, dx);

		}
	}
}



//*************************************************************************************************
//
//*************************************************************************************************

ScaleDimLyr::ScaleDimLyr()
{}

ScaleDimLyr::~ScaleDimLyr()
{}

void ScaleDimLyr::Gen_nearest(Layer &lyrOut, const Layer &lyrIn) const
{
	Vect3D<unsigned> dimIn  = lyrIn.GetDim();
	Vect3D<unsigned> dimOut = lyrOut.GetDim();

	DATA sX = (DATA)(dimIn.m_x-1) / (dimOut.m_x-1);
	DATA sY = (DATA)(dimIn.m_y-1) / (dimOut.m_y-1);
	DATA sZ = (DATA)(dimIn.m_z-1) / (dimOut.m_z-1);

	for(unsigned z=0; z<dimOut.m_z; z++)
	{
		unsigned zInv = (unsigned)(z*sZ + 0.5F);
		for(unsigned y=0; y<dimOut.m_y; y++)
		{
			unsigned yInv = (unsigned)(y*sY + 0.5F);
			for(unsigned x=0; x<dimOut.m_x; x++)
			{
				unsigned xInv = (unsigned)(x*sX + 0.5F);
				lyrOut.CellRef(x, y, z) = lyrIn.CellVal(xInv, yInv, zInv);
			}
		}
	}
}

void ScaleDimLyr::Gen_bilinear(Layer &lyrOut, const Layer &lyrIn) const
{
	Vect3D<unsigned> dimIn  = lyrIn.GetDim();
	Vect3D<unsigned> dimOut = lyrOut.GetDim();

	DATA sX = (DATA)(dimIn.m_x-1) / (dimOut.m_x-1);
	DATA sY = (DATA)(dimIn.m_y-1) / (dimOut.m_y-1);
	DATA sZ = (DATA)(dimIn.m_z-1) / (dimOut.m_z-1);

	for(unsigned z=0; z<dimOut.m_z; z++)
	{
		DATA zInv = z * sZ;
		unsigned zN = (unsigned)zInv;
		unsigned zF = zN + 1;
		if(zF >= dimIn.m_z)
		{
			zF = zN;
		}

		DATA zNRatio = zInv - zN;
		for(unsigned y=0; y<dimOut.m_y; y++)
		{
			DATA yInv = y * sY;
			unsigned yB = (unsigned)yInv;
			unsigned yT = yB + 1; //(unsigned)(yInv + 0.5F);
			if(yT >= dimIn.m_y)
			{
				yT = yB;
			}

			DATA yBRatio = yInv - yB;
			for(unsigned x=0; x<dimOut.m_x; x++)
			{
				DATA xInv = x * sX;
				unsigned xL = (unsigned)xInv;
				unsigned xR = xL + 1; //(unsigned)(xInv + 0.5F);
				if(xR >= dimIn.m_x)
				{
					xR = xL;
				}

				DATA xLRatio = xInv - xL;

				DATA lbn = lyrIn.CellVal(xL, yB, zN);
				DATA rbn = lyrIn.CellVal(xR, yB, zN);
				DATA ltn = lyrIn.CellVal(xL, yT, zN);
				DATA rtn = lyrIn.CellVal(xR, yT, zN);
				DATA lbf = lyrIn.CellVal(xL, yB, zF);
				DATA rbf = lyrIn.CellVal(xR, yB, zF);
				DATA ltf = lyrIn.CellVal(xL, yT, zF);
				DATA rtf = lyrIn.CellVal(xR, yT, zF);
				lyrOut.CellRef(x, y, z) = myMath.Interpolate_linear(
					lbn, rbn, ltn, rtn, 
					lbf, rbf, ltf, rtf, 
					xLRatio, yBRatio, zNRatio);
			}
		}
	}
}

/*
void ScaleDim::Gen_blockAvg(Mtx &mOut, const Mtx &mIn, Mtx *pmBuf) const
{
	Vect2D<unsigned> dimIn  = mIn.GetDim();
	Vect2D<unsigned> dimOut = mOut.GetDim();
	if(dimOut.m_x>=dimIn.m_x || dimOut.m_y>dimIn.m_y)
	{
		Gen_bilinear(mOut, mIn);
		return;
	}
	else {}

	mtxOp.sum.Gen(*pmBuf, mIn);

	DATA sX = (DATA)(dimIn.m_x-1) / (dimOut.m_x-1);
	DATA sY = (DATA)(dimIn.m_y-1) / (dimOut.m_y-1);
	
	unsigned xLen = (unsigned)(sX + 0.5F);
	int xL = -((int)xLen / 2); 
	unsigned xR = xL + xLen - 1;

	unsigned yLen = (unsigned)(sY + 0.5F);
	int yB = -((int)yLen / 2);
	unsigned yT = yB + yLen - 1;

	for(unsigned y=0; y<dimOut.m_y; y++)
	{
		unsigned yCnt = (unsigned)(y*sY + 0.5F);
		int yLocB = (int)yCnt + yB;
		if(yLocB < 0)
		{
			yLocB = 0;
		}
		else {}

		unsigned yLocT = yCnt + yT;
		if(yLocT >= dimIn.m_y)
		{
			yLocT = dimIn.m_y - 1;
		}
		else {}

		for(unsigned x=0; x<dimOut.m_x; x++)
		{
			unsigned xCnt = (unsigned)(x*sX + 0.5F);
			int xLocL = xCnt + xL;
			if(xLocL < 0)
			{
				xLocL = 0;
			}
			else {}

			unsigned xLocR = xCnt + xR;
			if(xLocR >= dimIn.m_x)
			{
				xLocR = dimIn.m_x - 1;
			}
			else {}

			unsigned num = (yLocT - yLocB + 1) * (xLocR - xLocL + 1);
			DATA sum = pmBuf->CellVal(xLocR, yLocT) - pmBuf->CellVal(xLocR, yLocB) - pmBuf->CellVal(xLocL, yLocT) + pmBuf->CellVal(xLocL, yLocB);
			mOut.CellRef(x, y) = sum / num;
		}
	}
}
*/

void ScaleDimLyr::Gen(Layer &lyrOut, const Layer &lyrIn, METHOD m, Layer *pLyrBuf) const
{
	switch(m)
	{
	case NEAREST:
		Gen_nearest(lyrOut, lyrIn);
		break;

	case TRLINEAR:
		Gen_bilinear(lyrOut, lyrIn);
		break;

	/*
	case BLOCK_AVG:
		Gen_blockAvg(mOut, mIn, pmBuf);
		break;
		*/

	default:
		assert(0);
	}
}

//*************************************************************************************************

ScaleDimMtxLyr::ScaleDimMtxLyr()
{}

ScaleDimMtxLyr::~ScaleDimMtxLyr()
{}

void ScaleDimMtxLyr::Gen(Layer &lyrOut, const Layer &lyrIn, METHOD m)
{
	unsigned num = lyrOut.GetDim().m_z;
	for(unsigned i=0; i<num; i++)
	{
		Mtx *pMIn = lyrIn.GetMtx(i);
		Mtx *pMOut = lyrOut.GetMtx(i);

		switch(m)
		{
		case NEAREST:
			mtxOp.scaleDim.Gen(*pMOut, *pMIn, mtxOp.scaleDim.NEAREST);
			break;

		case BILINEAR:
			mtxOp.scaleDim.Gen(*pMOut, *pMIn, mtxOp.scaleDim.BILINEAR);
			break;

		default:
			assert(0);
		}


		delete pMIn;
		delete pMOut;
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

Rgb2Hsv::Rgb2Hsv()
{}

Rgb2Hsv::~Rgb2Hsv()
{}

void Rgb2Hsv::Gen(Layer &lyrHsv, const Layer &lyrRgb, bool bNormH) const
{
	Vect3D<unsigned> dim = lyrRgb.GetDim();
	
	DATA aRgb[3], aHsv[3];
	for(unsigned y=0; y<dim.m_y; y++)
	{
		for(unsigned x=0; x<dim.m_x; x++)
		{
			if(dim.m_z >= 3)
			{
				aRgb[0] = lyrRgb.CellVal(x, y, 0);
				aRgb[1] = lyrRgb.CellVal(x, y, 1);
				aRgb[2] = lyrRgb.CellVal(x, y, 2);
			}
			else if(dim.m_z == 1)
			{
				aRgb[0] = lyrRgb.CellVal(x, y, 0);
				aRgb[1] = aRgb[0];
				aRgb[2] = aRgb[0];
			}

			myMath.RGB2HSV(aHsv, aRgb, bNormH);
			lyrHsv.CellRef(x, y, 0) = aHsv[0];
			lyrHsv.CellRef(x, y, 1) = aHsv[1];
			lyrHsv.CellRef(x, y, 2) = aHsv[2];
		}
	}
}

Gray::Gray()
{}

Gray::~Gray()
{}

void Gray::Gen(Mtx &mtxG, Layer const &lyrIn) const
{
	Vect3D<unsigned> dim = lyrIn.GetDim();

	if(dim.m_z == 1)
	{
		for(unsigned y=0; y<dim.m_y; y++)
		{
			for(unsigned x=0; x<dim.m_x; x++)
			{
				mtxG.CellRef(x, y) = lyrIn.CellVal(x, y, 0);
			}
		}
	}
	else if(dim.m_z==3 || dim.m_z==4)
	{
		for(unsigned y=0; y<dim.m_y; y++)
		{
			for(unsigned x=0; x<dim.m_x; x++)
			{
				DATA r = lyrIn.CellVal(x, y, 0);
				DATA g = lyrIn.CellVal(x, y, 1);
				DATA b = lyrIn.CellVal(x, y, 2);
				mtxG.CellRef(x, y) = 0.3F*r + 0.59F*g + 0.11F*b;
			}
		}
	}
	else assert(0);
}

//*************************************************************************************************
//
//*************************************************************************************************

HisLyr::HisLyr()
{
	for(unsigned i=0; i<3; i++)
	{
		m_aMaxV[i] = 1.F;
		m_aMinV[i] = 0;
		m_aVLen[i] = m_aMaxV[i] - m_aMinV[i];
		m_aInvVLen[i] = 1.F / m_aVLen[i];
	}
}

HisLyr::~HisLyr()
{}

//***********************************************
/*
void HisLyr::SetBound(DATA aMaxV[], DATA aMinV[])
{
	for(unsigned i=0; i<3; i++)
	{
		m_aMaxV[i] = aMaxV[i];
		m_aMinV[i] = aMinV[i];
		m_aVLen[i] = m_aMaxV[i] - m_aMinV[i];
		m_aInvVLen[i] = 1.F / m_aVLen[i];
	}
}

void HisLyr::Add(Layer &lyrHis, DATA aVal[]) const
{
	Vect3D<unsigned> dim = lyrHis.GetDim();
	unsigned aHisMax[] = {dim.m_x-1, dim.m_y-1, dim.m_z-1};

	int aLoc[3];
	for(unsigned c=0; c<3; c++)
	{
		aLoc[c] = (int)((aVal[c] - m_aMinV[c]) * m_aInvVLen[c] * aHisMax[c]);
		if(aLoc[c] < 0)
		{
			aLoc[c] = 0;
		}
		else if(aLoc[c] > (int)aHisMax[c])
		{
			aLoc[c] = aHisMax[c];
		}
		else {}
	}
	lyrHis.CellRef(aLoc[0], aLoc[1], aLoc[2]) = lyrHis.CellVal(aLoc[0], aLoc[1], aLoc[2]) + 1.F;
} 
*/
//***********************************************

void HisLyr::Gen(Mtx &mtxHis, const Layer &lyrIn, DATA minV, DATA maxV) const
{
	Vect2D<unsigned> dimHis = mtxHis.GetDim();
	assert(dimHis.m_x == 1);
	for(unsigned y=0; y<dimHis.m_y; y++)
	{
		mtxHis.CellRef(0, y) = 0;
	}

	assert(maxV > minV);
	DATA vLen = maxV - minV;
	DATA vLenRatio = 1.F / vLen;
	Vect3D<unsigned> dimIn = lyrIn.GetDim();
	for(unsigned z=0; z<dimIn.m_z; z++)
	{
		for(unsigned y=0; y<dimIn.m_y; y++)
		{
			for(unsigned x=0; x<dimIn.m_x; x++)
			{
				int hisLoc = (int)((lyrIn.CellVal(x, y, z)-minV) * vLenRatio * (dimHis.m_y-1) + 0.5F);
				if(hisLoc<0 || hisLoc>=(int)dimHis.m_y)
				{
					continue;
				}
				else {}
				/*
				if(hisLoc < 0)
				{
					hisLoc = 0;
				}
				else if(hisLoc >= (int)dimHis.m_y)
				{
					hisLoc = dimHis.m_y - 1;
				}
				else {}
				*/
				mtxHis.CellRef(0, hisLoc) = mtxHis.CellVal(0, hisLoc) + 1.F;
			}
		}
	}
}

void HisLyr::Gen(Mtx &mtxHis, const Layer &lyrX, const Layer &lyrY, 
		DATA minX, DATA maxX, DATA minY, DATA maxY) const
{
	Vect2D<unsigned> dimHis = mtxHis.GetDim();
	mtxOp.zero.Gen(mtxHis);

	assert(maxX > minX);
	assert(maxY > minY);
	DATA xLen = maxX - minX;	DATA xLenRatio = 1.F / xLen;
	DATA yLen = maxY - minY;	DATA yLenRatio = 1.F / yLen;

	Vect3D<unsigned> dimX = lyrX.GetDim();
	Vect3D<unsigned> dimY = lyrY.GetDim();
	assert(dimX.m_x==dimY.m_x && dimX.m_y==dimY.m_y && dimX.m_z==dimY.m_z);
	for(unsigned z=0; z<dimX.m_z; z++)
	{
		for(unsigned y=0; y<dimX.m_y; y++)
		{
			for(unsigned x=0; x<dimX.m_x; x++)
			{
				int hisX = (int)((lyrX.CellVal(x, y, z)-minX) * xLenRatio * (dimHis.m_x-1) + 0.5F);
				int hisY = (int)((lyrY.CellVal(x, y, z)-minY) * yLenRatio * (dimHis.m_y-1) + 0.5F);
				if(hisX<0 || hisX>=(int)dimHis.m_x ||
				   hisY<0 || hisY>=(int)dimHis.m_y)
				{
					continue;
				}
				else {}

				mtxHis.CellRef(hisX, hisY) = mtxHis.CellVal(hisX, hisY) + 1;
			}
		}
	}
}

/*
void HisLyr::Gen(Layer &lyrHis, const Layer &lyrIn, DATA aMaxV[], DATA aMinV[]) const
{
	assert(aMaxV[0]>aMinV[0] && aMaxV[1]>aMinV[1] && aMaxV[2]>aMinV[2]);

	Vect3D<unsigned> dimIn = lyrIn.GetDim();
	assert(dimIn.m_z == 3);

	DATA aVLen[3];
	DATA aInvVLen[3];
	for(unsigned i=0; i<3; i++)
	{
		aVLen[i] = aMaxV[i] - aMinV[i];
		aInvVLen[i] = 1.F / aVLen[i];
	}

	Vect3D<unsigned> dimHis = lyrHis.GetDim();
	unsigned aHisMax[] = {dimHis.m_x-1, dimHis.m_y-1, dimHis.m_z-1};

	lyrOp.zero.Gen(lyrHis);
	for(unsigned y=0; y<dimIn.m_y; y++)
	{
		for(unsigned x=0; x<dimIn.m_x; x++)
		{
			DATA aV[3];
			aV[0] = lyrIn.CellVal(x, y, 0);
			aV[1] = lyrIn.CellVal(x, y, 1);
			aV[2] = lyrIn.CellVal(x, y, 2);

			int aLoc[3];
			for(unsigned c=0; c<3; c++)
			{
				aLoc[c] = (int)((aV[c] - aMinV[c]) * aInvVLen[c] * aHisMax[c]);
				if(aLoc[c] < 0)
				{
					aLoc[c] = 0;
				}
				else if(aLoc[c] > (int)aHisMax[c])
				{
					aLoc[c] = aHisMax[c];
				}
				else {}
			}
			lyrHis.CellRef(aLoc[0], aLoc[1], aLoc[2]) = lyrHis.CellVal(aLoc[0], aLoc[1], aLoc[2]) + 1.F;
		}
	}
}
*/
//*************************************************************************************************
//
//*************************************************************************************************

MulLyr::MulLyr()
{}

MulLyr::~MulLyr()
{}

void MulLyr::Gen(Layer &lyr, DATA scl) const
{
	Vect3D<unsigned> dim = lyr.GetDim();
	for(unsigned i=0; i<dim.m_z; i++)
	{
		Mtx *pM = lyr.GetMtx(i);
		mtxOp.mul.Gen(*pM, scl);
	}
}

SubLyr::SubLyr()
{}

SubLyr::~SubLyr()
{}

void SubLyr::Gen(Layer &lyrA, const Layer &lyrB)
{
	Vect3D<unsigned> dimA = lyrA.GetDim();
	Vect3D<unsigned> dimB = lyrB.GetDim();
	assert(dimA.m_x==dimB.m_x && dimA.m_y==dimB.m_y && dimA.m_z==dimB.m_z);

	for(unsigned z=0; z<dimA.m_z; z++)
	{
		for(unsigned y=0; y<dimA.m_y; y++)
		{
			for(unsigned x=0; x<dimA.m_x; x++)
			{
				lyrA.CellRef(x, y, z) = lyrA.CellVal(x, y, z) - lyrB.CellVal(x, y, z);
			}
		}
	}
}

void SubLyr::Gen(Layer &lyr, DATA v)
{
	Vect3D<unsigned> dim = lyr.GetDim();

	for(unsigned z=0; z<dim.m_z; z++)
	{
		for(unsigned y=0; y<dim.m_y; y++)
		{
			for(unsigned x=0; x<dim.m_x; x++)
			{
				lyr.CellRef(x, y, z) = lyr.CellVal(x, y, z) - v;
			}
		}
	}
}

AddLyr::AddLyr()
{}

AddLyr::~AddLyr()
{}

void AddLyr::Gen(Layer &lyr, DATA v)
{
	Vect3D<unsigned> dim = lyr.GetDim();

	for(unsigned z=0; z<dim.m_z; z++)
	{
		for(unsigned y=0; y<dim.m_y; y++)
		{
			for(unsigned x=0; x<dim.m_x; x++)
			{
				lyr.CellRef(x, y, z) = lyr.CellVal(x, y, z) + v;
			}
		}
	}
}

CellMultiplyLyr::CellMultiplyLyr()
{}

CellMultiplyLyr::~CellMultiplyLyr()
{}

void CellMultiplyLyr::Gen(CplxLyr &clTo, const CplxLyr &clA, const CplxLyr &clB)
{
	Vect3D<unsigned> dimA  = clA.GetDim();
	Vect3D<unsigned> dimB  = clB.GetDim();
	Vect3D<unsigned> dimTo = clTo.GetDim();
	assert(dimA.m_x==dimB.m_x && dimA.m_y==dimB.m_y && dimA.m_z==dimB.m_z);
	assert(dimTo.m_x==dimA.m_x && dimTo.m_y==dimA.m_y && dimTo.m_z==dimA.m_z);

	for(unsigned z=0; z<dimA.m_z; z++)
	{
		for(unsigned y=0; y<dimA.m_y; y++)
		{
			for(unsigned x=0; x<dimA.m_x; x++)
			{
				Vect2D<DATA> valA = clA.CellVal(x, y, z);
				Vect2D<DATA> valB = clB.CellVal(x, y, z);
				Vect2D<DATA&> valTo = clTo.CellRef(x, y, z);
				valTo.m_x = valA.m_x*valB.m_x - valA.m_y*valB.m_y;
				valTo.m_y = valA.m_x*valB.m_y + valA.m_y*valB.m_x;
			}
		}
	}
}

void CellMultiplyLyr::Gen(Layer &lyrTo, const Layer &lyrA, const Layer &lyrB)
{
	Vect3D<unsigned> dimA  = lyrA.GetDim();
	Vect3D<unsigned> dimB  = lyrB.GetDim();
	Vect3D<unsigned> dimTo = lyrTo.GetDim();
	assert(dimA.m_x==dimB.m_x && dimA.m_y==dimB.m_y && dimA.m_z==dimB.m_z);
	assert(dimTo.m_x==dimA.m_x && dimTo.m_y==dimA.m_y && dimTo.m_z==dimA.m_z);

	for(unsigned z=0; z<dimA.m_z; z++)
	{
		for(unsigned y=0; y<dimA.m_y; y++)
		{
			for(unsigned x=0; x<dimA.m_x; x++)
			{
				lyrTo.CellRef(x, y, z) = lyrA.CellVal(x, y, z) * lyrB.CellVal(x, y, z);
			}
		}
	}
}


//*************************************************************************************************
//
//*************************************************************************************************

AvgLyr::AvgLyr()
{}

AvgLyr::~AvgLyr()
{}

DATA AvgLyr::Gen(const Layer &lyr) const
{
	Vect3D<unsigned> dim = lyr.GetDim();
	unsigned size = dim.m_x * dim.m_y * dim.m_z;

	DATA sum = 0;
	for(unsigned i=0; i<dim.m_z; i++)
	{
		Mtx *pM = lyr.GetMtx(i);
		sum += mtxOp.sum.Gen(*pM);
	}
	return sum/size;
}

VarLyr::VarLyr()
{}

VarLyr::~VarLyr()
{}

DATA VarLyr::Gen(const Layer &lyr) const
{
	Vect3D<unsigned> dim = lyr.GetDim();
	unsigned size = dim.m_x * dim.m_y * dim.m_z;

	DATA avg = lyrOp.avg.Gen(lyr);

	DATA var = 0;
	for(unsigned z=0; z<dim.m_z; z++)
	{
		for(unsigned y=0; y<dim.m_y; y++)
		{
			for(unsigned x=0; x<dim.m_x; x++)
			{
				DATA diff = lyr.CellVal(x, y, z) - avg;
				var += diff * diff;
			}
		}
	}
	return var/size;
}

RngLyr::RngLyr()
{}

RngLyr::~RngLyr()
{}

void RngLyr::Gen(DATA &vMin, DATA &vMax, const Layer &lyr) const
{
	Vect3D<unsigned> dim = lyr.GetDim();
	vMin = 1e10;
	vMax = -1e10;
	for(unsigned z=0; z<dim.m_z; z++)
	{
		for(unsigned y=0; y<dim.m_y; y++)
		{
			for(unsigned x=0; x<dim.m_x; x++)
			{
				DATA v = lyr.CellVal(x, y, z);
				if(v < vMin)
				{
					vMin = v;
				}
				else {}
				if(v > vMax)
				{
					vMax = v;
				}
				else {}
			}
		}
	}
}

//*************************************************************************************************
//
//*************************************************************************************************

FFTLyr::FFTLyr()
	:m_pFrom(0), m_pTo(0),   m_pCplxFrom(0), m_pCplxTo(0) 
	,m_fSize(0), m_tSize(0), m_cfSize(0),    m_ctSize(0)
	,m_tCase(R2C), m_planXDim(0), m_planYDim(0), m_planZDim(0), m_plan(0)
{}

FFTLyr::~FFTLyr()
{
	fftw_destroy_plan(m_plan);
	fftw_free(m_pFrom);	
	fftw_free(m_pTo);		
	fftw_free(m_pCplxFrom);
	fftw_free(m_pCplxTo);
}

void FFTLyr::NewFMem(unsigned size)
{
	if(m_fSize >= size)
	{
		return;
	}
	else
	{
		fftw_free(m_pFrom);
		m_fSize = size;
		m_pFrom = (double*)fftw_malloc(sizeof(double) * m_fSize);
	}
}
void FFTLyr::NewTMem(unsigned size)
{
	if(m_tSize >= size)
	{
		return;
	}
	else
	{
		fftw_free(m_pTo);
		m_tSize = size;
		m_pTo = (double*)fftw_malloc(sizeof(double) * m_tSize);
	}
}
void FFTLyr::NewCFMem(unsigned size)
{
	if(m_cfSize >= size)
	{
		return;
	}
	else
	{
		fftw_free(m_pCplxFrom);
		m_cfSize = size;
		m_pCplxFrom = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * m_cfSize);
	}
}
void FFTLyr::NewCTMem(unsigned size)
{
	if(m_ctSize >= size)
	{
		return;
	}	
	else
	{
		fftw_free(m_pCplxTo);
		m_ctSize = size;
		m_pCplxTo = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * m_ctSize);
	}
}

void FFTLyr::SetPlan(unsigned xDim, unsigned yDim, unsigned zDim, T_CASE tCase)
{
	if(xDim==m_planXDim && yDim==m_planYDim && zDim==m_planZDim && tCase ==m_tCase)	return;
	else {}

	fftw_destroy_plan(m_plan);

	unsigned rSize = xDim * yDim * zDim;
	unsigned cSize = (xDim/2 + 1) * yDim * zDim;  
	if(tCase == C2C)
	{
		NewCFMem(rSize);
		NewCTMem(rSize);
		m_plan = fftw_plan_dft_3d(zDim, yDim, xDim, m_pCplxFrom, m_pCplxTo, FFTW_FORWARD, FFTW_ESTIMATE);
	}
	else if(tCase == C2R)
	{
		NewCFMem(cSize);
		NewTMem (rSize);
		m_plan = fftw_plan_dft_c2r_3d(zDim, yDim, xDim, m_pCplxFrom, m_pTo, FFTW_ESTIMATE);
	}
	else if(tCase == R2C)
	{
		NewFMem (rSize);
		NewCTMem(cSize);
		m_plan = fftw_plan_dft_r2c_3d(zDim, yDim, xDim, m_pFrom, m_pCplxTo, FFTW_ESTIMATE);
	}
	else assert(0);
	m_tCase = tCase;
	m_planXDim = xDim;
	m_planYDim = yDim;
	m_planZDim = zDim;
}

void FFTLyr::Gen(CplxLyr &cLyrTo, Layer &lyrFrom, bool bDividNum)
{
	Vect3D<unsigned> dimFrom = lyrFrom.GetDim(); 
	Vect3D<unsigned> dimTo   = cLyrTo.GetDim();
	unsigned planCXDim = dimFrom.m_x/2 + 1;
	unsigned planCYDim = dimFrom.m_y;
	unsigned planCZDim = dimFrom.m_z;
	assert(dimTo.m_x==planCXDim && dimTo.m_y==planCYDim && dimTo.m_z==planCZDim);

	SetPlan(dimFrom.m_x, dimFrom.m_y, dimFrom.m_z, R2C);

	unsigned idx = 0;
	for(unsigned z=0; z<m_planZDim; z++)
	{
		for(unsigned y=0; y<m_planYDim; y++)
		{
			for(unsigned x=0; x<m_planXDim; x++)
			{
				m_pFrom[idx] = lyrFrom.CellVal(x, y, z);
				++idx;
			}
		}
	}

	fftw_execute(m_plan);

	idx = 0;
	for(unsigned z=0; z<planCZDim; z++)
	{
		for(unsigned y=0; y<planCYDim; y++)
		{
			for(unsigned x=0; x<planCXDim; x++)
			{
				Vect2D<DATA&> valTo = cLyrTo.CellRef(x, y, z);
				valTo.m_x = m_pCplxTo[idx][0] ;
				valTo.m_y = m_pCplxTo[idx][1];
				++idx;
			}
		}
	}

	if(bDividNum)
	{
		unsigned size = m_planXDim * m_planYDim * m_planZDim;
		DATA invN = 1.F / size;
		for(unsigned z=0; z<planCZDim; z++)
		{
			for(unsigned y=0; y<planCYDim; y++)
			{
				for(unsigned x=0; x<planCXDim; x++)
				{
					Vect2D<DATA&> val = cLyrTo.CellRef(x, y, z);
					val.m_x *= invN;
					val.m_y *= invN;
				}
			}
		}
	}
	else {}
}

void FFTLyr::Gen(Layer &lyrTo, CplxLyr &cLyrFrom, bool bDividNum)
{
	Vect3D<unsigned> dimFrom = cLyrFrom.GetDim(); 
	Vect3D<unsigned> dimTo   = lyrTo.GetDim();
	unsigned planCXDim = dimTo.m_x/2 + 1;
	unsigned planCYDim = dimTo.m_y;
	unsigned planCZDim = dimTo.m_z;
	assert(dimFrom.m_x==planCXDim && dimFrom.m_y==planCYDim && dimFrom.m_z==planCZDim) ;

	SetPlan(dimTo.m_x, dimTo.m_y, dimTo.m_z, C2R);

	unsigned idx = 0;
	for(unsigned z=0; z<planCZDim; z++)
	{
		for(unsigned y=0; y<planCYDim; y++)
		{
			for(unsigned x=0; x<planCXDim; x++)
			{
				Vect2D<DATA> valFrom = cLyrFrom.CellVal(x, y, z);
				m_pCplxFrom[idx][0] = valFrom.m_x;
				m_pCplxFrom[idx][1] = valFrom.m_y;
				++idx;
			}
		}
	}

	fftw_execute(m_plan);

	idx = 0;
	for(unsigned z=0; z<m_planZDim; z++)
	{
		for(unsigned y=0; y<m_planYDim; y++)
		{
			for(unsigned x=0; x<m_planXDim; x++)
			{
				lyrTo.CellRef(x, y, z) = m_pTo[idx];
				++idx;
			}
		}
	}

	if(bDividNum)
	{
		unsigned size = m_planXDim * m_planYDim * m_planZDim;
		DATA invN = 1.F / size;
		for(unsigned z=0; z<m_planZDim; z++)
		{
			for(unsigned y=0; y<m_planYDim; y++)
			{
				for(unsigned x=0; x<m_planXDim; x++)
				{
					lyrTo.CellRef(x, y, z) *= invN;
				}
			}
		}
	}
	else {}
}
//*************************************************************************************************

ConvLyr::ConvLyr()
{}

ConvLyr::~ConvLyr()
{}

void ConvLyr::Gen(Layer &lyrBase, const Layer &lyrKerl)
{
	Vect3D<unsigned> dimBase = lyrBase.GetDim();
	Vect3D<unsigned> dimKerl = lyrKerl.GetDim();
	Vect3D<unsigned> dimFft(dimBase.m_x + dimKerl.m_x - 1, 
						    dimBase.m_y + dimKerl.m_y - 1,
							dimBase.m_z + dimKerl.m_z - 1);

	cout << "convolve " << dimBase.m_x << "x" << dimBase.m_y << "x" << dimBase.m_z  
		 << " with "    << dimKerl.m_x << "x" << dimKerl.m_y << "x" << dimKerl.m_z << endl;
	
	unsigned cXDim = dimFft.m_x/2 + 1;
	unsigned cYDim = dimFft.m_y;
	unsigned cZDim = dimFft.m_z;
	Layer lBasePad(dimFft.m_x, dimFft.m_y, dimFft.m_z);		lyrOp.zero.Gen(lBasePad);
	Layer lKerlPad(dimFft.m_x, dimFft.m_y, dimFft.m_z);		lyrOp.zero.Gen(lKerlPad);
	Layer lBaseV(lBasePad, 0, 0, 0, dimBase.m_x, dimBase.m_y, dimBase.m_z);		lBaseV.CopyFrom(lyrBase);
	Layer lKerlV(lKerlPad, 0, 0, 0, dimKerl.m_x, dimKerl.m_y, dimKerl.m_z);		lKerlV.CopyFrom(lyrKerl);

	CplxLyr clBase(cXDim, cYDim, cZDim);	m_fft.Gen(clBase, lBasePad);
	CplxLyr clKerl(cXDim, cYDim, cZDim);	m_fft.Gen(clKerl, lKerlPad);

	CplxLyr clFft(cXDim, cYDim, cZDim);
	lyrOp.cellX.Gen(clFft, clBase, clKerl);

	Layer lFft(dimFft.m_x, dimFft.m_y, dimFft.m_z);
	m_ifft.Gen(lFft, clFft, true);
	
	unsigned kerlCntX = dimKerl.m_x / 2;
	unsigned kerlCntY = dimKerl.m_y / 2;
	unsigned kerlCntZ = dimKerl.m_z / 2;
	lyrBase.CopyFrom(Layer(lFft, kerlCntX, kerlCntY, kerlCntZ, dimBase.m_x, dimBase.m_y, dimBase.m_z));

	cout << "convolution end" << endl << endl;
}

/*
void ConvLyr::Gen(Layer &lyrOut, const Layer &lyrA, const Layer &lyrB)
{
	Vect3D<unsigned> dimOut = lyrOut.GetDim();
	Vect3D<unsigned> dimA = lyrA.GetDim();
	Vect3D<unsigned> dimB = lyrB.GetDim();
	assert(dimOut.m_x==dimA.m_x && dimOut.m_y==dimA.m_y && dimOut.m_z==dimA.m_z);

	for(unsigned z=0; z<dimOut.m_z; z++)
	{
		cout << z << " ";
		int zN = z - dimB.m_z/2;
		int zF = zN + dimB.m_z - 1;
		for(unsigned y=0; y<dimOut.m_y; y++)
		{
			int yB = y - dimB.m_y/2;
			int yT = yB + dimB.m_y - 1;
			for(unsigned x=0; x<dimOut.m_x; x++)
			{
				int xL = x - dimB.m_x/2;
				int xR = xL + dimB.m_x - 1;

				DATA sum = 0;
				for(int zz=zN; zz<=zF; zz++)
				{
					if(zz<0 || zz>=(int)dimOut.m_z)
					{
						continue;
					}
					else {}
					unsigned locBZ = zz - zN;

					for(int yy=yB; yy<=yT; yy++)
					{
						if(yy<0 || yy>=(int)dimOut.m_y)
						{
							continue;
						}
						else {}
						unsigned locBY = yy - yB;

						for(int xx=xL; xx<=xR; xx++)
						{
							if(xx<0 || xx>=(int)dimOut.m_x)
							{
								continue;
							}
							else {}
							unsigned locBX = xx - xL;
						
							sum += lyrA.CellVal(xx, yy, zz) * lyrB.CellVal(locBX, locBY, locBZ);
							//cout << lyrB.CellVal(locBX, locBY, locBZ) << " ";
						} // xx
					} // yy
				} // zz

				lyrOut.CellRef(x, y, z) = sum;
				//cout << sum << "";
			} // x
		} // y
	} // z
	cout << "convolution ok." << endl;
}
*/

DoGLyr::DoGLyr()
{}

DoGLyr::~DoGLyr()
{}

void DoGLyr::Gen(Layer &lyrIn, Layer &lyrGauss, DATA s, DATA sScl)
{
	Vect3D<unsigned> dimG = lyrGauss.GetDim();
	Layer lyrGTmp(dimG.m_x, dimG.m_y, dimG.m_z);
	lyrOp.Gauss3D.Gen(lyrGTmp, s*sScl);
	lyrOp.mul.Gen(lyrGTmp, 0.98F);

	lyrOp.Gauss3D.Gen(lyrGauss, s);
	lyrOp.sub.Gen(lyrGauss, lyrGTmp);

	lyrOp.conv.Gen(lyrIn, lyrGauss);
}

void DoGLyr::Gen(VolumeData &volIn, Layer &lyrGauss, DATA s, DATA sScl)
{
	VolumeData volGauss(lyrGauss, volIn.GetStep());
	VolumeData volGTmp(volGauss.GetDim(), volGauss.GetStep());

	lyrOp.Gauss3D.Gen(volGTmp, s*sScl);
	//lyrOp.mul.Gen(volGTmp.GetLyrRef(), 0.98F);

	lyrOp.Gauss3D.Gen(volGauss, s);
	lyrOp.sub.Gen(volGauss.GetLyrRef(), volGTmp.GetLyrRef());

	lyrOp.conv.Gen(volIn.GetLyrRef(), volGauss.GetLyrRef());
}

//***********************************************

SplattingLyr::SplattingLyr()
	: m_pLyrCount(0)
{}

SplattingLyr::~SplattingLyr()
{
	delete m_pLyrCount;
}

void SplattingLyr::New(Layer &lyrIn)
{
	Vect3D<unsigned> dimIn = lyrIn.GetDim();

	if(!m_pLyrCount)
	{
		delete m_pLyrCount;
		m_pLyrCount = new Layer(dimIn.m_x, dimIn.m_y, dimIn.m_z);
	}
	else 
	{
		Vect3D<unsigned> dimC = m_pLyrCount->GetDim();
		if(dimC.m_x<dimIn.m_x || 
		   dimC.m_y<dimIn.m_y || 
		   dimC.m_z<dimIn.m_z)
		{
			delete m_pLyrCount;
			m_pLyrCount = new Layer(dimIn.m_x, dimIn.m_y, dimIn.m_z);
		}
		else {}
	}

	lyrOp.zero.Gen(lyrIn);
	lyrOp.zero.Gen(*m_pLyrCount);
}

void SplattingLyr::AddPoint(Layer &lyrIn, DATA val, const Layer &lyrW, 
	unsigned x, unsigned y, unsigned z)
{
	Vect3D<unsigned> dimW = lyrW.GetDim();
 
	for(unsigned zz=0; zz<dimW.m_z; zz++)
	{
		int zLoc = (int)z + zz - dimW.m_z/2;
		if(!lyrIn.IsCInside(zLoc))	continue;	else {}
		for(unsigned yy=0; yy<dimW.m_y; yy++)
		{
			int yLoc = (int)y + yy - dimW.m_y/2;
			if(!lyrIn.IsYInside(yLoc))	continue;	else {}
			for(unsigned xx=0; xx<dimW.m_x; xx++)
			{
				int xLoc = (int)x + xx - dimW.m_x/2;
				if(!lyrIn.IsXInside(xLoc))	continue;	else {}

				//cout << xLoc << " " << yLoc << " " << zLoc << "aa" << endl;
				lyrIn.CellRef(xLoc, yLoc, zLoc) = 
					lyrIn.CellVal(xLoc, yLoc, zLoc) + val*lyrW.CellVal(xx, yy, zz);
				m_pLyrCount->CellRef(xLoc, yLoc, zLoc) = 
					m_pLyrCount->CellVal(xLoc, yLoc, zLoc) + lyrW.CellVal(xx, yy, zz);
			}
		}
	}
}

void SplattingLyr::Gen(Layer &lyrIn)
{
	Vect3D<unsigned> dimIn = lyrIn.GetDim();
	for(unsigned z=0; z<dimIn.m_z; z++)
	{
		for(unsigned y=0; y<dimIn.m_y; y++)
		{
			for(unsigned x=0; x<dimIn.m_x; x++)
			{
				if(m_pLyrCount->CellVal(x, y, z) > 0)
				{
					lyrIn.CellRef(x, y, z) = 
						lyrIn.CellVal(x, y, z) / m_pLyrCount->CellVal(x, y, z);
				}
				else
				{
					lyrIn.CellRef(x, y, z) = 0;
				}
			}
		}
	}
}