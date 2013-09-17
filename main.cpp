#include "denseMatrix.h"
#include "imageIO.h"
#include "matrixOperation.h"
#include "sift.h"

//#include <gsl/gsl_matrix.h>

DATA ComputeDis(Keypoint &kA, Keypoint &kB)
{
	DATA dis = 0;
	for(unsigned y=0; y<4; y++)
	{
		for(unsigned x=0; x<4; x++)
		{
			for(unsigned d=0; d<8; d++)
			{
				DATA dif = kA.m_aaaDspt[y][x][d] - kB.m_aaaDspt[y][x][d];
				dis += dif * dif;
			}
		}
	}
	dis = sqrt(dis);
	return dis;
}

void MarkKp(MyImg &img, Keypoint &kp, bool br)
{
	Vect3D<unsigned> dim = img.GetDim();
	
	unsigned u = 15;
	int xCnt = (unsigned)(kp.m_x + 0.5F);
	int yCnt = (unsigned)(kp.m_y + 0.5f);
	if(!br)
	{
		xCnt += (int)(u * cos(kp.m_ornt));
		yCnt += (int)(u * sin(kp.m_ornt));
	}
	else {}
	
	unsigned w = 10;
	for(unsigned i=0; i<w; i++)
	{
		int aXT[] = {-1,  1, 1, -1};
		int aYT[] = {-1, -1, 1,  1};
		for(unsigned j=0; j<4; j++)
		{
			int x = xCnt + aXT[j]*i;
			int y = yCnt + aYT[j]*i;
			if(x<0 || x>=(int)dim.m_x) continue;		else {}
			if(y<0 || y>=(int)dim.m_y) continue;		else {}

			if(br)
			{
				img.CellRef(x, y, 0) = 255.F;
				img.CellRef(x, y, 1) = 0;
				img.CellRef(x, y, 2) = 0;
			}
			else
			{
				img.CellRef(x, y, 0) = 0;
				img.CellRef(x, y, 1) = 255;
				img.CellRef(x, y, 2) = 100.;
			}
		}
	}
}

void main()
{
	//Mtx2D::Test();
	//ImgIO::Test();
	//Gauss2D::Test();
	//FFT::Test();
	//Conv::Test();
	//Inv::Test();
	//Mul::Test();
	//Sift::Test();
	
	/*
	string tName  = TEST_DATA_DIR + "\\t.bmp"; 
	MyImg *pTImg  = imgIO.Read(tName);
	Mtx *pTM  = pTImg->GetMtx(0);

	Sift sift;
	Mtx mT(pTM->GetDim());		mT.CopyFrom(*pTM);
	vector<Keypoint>kpT = sift.Gen(mT,  4);
	for(unsigned i=0; i<kpT.size(); i++)
	{
		MarkKp(*pTImg, kpT[i], true);			MarkKp(*pTImg, kpT[i], false);	
	}
	string outTName = OUTPUT_DIR + "\\outT.bmp";
	imgIO.Write(outTName, *pTImg);
	*/
	
	string markName  = TEST_DATA_DIR + "\\basmati.bmp"; 
	string sceneName = TEST_DATA_DIR + "\\boxes.bmp";
	MyImg *pMarkImg  = imgIO.Read(markName);
	MyImg *pSceneImg = imgIO.Read(sceneName);
	Mtx *pMarkM  = pMarkImg->GetMtx(0);
	Mtx *pSceneM = pSceneImg->GetMtx(0);

	Sift sift;
	Mtx mMark(pMarkM->GetDim());		mMark.CopyFrom(*pMarkM);
	Mtx mScene(pSceneM->GetDim());		mScene.CopyFrom(*pSceneM);
	vector<Keypoint>kpMark  = sift.Gen(mMark,  4);
	vector<Keypoint>kpScene = sift.Gen(mScene, 4);

	for(unsigned i=0; i<kpMark.size(); i++)
	{
		DATA fDis = 1e10;
		DATA sDis = 1e10;
		unsigned fId = 0;
		unsigned sId = 0;
		for(unsigned j=0; j<kpScene.size(); j++)
		{
			DATA dis = ComputeDis(kpMark[i], kpScene[j]);
			
			if(dis < fDis)
			{
				sId = fId;		sDis = fDis;
				fId = j;		fDis = dis;
			}
			else if(dis < sDis)
			{
				sId = j;		sDis = dis;
			}
			else {}
		}

		DATA ratioDis = fDis / sDis;
		if(ratioDis < 0.8F)
		{
			MarkKp(*pSceneImg, kpScene[fId], true);			//MarkKp(*pSceneImg, kpScene[fId], false);	
			MarkKp(*pMarkImg,  kpMark[i], true);			//MarkKp(*pMarkImg,  kpMark[i], false);
		}	
		else {}
	}

	string outSceneName = OUTPUT_DIR + "\\outScene.bmp";
	string outMarkName  = OUTPUT_DIR + "\\outMark.bmp";
	imgIO.Write(outSceneName, *pSceneImg);
	imgIO.Write(outMarkName,  *pMarkImg);
	
	

	system("pause");
}