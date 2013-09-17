#include "define.h"

const DATA DATA_NAN = numeric_limits<DATA>::quiet_NaN();

const string WORK_DIR			= "D:\\working\\visionComputing\\siftImplement";
const string TEST_DATA_DIR		= WORK_DIR + "\\testData";
const string OUTPUT_DIR			= WORK_DIR + "\\output";

//*************************************************************************************************

const DATA PI		= (DATA)3.14159265358979323846;
const DATA PI_2		= PI * 2.F;
const DATA R2D		= 180.F / PI;
const DATA D2R		= PI / 180.F;
const DATA LOG2		= log(2.F);
const DATA INV_LOG2 = 1.F / LOG2;
const DATA INV_255	= 1.F / 255.F;
const DATA INV_360	= 1.F / 360.F;

//*************************************************************************************************