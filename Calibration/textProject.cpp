#include "Calibration.h"
#include "../include/ICP/MIcpPointToPlane.h"

int main()
{
	//======================内部角点的个数
#ifdef calibration
	const std::string inFileName("data/test.txt");
	cv::Size broadSize = cv::Size(9,6);
	cv::Size roadSize2 = cv::Size(50,50);
	CCalibration a;
	//a.GetcameraData();
	CCalibration cal(inFileName,broadSize,roadSize2);
	cal.SingleCalibration();
#else
	CCalibration cal;
	cal.Icp();
#endif // calibration


		system("pause");
		return 0;


}