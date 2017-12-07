#include "Calibration.h"
#include <fstream>
void CUndistort::getMatrixParameters()
{
	cv::FileStorage ofs(srcImgPath , cv::FileStorage::READ);
	ofs["cameraMatrix"] >> cameraMatrix;
	ofs["distCoeffs"]>>distCoeffs;
	ofs.release();
}
void CUndistort::initUndistortRectifyMap()
{
	getMatrixParameters();
	cv::Mat srcImg = cv::imread(srcImgPath+"text",0);
	cv::pyrDown(srcImg,srcImg);
	cv::Mat R=cv::Mat::eye(cv::Size(3, 3),CV_32FC1);
	cv::Mat mapx, mapy;
	cv::initUndistortRectifyMap(cameraMatrix,distCoeffs,R,cameraMatrix,cv::Size(srcImg.cols,srcImg.rows),CV_32FC1,mapx,mapy);
	cv::remap(srcImg,srcImg,mapx,mapy,CV_INTER_LINEAR);
	cv::imshow("dst",srcImg);
	std::cout << "ÊäÈë°´¼ü¼ÌÐø" << std::endl;
	cv::waitKey(0);
}