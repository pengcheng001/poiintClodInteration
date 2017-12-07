#ifndef  _CALIBRATION_H_
#define  _CALIBRATION_H_

#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
//定义设备数量
#define  camereNum 2
////启动设备IR图，按R截取读片，并之后启动标定程序 
//#define IR
//#define calibration 
//启动设备重新获取点云yml数据，如果data下已经有数据，则不必重新生成并保持开关关闭
//#define ifNeedgetPointData 
class CCalibration
{
public:
	CCalibration(){};
	CCalibration(std::string SrcimgPath ,const cv::Size broadSize,const cv::Size realBroadSize)  //初始化类
	{  this->SrcimgPath = SrcimgPath; this->broadSize = broadSize;this->realBroadSize = realBroadSize;}
	void  findChessboardCorners(std::vector<cv::Point2f> &cornersPointBuf,
		std::vector<std::vector<cv::Point2f>> &all_cornersPointBuf) ;
	void SingleCalibration();
	void CalUnionCalibrationParameters(std::vector<std::vector<cv::Mat>> &rvecs,
	const std::vector<std::vector<cv::Mat>> tvecs);//输出联合标定的内外参结果
	void Calibration_errorRating(std::vector<std::vector<cv::Point3f>> objectPoints,
		std::vector<std::vector<cv::Point2f>>cornersPointBuf,cv::Mat cameraMatrix,cv::Mat distCoeffs,std::vector<cv::Mat> rvecs , std::vector<cv::Mat> tvecs);
	std::vector<cv::Mat>& GetcameraData();//获取标定图片时启用
	void Icp();
	~CCalibration(){};
	
private:
	void Cornerspoint_2D3DCoordinateAcquisition(const std::vector<std::vector<cv::Point2f>> all_corners2DPointBuf,std::vector<std::vector<cv::Point2f>> &one_corners2DPointBuf,
		std::vector<std::vector<cv::Point2f>>& two_corners2DPointBuf,std::vector<std::vector<cv::Point2f>>&three_corners2DPointBuf ,std::vector<std::vector<cv::Point2f>>&four_corners2DPointBuf ,
		std::vector<std::vector<cv::Point3f>>& one_object3DPoints,std::vector<std::vector<cv::Point3f>>& two_object3DPoints,std::vector<std::vector<cv::Point3f>>& three_object3DPoints,std::vector<std::vector<cv::Point3f>> &four_object3DPoints);
	void Cornerspoint_2D3DCoordinateAcquisition(const std::vector<std::vector<cv::Point2f>> all_corners2DPointBuf, std::vector<std::vector<cv::Point2f>> &one_corners2DPointBuf,
		std::vector<std::vector<cv::Point2f>>& two_corners2DPointBuf,std::vector<std::vector<cv::Point3f>>& one_object3DPoints, std::vector<std::vector<cv::Point3f>>& two_object3DPoints);
	std::string SrcimgPath;
    cv::Size broadSize;cv::Size firstImgSize;cv::Size realBroadSize;
	int Img_count,PicturegroupNum;
	std::vector<cv::Mat> AllpointCloud;
	
};

class CUndistort
{
public:
	CUndistort();
	CUndistort(const std::string srcImgPath , const std::string IntrinsicsPath)
	{
		this->srcImgPath = srcImgPath;
		this->IntrinsicsPath = IntrinsicsPath;
	};
	~CUndistort(){};
public:
	void initUndistortRectifyMap();
	void getMatrixParameters();
private:
	std::string srcImgPath;
	std::string IntrinsicsPath;
	cv::Mat cameraMatrix;
	cv::Mat distCoeffs;
};





#endif //Calibration.h
