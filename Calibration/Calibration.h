#ifndef  _CALIBRATION_H_
#define  _CALIBRATION_H_

#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
//�����豸����
#define  camereNum 2
////�����豸IRͼ����R��ȡ��Ƭ����֮�������궨���� 
//#define IR
//#define calibration 
//�����豸���»�ȡ����yml���ݣ����data���Ѿ������ݣ��򲻱��������ɲ����ֿ��عر�
//#define ifNeedgetPointData 
class CCalibration
{
public:
	CCalibration(){};
	CCalibration(std::string SrcimgPath ,const cv::Size broadSize,const cv::Size realBroadSize)  //��ʼ����
	{  this->SrcimgPath = SrcimgPath; this->broadSize = broadSize;this->realBroadSize = realBroadSize;}
	void  findChessboardCorners(std::vector<cv::Point2f> &cornersPointBuf,
		std::vector<std::vector<cv::Point2f>> &all_cornersPointBuf) ;
	void SingleCalibration();
	void CalUnionCalibrationParameters(std::vector<std::vector<cv::Mat>> &rvecs,
	const std::vector<std::vector<cv::Mat>> tvecs);//������ϱ궨������ν��
	void Calibration_errorRating(std::vector<std::vector<cv::Point3f>> objectPoints,
		std::vector<std::vector<cv::Point2f>>cornersPointBuf,cv::Mat cameraMatrix,cv::Mat distCoeffs,std::vector<cv::Mat> rvecs , std::vector<cv::Mat> tvecs);
	std::vector<cv::Mat>& GetcameraData();//��ȡ�궨ͼƬʱ����
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
