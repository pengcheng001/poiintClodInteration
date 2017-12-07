#include "Calibration.h"
#include "../include/ICP/MIcpPointToPlane.h"
#include <fstream>
#include <vector>

bool flag[20] = {false};

void CCalibration::SingleCalibration()
{



	std::vector<cv::Point2f> cornersPointBuf;
	std::vector<std::vector<cv::Point2f>> all_corners2DPointBuf;
	CCalibration::findChessboardCorners(cornersPointBuf,all_corners2DPointBuf);
	//=========================将角点转为三维坐标存到all_cornersPointBuf中


	//根据真实长度的棋盘，来生成建立在棋盘上的世界坐标中角点的三维坐标
	Img_count = all_corners2DPointBuf.size();
	std::vector<cv::Mat> cameraMatrix(camereNum);
	std::vector<cv::Mat> distCoeffs(camereNum);
	//cv::Mat cameraMatrix = cv::Mat(3,3,CV_32FC1,cv::Scalar(0,0,0));
	//cv::Mat distCoeffs = cv::Mat(1,5,CV_32FC1,cv::Scalar(0,0,0));

	//四个相机的外参
	std::vector<std::vector<cv::Mat>> rvecs(camereNum);
	std::vector<std::vector<cv::Mat>> tvecs(camereNum);



	//相机标定程序
	std::cout << "开始相机标定程序："<<std::endl;
	if (camereNum == 2)
	{
		std::vector<std::vector<cv::Point3f>> one_object3DPoints;
		std::vector<std::vector<cv::Point3f>> two_object3DPoints;
		std::vector<std::vector<cv::Point2f>> one_corners2DPointBuf;
		std::vector<std::vector<cv::Point2f>> two_corners2DPointBuf;
		Cornerspoint_2D3DCoordinateAcquisition(all_corners2DPointBuf, one_corners2DPointBuf, two_corners2DPointBuf,
			one_object3DPoints, two_object3DPoints);
		cv::calibrateCamera(one_object3DPoints, one_corners2DPointBuf, firstImgSize, cameraMatrix[0], distCoeffs[0], rvecs[0], tvecs[0], 0);
		cv::calibrateCamera(two_object3DPoints, two_corners2DPointBuf, firstImgSize, cameraMatrix[1], distCoeffs[1], rvecs[1], tvecs[1], 0);
		std::cout << "一号模组标定误差评定:" << std::endl;
		Calibration_errorRating(one_object3DPoints, one_corners2DPointBuf, cameraMatrix[0], distCoeffs[0], rvecs[0], tvecs[0]);
		std::cout << std::endl;
		std::cout << "二号模组标定误差评定:" << std::endl;
		Calibration_errorRating(two_object3DPoints, two_corners2DPointBuf, cameraMatrix[1], distCoeffs[1], rvecs[1], tvecs[1]);
		std::cout << std::endl;
		std::cout << "评价完成，正在保存标定参数，请稍后......" << std::endl;
		//保存相机联合标定外参，即相机之间的RT关系
		CalUnionCalibrationParameters(rvecs, tvecs);
	}
	else
	{
		std::vector<std::vector<cv::Point3f>> one_object3DPoints;
		std::vector<std::vector<cv::Point3f>> two_object3DPoints;
		std::vector<std::vector<cv::Point3f>> three_object3DPoints;
		std::vector<std::vector<cv::Point3f>> four_object3DPoints;
		std::vector<std::vector<cv::Point2f>> one_corners2DPointBuf;
		std::vector<std::vector<cv::Point2f>> two_corners2DPointBuf;
		std::vector<std::vector<cv::Point2f>> three_corners2DPointBuf;
		std::vector<std::vector<cv::Point2f>> four_corners2DPointBuf;
		Cornerspoint_2D3DCoordinateAcquisition(all_corners2DPointBuf, one_corners2DPointBuf, two_corners2DPointBuf, three_corners2DPointBuf, four_corners2DPointBuf,
			one_object3DPoints, two_object3DPoints, three_object3DPoints, four_object3DPoints);
		cv::calibrateCamera(one_object3DPoints, one_corners2DPointBuf, firstImgSize, cameraMatrix[0], distCoeffs[0], rvecs[0], tvecs[0], 0);
		cv::calibrateCamera(two_object3DPoints, two_corners2DPointBuf, firstImgSize, cameraMatrix[1], distCoeffs[1], rvecs[1], tvecs[1], 0);
		cv::calibrateCamera(three_object3DPoints, three_corners2DPointBuf, firstImgSize, cameraMatrix[2], distCoeffs[2], rvecs[2], tvecs[2], 0);
		cv::calibrateCamera(four_object3DPoints, four_corners2DPointBuf, firstImgSize, cameraMatrix[3], distCoeffs[3], rvecs[3], tvecs[3], 0);
		//各相机角点误差评定，推荐Totalerr小于0.15时使用
		std::cout << "一号模组标定误差评定:" << std::endl;
		Calibration_errorRating(one_object3DPoints, one_corners2DPointBuf, cameraMatrix[0], distCoeffs[0], rvecs[0], tvecs[0]);
		std::cout << std::endl;
		std::cout << "二号模组标定误差评定:" << std::endl;
		Calibration_errorRating(two_object3DPoints, two_corners2DPointBuf, cameraMatrix[1], distCoeffs[1], rvecs[1], tvecs[1]);
		std::cout << std::endl;
		std::cout << "三号模组标定误差评定:" << std::endl;
		Calibration_errorRating(three_object3DPoints, three_corners2DPointBuf, cameraMatrix[2], distCoeffs[2], rvecs[2], tvecs[2]);
		std::cout << std::endl;
		std::cout << "四号模组标定误差评定:" << std::endl;
		Calibration_errorRating(four_object3DPoints, four_corners2DPointBuf, cameraMatrix[3], distCoeffs[3], rvecs[3], tvecs[3]);
		std::cout << std::endl;
		std::cout << "评价完成，正在保存标定参数，请稍后......" << std::endl;
		//保存相机联合标定外参，即相机之间的RT关系
		CalUnionCalibrationParameters(rvecs, tvecs);
	}



	//保存四个ir相机的平均内参
	cv::Mat AvecameraMatrix(cv::Size(3,3),CV_64FC1,cv::Scalar(0));
	for (int i=0;i<camereNum;i++)
	{
		AvecameraMatrix += cameraMatrix[i];
	}
	AvecameraMatrix /= camereNum;
	cv::FileStorage fs("data/intrinsics.yml",cv::FileStorage::WRITE);
	//fs << "cameraMatrix" <<cameraMatrix << "distCoeffs"<<distCoeffs ;
	fs << "cameraMatrix" <<AvecameraMatrix ;
	fs.release();
}

void CCalibration::Icp()
{
#ifdef ifNeedgetPointData  //需要获取时定义该开关，按c获取点云数据

	GetcameraData();
	float** Mod = (float**)malloc(4 * sizeof(float*));
	int cloudNum[camereNum] = { 0 };
	char txt[10];
	for (int i = 0; i < camereNum; i++)
	{
		std::sprintf(txt, "out%d.txt", i);
		std::ofstream out(txt);
		Mod[i] = (float*)calloc(3 * 640 * 480, sizeof(float));
		int n = 0, cloudnum = 0;
		for (int j = 0; j < AllpointCloud[i].rows; j++)
		{
			for (int h = 0; h < AllpointCloud[i].cols; h++)
			{
				if (_isnan(AllpointCloud[i].at<cv::Vec3f>(j, h)[0]))
					continue;
				Mod[i][n * 3 + 0] = AllpointCloud[i].at<cv::Vec3f>(j, h)[0];
				Mod[i][n * 3 + 1] = AllpointCloud[i].at<cv::Vec3f>(j, h)[1];
				Mod[i][n * 3 + 2] = AllpointCloud[i].at<cv::Vec3f>(j, h)[2];
				out << Mod[i][n * 3 + 0] << " " << Mod[i][n * 3 + 1] << " " << Mod[i][n * 3 + 2] << std::endl;
				n++; cloudnum++;
			}
		}
		out.close();
		cloudNum[i] = cloudnum;
	}
	//读取外参赋值icp初始参数,如果是两个参数，则不用首尾循环
	if (camereNum == 2)
	{
		cv::FileStorage fs2("data/extrinsic.yml", cv::FileStorage::READ);
		char ss[50], ss1[50], ss2[50];
		cv::Mat R;
		cv::Mat	T;
		Matrix Ri(3, 3), Ti(3, 1);


		std::sprintf(ss1, "trance_matrix-%dTo%d", 0, 1);
		std::sprintf(ss, "rotation_matrix-%dTo%d", 0, 1);
		std::sprintf(ss2, "ICP_RT-%dTo%d", 0, 1);
		fs2[ss] >> R; fs2[ss1] >> T; fs2.release();
		for (int h = 0; h < 3; h++)
			for (int g = 0; g < 3; g++)
				Ri.val[h][g] = R.at<double>(h, g);
		Ti.val[0][0] = T.at<double>(0, 0);
		Ti.val[1][0] = T.at<double>(1, 0);
		Ti.val[2][0] = T.at<double>(2, 0);
		int32_t numMod = cloudNum[1];
		int32_t numTep = cloudNum[0];
		MIcpPointToPlane icp(Mod[1], numMod, 3);
		icp.fit(Mod[0], numTep, Ri, Ti, -1);
		std::ofstream out(ss2);
		out << Ri;  out << std::endl;
		out << Ti;  out << std::endl;
		out.close();
	}
	else
	{
		for (int i = 0; i < AllpointCloud.size(); i++)
		{
			cv::FileStorage fs2("data/extrinsic.yml", cv::FileStorage::READ);
			char ss[50], ss1[50], ss2[50];
			cv::Mat R;
			cv::Mat	T;
			Matrix Ri(3, 3), Ti(3, 1);
			if (i == 3)
			{
				std::sprintf(ss1, "trance_matrix-%dTo%d", 0, 3);
				std::sprintf(ss, "rotation_matrix-%dTo%d", 0, 3);
				std::sprintf(ss2, "ICP_RT-%dTo%d", 0, 3);
				fs2[ss] >> R; fs2[ss1] >> T; fs2.release();
				//赋值旋转平移
				for (int h = 0; h < 3; h++)
					for (int g = 0; g < 3; g++)
						Ri.val[h][g] = R.at<double>(h, g);
				Ti.val[0][0] = T.at<double>(0, 0);
				Ti.val[1][0] = T.at<double>(1, 0);
				Ti.val[2][0] = T.at<double>(2, 0);
				//每个点云的数量
				int32_t numMod = cloudNum[3];
				int32_t numTep = cloudNum[0];
				MIcpPointToPlane icp(Mod[3], numMod, 3);
				icp.fit(Mod[0], numTep, Ri, Ti, -1);
				std::ofstream out(ss2);
				out << Ri;  out << std::endl;
				out << Ti;  out << std::endl;
				out.close();
			}
			else
			{
				std::sprintf(ss1, "trance_matrix-%dTo%d", i, i + 1);
				std::sprintf(ss, "rotation_matrix-%dTo%d", i, i + 1);
				std::sprintf(ss2, "ICP_RT-%dTo%d", 0, 3);
				fs2[ss] >> R; fs2[ss1] >> T; fs2.release();
				for (int h = 0; h < 3; h++)
					for (int g = 0; g < 3; g++)
						Ri.val[h][g] = R.at<double>(h, g);
				Ti.val[0][0] = T.at<double>(0, 0);
				Ti.val[1][0] = T.at<double>(1, 0);
				Ti.val[2][0] = T.at<double>(2, 0);
				int32_t numMod = cloudNum[i + 1];
				int32_t numTep = cloudNum[i];
				MIcpPointToPlane icp(Mod[i + 1], numMod, 3);
				icp.fit(Mod[i], numTep, Ri, Ti, -1);
				std::ofstream out(ss2);
				out << Ri;  out << std::endl;
				out << Ti;  out << std::endl;
				out.close();
			}
		}
	}

	for (int i = 0; i < 4; i++)
		free(Mod[i]);
	free(Mod);
#else
	//避开从相机获取数据，直接从之前保存的yml文件读取
	float** Mod = (float**)malloc(4 * sizeof(float*));
	char win[50];
	char ss[50];
	std::vector<cv::Mat> DataVMat(camereNum);
	for (int deviceNum = 0; deviceNum < camereNum; deviceNum++)
	{
		std::sprintf(win, "data/pointcloud(%d).yml", deviceNum);
		std::sprintf(ss, "cloudData");
		cv::FileStorage fs(win, cv::FileStorage::READ);
		fs[ss] >> DataVMat[deviceNum];
		fs.release();
	}
	int cloudNum[camereNum] = { 0 };
	for (int i = 0; i < camereNum; i++)
	{
		Mod[i] = (float*)calloc(3 * 640 * 480, sizeof(float));
		int n = 0, cloudnum = 0;
		for (int j = 0; j < DataVMat[i].rows; j++)
		{
			for (int h = 0; h < DataVMat[i].cols; h++)
			{
				if (_isnan(DataVMat[i].at<cv::Vec3f>(j, h)[0]))
					continue;
				Mod[i][n * 3 + 0] = DataVMat[i].at<cv::Vec3f>(j, h)[0];
				Mod[i][n * 3 + 1] = DataVMat[i].at<cv::Vec3f>(j, h)[1];
				Mod[i][n * 3 + 2] = DataVMat[i].at<cv::Vec3f>(j, h)[2];
				n++; cloudnum++;
			}
		}
		cloudNum[i] = cloudnum;
	}
	if (camereNum == 2)
	{
		cv::FileStorage fs2("data/extrinsic.yml", cv::FileStorage::READ);
		char ss[50], ss1[50], ss2[50];
		cv::Mat R;
		cv::Mat T;
		Matrix Ri(3, 3), Ti(3, 1);

		std::sprintf(ss1, "trance_matrix-%dTo%d", 0, 1);
		std::sprintf(ss, "rotation_matrix-%dTo%d", 0, 1);
		std::sprintf(ss2, "ICP_RT-%dTo%d", 0, 1);
		fs2[ss] >> R; fs2[ss1] >> T;
		fs2.release();
		for (int h = 0; h < 3; h++)
			for (int g = 0; g < 3; g++)
				Ri.val[h][g] = R.at<double>(h, g);
		Ti.val[0][0] = T.at<double>(0, 0);
		Ti.val[1][0] = T.at<double>(1, 0);
		Ti.val[2][0] = T.at<double>(2, 0);
		int32_t numMod = cloudNum[1];
		int32_t numTep = cloudNum[0];
		MIcpPointToPlane icp(Mod[1], numMod, 3);
		icp.fit(Mod[0], numTep, Ri, Ti, -1);
		std::ofstream out(ss2);
		out << Ri;  out << std::endl;
		out << Ti;  out << std::endl;
		out.close();
	}
	else
	{
		for (int i = 0; i < DataVMat.size(); i++)
		{
			cv::FileStorage fs2("data/extrinsic.yml", cv::FileStorage::READ);
			char ss[50], ss1[50], ss2[50];
			cv::Mat R;
			cv::Mat T;
			Matrix Ri(3, 3), Ti(3, 1);
			if (i == 3)
			{
				std::sprintf(ss1, "trance_matrix-%dTo%d", 0, 3);
				std::sprintf(ss, "rotation_matrix-%dTo%d", 0, 3);
				std::sprintf(ss2, "ICP_RT-%dTo%d", 0, 3);
				fs2[ss] >> R; fs2[ss1] >> T;
				fs2.release();
				for (int h = 0; h < 3; h++)
					for (int g = 0; g < 3; g++)
						Ri.val[h][g] = R.at<double>(h, g);
				Ti.val[0][0] = T.at<double>(0, 0);
				Ti.val[1][0] = T.at<double>(1, 0);
				Ti.val[2][0] = T.at<double>(2, 0);
				int32_t numMod = cloudNum[3];
				int32_t numTep = cloudNum[0];
				MIcpPointToPlane icp(Mod[3], numMod, 3);
				icp.fit(Mod[0], numTep, Ri, Ti, -1);
				std::ofstream out(ss2);
				out << Ri;  out << std::endl;
				out << Ti;  out << std::endl;
				out.close();
			}
			else
			{
				std::sprintf(ss1, "trance_matrix-%dTo%d", i, i + 1);
				std::sprintf(ss, "rotation_matrix-%dTo%d", i, i + 1);
				std::sprintf(ss2, "ICP_RT-%dTo%d", i, i + 1);
				fs2[ss] >> R; fs2[ss1] >> T;
				fs2.release();
				for (int h = 0; h < 3; h++)
					for (int g = 0; g < 3; g++)
						Ri.val[h][g] = R.at<double>(h, g);
				Ti.val[0][0] = T.at<double>(0, 0);
				Ti.val[1][0] = T.at<double>(1, 0);
				Ti.val[2][0] = T.at<double>(2, 0);
				int32_t numMod = cloudNum[i + 1];
				int32_t numTep = cloudNum[i];
				MIcpPointToPlane icp(Mod[i + 1], numMod, 3);
				icp.fit(Mod[i], numTep, Ri, Ti, -1);
				std::ofstream out(ss2);
				out << Ri;  out << std::endl;
				out << Ti;  out << std::endl;
				out.close();
			}
		}
	}

	for (int i = 0; i < 4; i++)
		free(Mod[i]);
	free(Mod);
#endif

	//	IcpPointToPlane icp1;
	//	IcpPointToPoint icp2;
}

//void Calibration_errorRating(std::vector<std::vector<cv::Point3f>> all_object3DPoints,
	//std::vector<std::vector<cv::Point2f>> all_corners2DPointBuf,cv::Mat cameraMatrix,cv::Mat distCoeffs,std::vector<cv::Mat> rvecs , std::vector<cv::Mat> tvecs);
void  CCalibration::findChessboardCorners(std::vector<cv::Point2f> &cornersPointBuf,
										 std::vector<std::vector<cv::Point2f>> &all_cornersPointBuf) 
{
	
	std::ifstream fin(SrcimgPath); // 静态成员函数调用私有成员函数的方法
	std::cout << "开始提取角点，请稍等"<<std::endl;
	std::string filename;
	int image_count=0;
	while (getline(fin , filename))
	{
		cv::Mat Srcimg = cv::imread(filename);
		cv::cvtColor(Srcimg,Srcimg,CV_BGR2GRAY);
		image_count++;
		if (image_count==1)
		{
			firstImgSize.width = Srcimg.cols;
			firstImgSize.height = Srcimg.rows;
			std::cout << "imgSize.width="<<firstImgSize.width<<std::endl;
			std::cout << "imgSize.height="<<firstImgSize.height<<std::endl;
		}
		if (!cv::findChessboardCorners(Srcimg , broadSize , cornersPointBuf))
		{
			std::cout <<filename<< " : 找不到角点信息---------->"<<"自动过滤成功"<<std::endl;
			all_cornersPointBuf.push_back(cornersPointBuf);
			continue;
			//exit(1);
		}
		else
		{
			flag[image_count-1] = true;
			cv::cornerSubPix(Srcimg , cornersPointBuf , cv::Size(5,5),cv::Size(-1,-1),
					cv::TermCriteria(CV_TERMCRIT_EPS+CV_TERMCRIT_ITER,30,0.1));
			all_cornersPointBuf.push_back(cornersPointBuf);//存取每幅图片的角点
			cv::drawChessboardCorners(Srcimg,CCalibration::broadSize,cornersPointBuf,false);//绘制绘制每幅图的角点
			cv::imshow("drawChessboardCorners",Srcimg);
			cv::waitKey(50);
		}
	}
	std::cout << "每张图片的角点数量："<<CCalibration::broadSize.width * CCalibration::broadSize.height <<std::endl;
	std::cout << "角点提取成功！"<<std::endl;
}

void CCalibration::Calibration_errorRating(std::vector<std::vector<cv::Point3f>>objectPoints,
	std::vector<std::vector<cv::Point2f>>cornersPointBuf,cv::Mat cameraMatrix,cv::Mat distCoeffs,std::vector<cv::Mat> rvecs , std::vector<cv::Mat> tvecs)
{
	
	double  err=0.0;
	double total_err = 0.0;
	for (int i=0;i<objectPoints.size();i++)
	{
		std::vector<cv::Point2f> tempPoints;
		cv::projectPoints(objectPoints[i],rvecs[i],tvecs[i],cameraMatrix,distCoeffs,tempPoints);
		err = cv::norm(tempPoints,cornersPointBuf[i],cv::NORM_L2);
		err/= tempPoints.size();
		total_err+=err;
	}
	std::cout <<"该模组角点的累计误差为："<< total_err <<std::endl;
	std::cout <<"该模组角点的平均误差为："<< total_err/objectPoints.size() <<std::endl;
	//std::cout << "评价完成，正在保存标定参数，请稍后"<<std::endl;


	//cv::FileStorage fs("intrinsics.yml",cv::FileStorage::WRITE);
	//fs << "cameraMatrix" <<cameraMatrix << "distCoeffs"<<distCoeffs ;
	//fs.release();


	
}

//换算真实三维坐标
void CCalibration::Cornerspoint_2D3DCoordinateAcquisition(const std::vector<std::vector<cv::Point2f>> all_corners2DPointBuf,std::vector<std::vector<cv::Point2f>> &one_corners2DPointBuf,std::vector<std::vector<cv::Point2f>>& two_corners2DPointBuf, std::vector<std::vector<cv::Point2f>>&three_corners2DPointBuf ,std::vector<std::vector<cv::Point2f>>&four_corners2DPointBuf , std::vector<std::vector<cv::Point3f>>& one_object3DPoints,std::vector<std::vector<cv::Point3f>>& two_object3DPoints, std::vector<std::vector<cv::Point3f>>& three_object3DPoints,std::vector<std::vector<cv::Point3f>> &four_object3DPoints)
{
	std::vector<cv::Point3f> temp;

	for (int i=0;i<CCalibration::broadSize.height;i++)
	{
		for (int j=0;j<CCalibration::broadSize.width;j++)
		{
			cv::Point3f  RealPoint;
			RealPoint.x = i * realBroadSize.width;
			RealPoint.y = j * realBroadSize.height;
			RealPoint.z = 0;
			temp.push_back(RealPoint);
		}
	}


	for (int i = 0; i < Img_count; i = i + camereNum)
	{
		//如果该组中有一幅图片未检测到角点，忽略本组图片。
		if (flag[i] == false||flag[i+1] == false||flag[i+2] == false||flag[i+3] == false)
			continue;
		else
		{
			one_corners2DPointBuf.push_back(all_corners2DPointBuf[i]);
			two_corners2DPointBuf.push_back(all_corners2DPointBuf[i+1]);
			three_corners2DPointBuf.push_back(all_corners2DPointBuf[i+2]);
			four_corners2DPointBuf.push_back(all_corners2DPointBuf[i+3]);

			one_object3DPoints.push_back(temp);
			two_object3DPoints.push_back(temp);
			three_object3DPoints.push_back(temp);
			four_object3DPoints.push_back(temp);
		}
	}
	PicturegroupNum = one_object3DPoints.size();
}

void CCalibration::Cornerspoint_2D3DCoordinateAcquisition(const std::vector<std::vector<cv::Point2f>> all_corners2DPointBuf, std::vector<std::vector<cv::Point2f>> &one_corners2DPointBuf, std::vector<std::vector<cv::Point2f>>& two_corners2DPointBuf,
	std::vector<std::vector<cv::Point3f>>& one_object3DPoints, std::vector<std::vector<cv::Point3f>>& two_object3DPoints)
{
	std::vector<cv::Point3f> temp;

	for (int i = 0; i < CCalibration::broadSize.height; i++)
	{
		for (int j = 0; j < CCalibration::broadSize.width; j++)
		{
			cv::Point3d  RealPoint;
			RealPoint.x = i * realBroadSize.width;
			RealPoint.y = j * realBroadSize.height;
			RealPoint.z = 0;
			temp.push_back(RealPoint);
		}
	}


	for (int i = 0; i < Img_count; i = i + camereNum)
	{
		//如果该组中有一幅图片未检测到角点，忽略本组图片。
		if (flag[i] == false || flag[i + 1] == false )
			continue;
		else
		{
			one_corners2DPointBuf.push_back(all_corners2DPointBuf[i]);
			two_corners2DPointBuf.push_back(all_corners2DPointBuf[i + 1]);
			one_object3DPoints.push_back(temp);
			two_object3DPoints.push_back(temp);
		}
	}
	PicturegroupNum = one_object3DPoints.size();
}

void CCalibration::CalUnionCalibrationParameters(std::vector<std::vector<cv::Mat>> &rvecs, const std::vector<std::vector<cv::Mat>> tvecs)
{
	std::vector<cv::Mat> CrossRvec(camereNum);//1->4 , 2->3 , 3->4 
	std::vector<cv::Mat> CrossTvec(camereNum);
	cv::Mat SumtempR(cv::Size(3,3),CV_64FC1,cv::Scalar(0));
	cv::Mat SumtempT(cv::Size(1,3),CV_64FC1,cv::Scalar(0));
	for (int i = 0; i<camereNum; i++)
		for(int j=0;j<PicturegroupNum;j++)
			cv::Rodrigues(rvecs[i][j],rvecs[i][j]);

	//Img_count代表图片组数
	if (camereNum == 2)
	{
		for (int j = 0; j < PicturegroupNum; j++)
		{
			cv::Mat tempR(cv::Size(3, 3), CV_64FC1, cv::Scalar(0));
			cv::Mat tempT(cv::Size(1, 3), CV_64FC1, cv::Scalar(0));

			tempR = rvecs[1][j] * rvecs[0][j].t();
			SumtempR += tempR;
			tempT = tvecs[1][j] - tempR * tvecs[0][j];
			SumtempT += tempT;
		}
		CrossRvec[0] = SumtempR / PicturegroupNum;
		CrossTvec[0] = SumtempT / PicturegroupNum;
		cv::FileStorage fs2("data/extrinsic.yml", cv::FileStorage::WRITE);
		char ss[50];
		char ss1[50];
		std::sprintf(ss1, "trance_matrix-%dTo%d", 0, 1);
		std::sprintf(ss, "rotation_matrix-%dTo%d", 0, 1);
		fs2 << ss << CrossRvec[0] << ss1 << CrossTvec[0];
		fs2.release();
	}
	else
	{
		for (int i = 0; i < camereNum; i++)
		{
			/*	for (int j=0;j<PicturegroupNum;j++)
			{
			if (i==3)
			{
			tempR = rvecs[0][j] * rvecs[i][j].t();
			SumtempR += tempR;
			tempT = tvecs[0][j] - tempR * tvecs[i][j];
			SumtempT += tempT;
			std::cout << tempR<<std::endl;
			std::cout << tempT<<std::endl;
			}
			else
			{
			tempR = rvecs[i][j] * rvecs[i+1][j].t();

			SumtempR += tempR;
			tempT = tvecs[i][j] - tempR * tvecs[i+1][j];
			SumtempT += tempT;
			std::cout << tempR<<std::endl;
			std::cout << tempT<<std::endl;
			}

			}*/
			for (int j = 0; j < PicturegroupNum; j++)
			{
				cv::Mat tempR(cv::Size(3, 3), CV_64FC1, cv::Scalar(0));
				cv::Mat tempT(cv::Size(1, 3), CV_64FC1, cv::Scalar(0));
				if (i == 3)
				{
					tempR = rvecs[3][j] * rvecs[0][j].t();
					SumtempR += tempR;
					tempT = tvecs[3][j] - tempR * tvecs[0][j];
					SumtempT += tempT;
				}
				else
				{
					tempR = rvecs[i + 1][j] * rvecs[i][j].t();

					SumtempR += tempR;
					tempT = tvecs[i + 1][j] - tempR * tvecs[i][j];
					SumtempT += tempT;
				}

			}
			CrossRvec[i] = SumtempR / PicturegroupNum;
			CrossTvec[i] = SumtempT / PicturegroupNum;
		}
		cv::FileStorage fs2("data/extrinsic.yml", cv::FileStorage::WRITE);
		char ss[50];
		char ss1[50];
		for (int i = 0; i < camereNum; i++)
		{
			if (i == 3)
			{
				std::sprintf(ss1, "trance_matrix-%dTo%d", 0, 3);
				std::sprintf(ss, "rotation_matrix-%dTo%d", 0, 3);
				fs2 << ss << CrossRvec[i] << ss1 << CrossTvec[i];
			}
			else
			{
				std::sprintf(ss1, "trance_matrix-%dTo%d", i, i + 1);
				std::sprintf(ss, "rotation_matrix-%dTo%d", i, i + 1);
				fs2 << ss << CrossRvec[i] << ss1 << CrossTvec[i];
			}
		}
		fs2.release();
	}

	
	std::cout << "保存完成，请到目录查看"<<std::endl;
}

