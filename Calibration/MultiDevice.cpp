#include "../include/ty/common.hpp"
#include "Calibration.h"
#include <fstream>
static char buffer[1024*1024];


struct CamInfo
{
    char                sn[32];
    TY_DEV_HANDLE       hDev;
    char*               fb[2];
    TY_FRAME_DATA       frame;
    int                 idx;
    DepthRender         render;

    CamInfo() : hDev(0), idx(0) {fb[0]=0; fb[1]=0;}
};


/*void frameHandler(TY_FRAME_DATA* frame, void* userdata)
{
    CamInfo* pData = (CamInfo*) userdata;

    cv::Mat depth, irl, irr, color;
    parseFrame(*frame, &depth, &irl, &irr, &color, 0);

    char win[64];
    if(!depth.empty()){
        cv::Mat colorDepth = pData->render.Compute(depth);
        sprintf(win, "depth-%s", pData->sn);
		std::cout << pData->sn;
        cv::imshow(win, colorDepth);
    }
    if(!irl.empty()){
        sprintf(win, "LeftIR-%s", pData->sn);
        cv::imshow(win, irl);
    }
    if(!irr.empty()){
        sprintf(win, "RightIR-%s", pData->sn);
        cv::imshow(win, irr);
    }
    if(!color.empty()){
        sprintf(win, "color-%s", pData->sn);
        cv::imshow(win, color);
    }

    pData->idx++;
    //LOGD("=== Callback: Re-enqueue buffer(%p, %d)", frame->userBuffer, frame->bufferSize);
    ASSERT_OK( TYEnqueueBuffer(pData->hDev, frame->userBuffer, frame->bufferSize) );
}*/

std::vector<cv::Mat>& CCalibration::GetcameraData()
{
    LOGD("=== Init lib");
    ASSERT_OK( TYInitLib() );
    TY_VERSION_INFO* pVer = (TY_VERSION_INFO*)buffer;
    ASSERT_OK( TYLibVersion(pVer) );
    LOGD("     - lib version: %d.%d.%d", pVer->major, pVer->minor, pVer->patch);

    LOGD("=== Get device info");
    int n = camereNum ;
    LOGD("     - device number %d", n);

    TY_DEVICE_BASE_INFO* pBaseInfo = (TY_DEVICE_BASE_INFO*)buffer;
    ASSERT_OK( TYGetDeviceList(pBaseInfo, 100, &n) );

    if(n < 1){
        LOGD("=== Need more than 1 devices");
    }

    std::vector<CamInfo> cams(n);
    for(int i = 0; i < n; i++){
        LOGD("=== Open device %d (id: %s)", i, pBaseInfo[i].id);
        strncpy(cams[i].sn, pBaseInfo[i].id, sizeof(cams[i].sn));

        ASSERT_OK( TYOpenDevice(pBaseInfo[i].id, &cams[i].hDev) );

        int32_t allComps;
        ASSERT_OK( TYGetComponentIDs(cams[i].hDev, &allComps) );
        if(0 && allComps & TY_COMPONENT_RGB_CAM){
            LOGD("=== Has RGB camera, open RGB cam");
            ASSERT_OK( TYEnableComponents(cams[i].hDev, TY_COMPONENT_RGB_CAM) );
        }

        LOGD("=== Configure components, open depth cam");
#ifdef IR
		int32_t componentIDs =  TY_COMPONENT_IR_CAM_LEFT;
#else
		int32_t componentIDs = TY_COMPONENT_POINT3D_CAM | TY_COMPONENT_RGB_CAM;
#endif // IR

        
    	ASSERT_OK( TYEnableComponents(cams[i].hDev, componentIDs) );
        LOGD("=== Configure feature, set resolution to 640x480.");
        int err = TYSetEnum(cams[i].hDev, TY_COMPONENT_DEPTH_CAM, TY_ENUM_IMAGE_MODE, TY_IMAGE_MODE_640x480);
		ASSERT(err == TY_STATUS_OK || err == TY_STATUS_NOT_PERMITTED);

        LOGD("=== Prepare image buffer");
        int32_t frameSize;
        ASSERT_OK( TYGetFrameBufferSize(cams[i].hDev, &frameSize) );
        LOGD("     - Get size of framebuffer, %d", frameSize);
        ASSERT( frameSize >= 640*480*2 );

        LOGD("     - Allocate & enqueue buffers");
        cams[i].fb[0] = new char[frameSize];
        cams[i].fb[1] = new char[frameSize];
        LOGD("     - Enqueue buffer (%p, %d)", cams[i].fb[0], frameSize);
        ASSERT_OK( TYEnqueueBuffer(cams[i].hDev, cams[i].fb[0], frameSize) );
        LOGD("     - Enqueue buffer (%p, %d)", cams[i].fb[1], frameSize);
        ASSERT_OK( TYEnqueueBuffer(cams[i].hDev, cams[i].fb[1], frameSize) );

#ifdef IR
		bool triggerMode = false;
#else
		bool triggerMode = false;
#endif //
        
        LOGD("=== Set trigger mode %d", triggerMode);
        ASSERT_OK( TYSetBool(cams[i].hDev, TY_COMPONENT_DEVICE, TY_BOOL_TRIGGER_MODE, triggerMode) );
        LOGD("=== Start capture");
        ASSERT_OK( TYStartCapture(cams[i].hDev) );
    }

    LOGD("=== While loop to fetch frame");
    bool exit_main = false;
#ifdef IR
	cv::namedWindow("ty",0);
	//存储各个ir的图片，按r触发
	int whileNum = 1;
	std::ofstream outFile("data/test.txt", std::ios_base::out); 
#else
	cv::namedWindow("pointCloud",0);
#endif // IR

	std::cout<<"command:\n\tq:exit\n\tr:cut off a set of pictures\n\tc:save pointCloud data"<<std::endl;
    while(!exit_main){
		
#ifdef IR
		cv::Mat mergeImg(cv::Size(640*camereNum,480),CV_8UC1,cv::Scalar(255));
#else
		cv::Mat mergePointImg(cv::Size(640*camereNum,480),CV_32FC3,cv::Scalar(255));
		std::vector<cv::Mat> tempSwp;
		cv::Mat pointcloutTemp;
		AllpointCloud.swap(tempSwp);
#endif // IR

        for(int i = 0; i < (cams.size()); i++) {
#ifdef IR
			int err = TYFetchFrame(cams[i].hDev, &cams[i].frame, 1000);
#else
			int err = TYFetchFrame(cams[i].hDev, &cams[i].frame, -1);
#endif // IR    
            if( err != TY_STATUS_OK ){
                LOGD("cam %s %d ... Drop one frame", cams[i].sn, cams[i].idx);
                continue;
            }
			//拼接图片
#ifdef IR
			cv::Mat  irl;
			TY_FRAME_DATA *frame = &cams[i].frame;
			parseFrame(*frame, 0,&irl,0,0,0);
			cams[i].idx++;
			ASSERT_OK( TYEnqueueBuffer(cams[i].hDev, cams[i].frame.userBuffer, cams[i].frame.bufferSize) );
			cv::resize(irl,irl,cv::Size(640,480));
			irl.copyTo(mergeImg(cv::Rect(i*640,0,640,480)));
#else//拼接点云图片，并存储到容器
	cv::Mat  pointCloud;
	TY_FRAME_DATA *frame = &cams[i].frame;
	parseFrame(*frame, 0,0,0,0,&pointCloud);// pointCloud : CV_32FC3
	pointcloutTemp = cv::Mat(pointCloud.size(),CV_32FC3);
	pointcloutTemp = pointCloud.clone();
	cams[i].idx++;
	ASSERT_OK( TYEnqueueBuffer(cams[i].hDev, cams[i].frame.userBuffer, cams[i].frame.bufferSize));
	AllpointCloud.push_back(pointcloutTemp);
	pointcloutTemp.copyTo(mergePointImg(cv::Rect(i*640,0,640,480)));

#endif
   }
		
#ifdef IR
		cv::imshow("ty",mergeImg);
#else
		cv::imshow("pointCloud",mergePointImg);
#endif // IR

        int key = cv::waitKey(1);
        switch(key & 0xff){
            case 0xff:
                break;
            case 'q':
                exit_main = true;
                break;
			case 'c':
				char win[50];
				char ss[50];
				for (int deviceNum = 0; deviceNum<camereNum; deviceNum++)
				{
					std::sprintf(win, "data/pointcloud(%d).yml",deviceNum);
					std::sprintf(ss,"cloudData");
					cv::FileStorage fs(win,cv::FileStorage::WRITE);
					fs << ss <<AllpointCloud[deviceNum];
					fs.release();
				}
				std::cout << "======================================" << std::endl;
				std::cout << "     保存点云数据成功,请到data目录查看    " << std::endl;
				std::cout << "======================================" << std::endl;
				exit_main = true;
				break;
			case 'r':
#ifdef IR
				char winJPG[50];
				char winTXT[50];
				for (int photoNum = 0; photoNum<camereNum; photoNum++)
				{
					sprintf(winJPG, "data/LeftIR-%d(%d).jpg",whileNum,photoNum);
					sprintf(winTXT, "data/LeftIR-%d(%d).jpg",whileNum,photoNum);
					outFile << winTXT;
					outFile << std::endl;
					cv::Mat photoMat = mergeImg(cv::Rect(photoNum*640,0,640,480));
					cv::imwrite(winJPG,photoMat);
					std::cout << "save ok:" << winJPG <<std::endl;
				}
				whileNum++;
#endif //
				break;

        }
    }

    for(int i = 0; i < (cams.size()); i++){
        ASSERT_OK( TYStopCapture(cams[i].hDev) );
        ASSERT_OK( TYCloseDevice(cams[i].hDev) );
        // MSLEEP(10); // sleep to ensure buffer is not used any more
        delete cams[i].fb[0];
        delete cams[i].fb[1];
    }
    ASSERT_OK( TYDeinitLib() );

    LOGD("=== Main done!");
#ifdef IR
		outFile.close();
#endif // IR
		return AllpointCloud;
}
