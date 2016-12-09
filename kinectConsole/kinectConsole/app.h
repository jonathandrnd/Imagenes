#include <Windows.h>
#include <Kinect.h>
#include <opencv2/opencv.hpp>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <vector>
#include <wrl/client.h>

using namespace Microsoft::WRL;
#define SCRWIDTH 512
#define SCRHEIGHT 424
// 1920 1080 para ver sin bordes   512 424 sin color

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef char int8;
typedef short int16;
typedef int int32;
typedef long long int64;

template<typename T>
void SafeRelease(T& ptr){ if (ptr){ ptr->Release(); ptr = nullptr; } }

class App{
	public:
		void Init();
		void Tick(float deltaTime);
		void Shutdown();

		void SetPixelBuffer(uint32* pixelBuffer){ m_pixelBuffer = pixelBuffer; }
		void Plot(int x, int y, uint32 color){
			if (x < 0 || x >= SCRWIDTH || y < 0 || y >= SCRHEIGHT)
				return;
			m_pixelBuffer[x + y * SCRWIDTH] = color;
		}

	private:
		uint32* m_pixelBuffer = nullptr;
		IKinectSensor* m_sensor = nullptr;
		IDepthFrameReader* m_depthFrameReader = nullptr;
		IColorFrameReader* m_colorFrameReader = nullptr;
		ICoordinateMapper* m_coordinateMapper = nullptr;

		uint16* m_depthBuffer = nullptr;
		uint32* m_colorBuffer = nullptr;

		ColorSpacePoint* m_colorSpacePoints = nullptr;
		int m_depthWidth = 0, m_depthHeight = 0;

};