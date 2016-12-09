#include "stdafx.h"
#include "app.h"

void App::Init(){
	HRESULT hr;
	hr = GetDefaultKinectSensor(&m_sensor);
	if (FAILED(hr)){
		printf("Failed to find the kinect sensor");
		exit(10);
	}

	m_sensor->Open(); //hace que se prenda el kinect
	//obtener la profundidad del source
	IDepthFrameSource* depthFrameSource;
	hr=m_sensor->get_DepthFrameSource(&depthFrameSource);

	if (FAILED(hr)){
		printf("Failed to get the depth frame source\n");
		exit(10);
	}

	//get depth frame description
	IFrameDescription *frameDesc;
	depthFrameSource->get_FrameDescription(&frameDesc);
	frameDesc->get_Width(&m_depthWidth);
	frameDesc->get_Height(&m_depthHeight);
	
	//get the depth frame reader
	hr = depthFrameSource->OpenReader(&m_depthFrameReader);
	if (FAILED(hr)){
		printf("Failed to open the depth frame reader!\n");
		exit(10);
	}

	//release depth frame source
	SafeRelease(depthFrameSource);

	//allocate depth buffer
	m_depthBuffer = new uint16[m_depthWidth * m_depthHeight];
	
	//get color frame source
	IColorFrameSource* colorFrameSource;
	hr = m_sensor->get_ColorFrameSource(&colorFrameSource);
	if (FAILED(hr)){
		printf("Failed to get color frame source!\n");
		exit(10);
	}

	//get color frame reader
	hr = colorFrameSource->OpenReader(&m_colorFrameReader);
	if (FAILED(hr)){
		printf("Failed to open color frame reader!\n");
		exit(10);
	}

	//release the color frame source
	SafeRelease(colorFrameSource);

	//allocate color buffer
	m_colorBuffer = new uint32[1920 * 1080];
	
	//get the coordinate mapper
	hr = m_sensor->get_CoordinateMapper(&m_coordinateMapper);
	if (FAILED(hr)){
		printf("Failed to get coordinate mapper!\n");
		exit(10);
	}
	
	//allocate a buffer of color space points
	m_colorSpacePoints = new ColorSpacePoint[m_depthWidth * m_depthHeight];
}

void App::Tick(float deltaTime){
	//put update and drawing stuff here
	HRESULT hr;

	//depth stuff
	IDepthFrame* depthFrame;
	hr = m_depthFrameReader->AcquireLatestFrame(&depthFrame);
	if (FAILED(hr)) return;
	
	printf("Copying data!\n");
	hr = depthFrame->CopyFrameDataToArray(m_depthWidth * m_depthHeight, m_depthBuffer);

	if (FAILED(hr)){
		SafeRelease(depthFrame);
		printf("oh no, something went wrong while copying!\n");
		return;
	}
	
	SafeRelease(depthFrame);

	// Copia la profundidad de la data al screen  -- me veo azul
	int maximo = 0; int minimo = 1 << 30;
	for (int i = 0; i < m_depthWidth*m_depthHeight; i++){
		
		maximo = std::max(maximo, (int)m_depthBuffer[i]);
		minimo = std::min(minimo, (int)m_depthBuffer[i]);
		if (m_depthBuffer[i] > 800){
			m_pixelBuffer[i] = 0;
		}
		else{
			m_pixelBuffer[i] = m_depthBuffer[i];
		}
		
		/*
		uint16 d = m_depthBuffer[i];
		m_pixelBuffer[i] = d;
		uint8 a = d & 0xff;
		m_pixelBuffer[i] = (a << 16) | (a << 8) | a;// ahora grises
		*/
	}
	std::cout << "min max " << minimo << " " << maximo << std::endl;

	/*
	m_depthBuffer = new uint16[m_depthWidth*m_depthHeight];
	IColorFrameSource* colorFrameSource;
	hr=m_sensor->get_ColorFrameSource(&colorFrameSource);
	if (FAILED(hr)){
		printf("Failed to get Color from source\n");
		exit(10);
	}

	SafeRelease(colorFrameSource);
	*/

	
	//color stuff  le metemos color
	/*
	IColorFrame* colorFrame;
	hr = m_colorFrameReader->AcquireLatestFrame(&colorFrame);
	if (FAILED(hr))return;

	hr = colorFrame->CopyConvertedFrameDataToArray(1920 * 1080 * 4, (BYTE*)m_colorBuffer, ColorImageFormat_Bgra);
	if (FAILED(hr))return;
	
	SafeRelease(colorFrame);
	// Con esto genial ya tenemos a color la imagen :D
	for (int i = 0; i < 1920 * 1080; i++){
		m_pixelBuffer[i] = m_colorBuffer[i];
	}
	*/


	
	/*
	hr = m_coordinateMapper->MapDepthFrameToColorSpace(
		m_depthWidth * m_depthHeight, m_depthBuffer,
		m_depthWidth * m_depthHeight, m_colorSpacePoints);
	if (FAILED(hr))
	{
		printf("Oh no! Failed map the depth frame to color space!\n");
		return;
	}

	
	//copy depth data to the screen
	for (int i = 0; i < m_depthWidth * m_depthHeight; ++i)
	{
		ColorSpacePoint csp = m_colorSpacePoints[i];
		int ix = (int)csp.X;
		int iy = (int)csp.Y;

		if (ix >= 0 && ix < 1920 && iy >= 0 && iy < 1080)
			m_pixelBuffer[i] = m_colorBuffer[ix + iy * 1920];
		else
			m_pixelBuffer[i] = 0xff0000;

		//uint16 d = m_depthBuffer[i];
		//uint8 a = d & 0xff;
		//m_pixelBuffer[i] = (a << 16) | (a << 8) | a;
	}*/
	/**/
}

void App::Shutdown(){
	delete[] m_colorBuffer;
	SafeRelease(m_colorFrameReader);
	
	delete[] m_depthBuffer;
	SafeRelease(m_depthFrameReader);
	SafeRelease(m_sensor);
}