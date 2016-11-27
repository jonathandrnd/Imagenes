// EdgeDetection_CPU.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <omp.h>
#include <ctime>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string>
#include <ctime>
#include <omp.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace std;
using namespace cv;

float convolutionKernelStore[256];

void convolve(unsigned char *source, int width, int height, int paddingX, int paddingY, unsigned int kOffset, int kWidth, int kHeight, unsigned char *destination){
	// Calculate our pixel's location
	//int x = (blockIdx.x * blockDim.x) + threadIdx.x;
	//int y = (blockIdx.y * blockDim.y) + threadIdx.y;
	int   pWidth = kWidth / 2;
	int   pHeight = kHeight / 2;

	#pragma omp parallel for
	for (int i = pWidth; i < width - pWidth; i++)
		for (int j = pHeight; j < height - pHeight; j++){
			float sum = 0.0;
			for (int dj = -pHeight; dj <= pHeight; dj++){
				for (int di = -pWidth; di <= pWidth; di++){
					int ki = (di + pWidth);
					int kj = (dj + pHeight);
					float w = convolutionKernelStore[(kj * kWidth) + ki + kOffset];
					sum += w * float(source[((j + dj) * width) + (i + di)]);
				}
			}
			destination[(j * width) + i] = (unsigned char)sum;
		}
}



void pythagoras(unsigned char *a, unsigned char *b, unsigned char *c,int width,int height){	
	#pragma omp parallel for
	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++){
			float af = float(a[j*width+i]);
			float bf = float(b[j*width+i]);
			c[j*width + i] = (unsigned char)min(sqrtf(af*af + bf*bf),(float)255);
		}
}

int main(){
	omp_set_num_threads(8);
	cv::Mat          frame;
	// Create the capture windows
	const float gaussianKernel5x5[25] = {
		2.f / 159.f, 4.f / 159.f, 5.f / 159.f, 4.f / 159.f, 2.f / 159.f,
		4.f / 159.f, 9.f / 159.f, 12.f / 159.f, 9.f / 159.f, 4.f / 159.f,
		5.f / 159.f, 12.f / 159.f, 15.f / 159.f, 12.f / 159.f, 5.f / 159.f,
		4.f / 159.f, 9.f / 159.f, 12.f / 159.f, 9.f / 159.f, 4.f / 159.f,
		2.f / 159.f, 4.f / 159.f, 5.f / 159.f, 4.f / 159.f, 2.f / 159.f,
	};

	const unsigned int gaussianKernel5x5Offset = 0;

	// Sobel gradient kernels
	const float sobelGradientX[9] = {
		-1.f, 0.f, 1.f,
		-2.f, 0.f, 2.f,
		-1.f, 0.f, 1.f,
	};
	const float sobelGradientY[9] = {
		1.f, 2.f, 1.f,
		0.f, 0.f, 0.f,
		-1.f, -2.f, -1.f,
	};

	const unsigned int sobelGradientXOffset = sizeof(gaussianKernel5x5) / sizeof(float);
	const unsigned int sobelGradientYOffset = sizeof(sobelGradientX) / sizeof(float) + sobelGradientXOffset;
	frame = cvLoadImage("C:\\Users\\Lenovo\\Desktop\\MAESTRIA\\lenna.jpg");

	int width = frame.size().width;
	int height = frame.size().height;

	unsigned char *sourceDataDevice = (unsigned char *)malloc(width*height);
	unsigned char *blurredDataDevice = (unsigned char *)malloc(width*height);
	unsigned char *edgesDataDevice = (unsigned char *)malloc(width*height);
	
	unsigned char *deviceGradientX, *deviceGradientY;

	//unsigned char *sourceDataDevice, *blurredDataDevice, *edgesDataDevice;
	cv::Mat source(frame.size(), CV_8U,  sourceDataDevice);
	cv::Mat blurred(frame.size(), CV_8U,  blurredDataDevice);
	cv::Mat edges(frame.size(), CV_8U, edgesDataDevice);
	deviceGradientX = (unsigned char *)malloc(width*height);
	deviceGradientY = (unsigned char *)malloc(width*height);

	for (int i = 0; i < 5; i++)for (int j = 0; j < 5; j++)convolutionKernelStore[i * 5 + j] = gaussianKernel5x5[i * 5 + j];
	for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++)convolutionKernelStore[i * 3 + j + sobelGradientXOffset] = sobelGradientX[i * 3 + j];
	for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++)convolutionKernelStore[i * 3 + j + sobelGradientYOffset] = sobelGradientY[i * 3 + j];

	cvtColor(frame, source, CV_BGR2GRAY);

	clock_t ini = clock();

	convolve(sourceDataDevice, width, height, 0, 0, gaussianKernel5x5Offset, 5, 5,blurredDataDevice);
	convolve(blurredDataDevice, width, height, 2, 2, sobelGradientXOffset, 3, 3, deviceGradientX);
	convolve(blurredDataDevice, width, height, 2, 2, sobelGradientYOffset, 3, 3, deviceGradientY);
	pythagoras(deviceGradientX, deviceGradientY, edgesDataDevice, width, height);

	clock_t fin = clock();
	std::cout << "Elapsed CPU time: " << fin-ini << " milliseconds" << std::endl;


	imshow("Source", frame);
	imshow("Greyscale", source);
	imshow("Blurred", blurred);
	imshow("Sobel", edges);

	waitKey(0);
	return 0;
}

