// HandDetectionOpenCV.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <opencv/cv.h> 
#include <opencv/highgui.h>
#include <opencv/ml.h>
#include <stdio.h>
#include <iostream>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <vector>
#include <sstream>
using namespace cv;
using namespace std;
int findBiggestContour(vector<vector<Point> >);

char ch[30];

//--------Using SURF as feature extractor and FlannBased for assigning a new point to the nearest one in the dictionary
Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("FlannBased");
Ptr<DescriptorExtractor> extractor = new SurfDescriptorExtractor();
SurfFeatureDetector detector(500);
//---dictionary size=number of cluster's centroids
int dictionarySize = 1500;
TermCriteria tc(CV_TERMCRIT_ITER, 10, 0.001);
int retries = 1;
int flags = KMEANS_PP_CENTERS;
BOWKMeansTrainer bowTrainer(dictionarySize, tc, retries, flags);
BOWImgDescriptorExtractor bowDE(extractor, matcher);
bool visited[5][300];
int tamtrain[] = { 50, 49, 68, 224 };
//int tameval[] = { 7, 9, 11, 21 };

void collectclasscentroids() {
	IplImage *img;
	int i, j;
	for (j = 1; j <= 4; j++)
		for (i = 0; i <= tamtrain[j - 1]; i++){
			cout << i << " " << j << endl;
			sprintf(ch, "%s%d%s%d%s", "train/", j, "_", i, ".png");
			const char* imageName = ch;
			img = cvLoadImage(imageName, 0);
			vector<KeyPoint> keypoint;
			detector.detect(img, keypoint);
			Mat features;
			extractor->compute(img, keypoint, features);
			try {

				bowTrainer.add(features);
			}
			catch (cv::Exception& e) { visited[j][i] = 1; }

		}
	return;
}



string tos(int t){
	stringstream st;
	st << t;
	return st.str();
}

int main(){

	IplImage *img2;
	cout << "Vector quantization..." << endl;
	memset(visited, 0, sizeof(visited));
	collectclasscentroids();


	vector<Mat> descriptors = bowTrainer.getDescriptors();
	int count = 0;
	for (vector<Mat>::iterator iter = descriptors.begin(); iter != descriptors.end(); iter++){
		count += iter->rows;
	}

	cout << "Clustering " << count << " features" << endl;
	//choosing cluster's centroids as dictionary's words
	Mat dictionary = bowTrainer.cluster();
	bowDE.setVocabulary(dictionary);
	cout << "extracting histograms in the form of BOW for each image " << endl;
	Mat labels(0, 1, CV_32FC1);
	Mat trainingData(0, dictionarySize, CV_32FC1);
	int k = 0;
	vector<KeyPoint> keypoint1;
	Mat bowDescriptor1;
	//extracting histogram in the form of bow for each image 
	for (int j = 1; j <= 4; j++)
		for (int i = 0; i <= tamtrain[j - 1]; i++){
			if (visited[j][i])continue;
			sprintf(ch, "%s%d%s%d%s", "train/", j, "_", i, ".png");
			const char* imageName = ch;
			img2 = cvLoadImage(imageName, 0);

			detector.detect(img2, keypoint1);
			bowDE.compute(img2, keypoint1, bowDescriptor1);
			trainingData.push_back(bowDescriptor1);
			labels.push_back((float)j);
		}



	//Setting up SVM parameters
	CvSVMParams params;
	params.kernel_type = CvSVM::RBF;
	params.svm_type = CvSVM::C_SVC;
	params.gamma = 0.50625000000000009;
	params.C = 312.50000000000000;
	params.term_crit = cvTermCriteria(CV_TERMCRIT_ITER, 100, 0.000001);
	CvSVM svm;
	printf("%s\n", "Training SVM classifier");

	bool res = svm.train(trainingData, labels, cv::Mat(), cv::Mat(), params);
	cout << "Processing evaluation data..." << endl;

	Mat groundTruth(0, 1, CV_32FC1);
	Mat evalData(0, dictionarySize, CV_32FC1);
	k = 0;
	vector<KeyPoint> keypoint2;
	Mat bowDescriptor2;

	VideoCapture cap(0);
	if (!cap.isOpened()){ cout << "Error en camara " << endl; return -1; }
	
	Mat src;// = imread("hand.png");
	Mat edges;
	//namedWindow("hand", 1);
	int cont = 0;

	while (true){
		cap >> src;
		//src = imread("hand.png");
		blur(src, src, Size(3, 3));
		cvtColor(src, edges, CV_BGR2GRAY);
		GaussianBlur(edges, edges, Size(7, 7), 1.5, 1.5);
		Canny(edges, edges, 0, 30, 3);


		int x = 150,
			y = 150,
			w = 180,
			h = 250;

		//Mat img1, img2;
		//img1 = imread("Lenna.png");
		Point P1(x, y);
		Point P2(x + w, y + h);

		// draw the box

		rectangle(src, P1, P2, Scalar(0, 255, 0));
		Rect roi(P1, P2);

		imshow("src", src);
		imshow("edges", edges);

		Mat smallImage = cv::Mat(edges, cv::Rect(x, y, w, h));
		imshow("img2", smallImage);

		if (cont > 10){

			detector.detect(smallImage, keypoint2);
			bowDE.compute(smallImage, keypoint2, bowDescriptor2);

			evalData.push_back(bowDescriptor2);
			float response = 4;
			try {
				response = svm.predict(bowDescriptor2);
			}
			catch (cv::Exception& e) { response = 4; }

			if (response == 1.0){
				cout << "DETENTE" << endl;
			}

			if (response == 2.0){
				cout << "GIRA DERECHA" << endl;
			}

			if (response == 3.0){
				cout << "AVANZA" << endl;
			}

			if (response == 4.0){
				cout << "NADA" << endl;
			}

			/*
			string namefile= "train/5/5_";
			namefile += tos( (cont-10) / 10);
			namefile += ".png";
			imwrite(namefile, smallImage);
			*/
		
		}

		/*
		Mat hsv;
		cvtColor(src, hsv, CV_BGR2HSV);

		Mat bw;
		inRange(hsv, Scalar(0, 10, 60), Scalar(20, 150, 255), bw);
		imshow("src", src);
		imshow("dst", bw);

		Mat canny_output;
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;

		findContours(bw, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
		int s = findBiggestContour(contours);

		Mat drawing = Mat::zeros(src.size(), CV_8UC1);
		drawContours(drawing, contours, s, Scalar(255), -1, 8, hierarchy, 0, Point());

		imshow("drw", drawing);
		*/
		cont++;
		if (waitKey(30) >= 0)break;
	}

	return 0;
}

int findBiggestContour(vector<vector<Point> > contours){
	int indexOfBiggestContour = -1;
	int sizeOfBiggestContour = 0;
	for (int i = 0; i < contours.size(); i++){
		if (contours[i].size() > sizeOfBiggestContour){
			sizeOfBiggestContour = contours[i].size();
			indexOfBiggestContour = i;
		}
	}
	return indexOfBiggestContour;
}