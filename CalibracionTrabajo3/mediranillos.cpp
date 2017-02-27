// PatternGrid.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include <iostream>
#include <omp.h>
#include <queue>
#include <time.h>
#include <fstream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <time.h>
#include <fstream>
using namespace std;
using namespace cv;
Mat src, grey;
int thresh = 100;
bool visited[1001];
Size imageSize;
Size boardSize(5, 6);
double squareSize = 25.5;
string nombrevideo;
int numframe = 0;
vector<double>dist;
bool getthreshold = 0;
Mat H;
Mat Hinv;
Mat finalrectify;

// estructura que almacena la informacion de contorno
struct nodo{
	int id;// id del contorno
	int visit;// si es 1 entonces es un contorno fallido y se debe borrar con el metodo limpiar(), 0 contorno valido
	float cx;// centro x
	float cy;// centro y
	nodo(int _id, int _visit, float _cx, float _cy){
		id = _id; visit = _visit;
		cx = _cx; cy = _cy;
	}
};

// define el metodo de ordenamiento sort en un nodo
bool operator<(nodo a, nodo b){
	if (a.cx != b.cx)return a.cx < b.cx;//los menores centros en x
	if (a.cy != b.cy)return a.cy < b.cy;//sino los menores centros en y
	if (a.visit != b.visit)return a.visit < b.visit;//prioridad los .visit=0 es decir los validos
	return a.id < b.id;
}

vector<int> get(int c[]){
	vector<int>dev;
	for (int i = 0; i < 75; i++)dev.push_back(c[i]);
	sort(dev.begin(), dev.end());
	return dev;
}

vector<int>takeframes(string namefile, int cant){
	vector<int>values;
	if (namefile == "PadronAnillos_01.avi"){
		int p1[75] = { 2, 3, 5, 6, 7, 24, 134, 140, 141, 144, 146, 148, 149, 150, 151, 153, 156, 159, 162, 163, 166, 170, 171, 196, 209, 217, 432, 460, 462, 488, 504, 508, 511, 517, 520, 524, 529, 537, 545, 594, 605, 628, 631, 643, 644, 656, 659, 670, 692, 702, 728, 739, 742, 746, 747, 750, 754, 793, 796, 798, 808, 843, 855, 875, 884, 897, 903, 956, 991, 995, 1036, 1039, 1046, 1078, 1085 };
		values = get(p1);
	}

	if (namefile == "PadronAnillos_02.avi"){
		int p1[75] = { 1, 6, 19, 22, 24, 27, 30, 33, 40, 46, 57, 61, 79, 84, 101, 125, 136, 146, 147, 153, 157, 161, 241, 261, 266, 271, 284, 298, 376, 406, 407, 416, 512, 517, 530, 548, 560, 592, 713, 722, 730, 746, 763, 768, 775, 787, 796, 803, 811, 912, 927, 994, 995, 1030, 1158, 1188, 1194, 1200, 1203, 1216, 1263, 1310, 1318, 1329, 1363, 1383, 1403, 1404, 1461, 1486, 1540, 1565, 1603, 1689, 1727 };
		values = get(p1);
	}

	if (namefile == "PadronAnillos_03.avi"){
		int p1[75] = { 36, 41, 89, 95, 154, 174, 222, 239, 284, 293, 368, 397, 436, 459, 472, 495, 533, 546, 595, 600, 630, 657, 677, 702, 725, 751, 809, 827, 895, 946, 977, 997, 1038, 1043, 1097, 1114, 1152, 1174, 1231, 1257, 1291, 1310, 1390, 1412, 1471, 1504, 1541, 1555, 1596, 1606, 1644, 1702, 1723, 1749, 1794, 1807, 1849, 1898, 1902, 1922, 1990, 2007, 2035, 2073, 2096, 2102, 2157, 2193, 2203, 2248, 2293, 2346, 2366, 2401, 2427 };
		values = get(p1);
	}


	if (namefile == "calibrar_anillo_nuevo_1280x720.wmv"){
		//int p1[75] = { 1, 3, 5, 7, 9, 12, 15, 18, 20, 22, 24, 26, 28, 30, 33, 36, 38, 40, 42, 44, 46, 50, 54, 58, 102, 191, 195, 198, 202, 205, 208, 212, 214, 217, 221, 230, 239, 243, 245, 248, 249, 250, 254, 255, 256, 257, 258, 259, 260, 261, 264, 318, 320, 392, 393, 394, 418, 419, 431, 433, 434, 450, 453, 455, 474, 477, 479, 480, 482, 486, 508, 519, 543, 565, 567};
		int p1[91] = { 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 28, 30, 32, 35, 38, 40, 42, 44, 46, 48, 50, 52, 54, 58, 102, 191, 193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223, 230, 238, 239, 241, 243, 245, 247, 249, 250, 253, 255, 257, 259, 262, 264, 318, 320, 321, 391, 393, 394, 396, 416, 419, 431, 433, 435, 450, 452, 454, 456, 474, 477, 479, 481, 482, 484, 486, 506, 508, 509, 519, 540, 543, 564, 567 };
		values = get(p1);
	}

	if (namefile == "calibrar_anillo_nuevo_640x360.wmv"){
		int p1[75] = { 2, 7, 13, 19, 32, 36, 49, 76, 95, 107, 129, 130, 140, 165, 172, 181, 203, 227, 247, 261, 282, 299, 301, 302, 312, 320, 344, 347, 353, 371, 375, 385, 407, 429, 442, 456, 466, 485, 497, 513, 525, 541, 551, 562, 574, 585, 596, 618, 619, 636, 668, 670, 688, 733, 783, 791, 815, 835, 861, 868, 900, 994, 1002, 1016, 1114, 1118, 1211, 1330, 1368, 1382, 1414, 1461, 1481, 1488, 1491 };
		values = get(p1);
	}
	
	if (namefile == "medir_1280x720_anillos.wmv"){
		//int p1[75] = { 0, 1, 3, 4, 7, 10, 11, 14, 16, 18, 19, 20, 22, 24, 25, 26, 30, 31, 32, 65, 66, 70, 74, 76, 77, 108, 109, 110, 122, 125, 126, 127, 128, 133, 137, 140, 142, 143, 148, 153, 154, 156, 160, 166, 168, 170, 173, 175, 178, 179, 186, 189, 190, 191, 193, 196, 199, 201, 206, 246, 248, 250, 251, 253, 257, 259, 260, 261, 263, 264, 266, 271, 275, 276, 278};
		int p1[75] = { 1, 3, 5, 9, 12, 14, 16, 25, 29, 30, 33, 35, 61, 62, 64, 69, 70, 71, 74, 75, 76, 77, 78, 108, 110, 119, 120, 121, 122, 123, 126, 128, 132, 135, 138, 144, 145, 146, 147, 148, 152, 153, 154, 155, 159, 164, 166, 169, 170, 173, 174, 175, 176, 179, 181, 184, 186, 188, 191, 194, 197, 200, 201, 202, 205, 241, 243, 247, 251, 253, 260, 264, 269, 273, 274 };
		values = get(p1);
	}

	if (namefile == "medir_640x360_anillos.wmv"){
		int p1[75] = { 0, 29, 47, 55, 69, 74, 75, 83, 89, 91, 96, 97, 99, 101, 108, 113, 115, 126, 129, 132, 135, 140, 146, 152, 158, 164, 168, 169, 182, 188, 200, 204, 208, 212, 222, 228, 237, 251, 252, 264, 286, 290, 363, 367, 376, 378, 386, 390, 397, 400, 401, 405, 413, 419, 423, 430, 439, 446, 453, 461, 467, 476, 479, 489, 497, 500, 506, 512, 517, 521, 533, 541, 548, 607, 613 };
		values = get(p1);
	}

	vector<int>dev;
	if (cant == 75){
		dev = values;
	}
	else if (cant == 50){
		for (int i = 0; i < 75; i++)if (i % 3 == 0 || i % 3 == 2)dev.push_back(values[i]);
	}
	else if (cant == 25){
		for (int i = 0; i < 75; i++)if (i % 3 == 1)dev.push_back(values[i]);
	}
	else{
		for (int i = 0; i < cant; i++)dev.push_back(values[i]);
	}
	return dev;
}

vector<int> primeroscantframes(int cant){
	vector<int>dev;
	for (int i = 0; i < cant; i++)dev.push_back(i);
	return dev;
}

// consideramos solo los nodos considerados positivos
vector<nodo>limpiar(vector<nodo>v){
	vector<nodo>ans;
	for (int i = 0; i < v.size(); i++)
		if (v[i].visit == 0)
			ans.push_back(v[i]);
	return ans;
}

vector<nodo>v;


//hallar el area de un triangulo 
double area(pair<double, double>a, pair<double, double>b, pair<double, double>c){
	double x1 = a.first, y1 = a.second;
	double x2 = b.first, y2 = b.second;
	double x3 = c.first, y3 = c.second;
	return abs(x2*y3 + x1*y2 + y1*x3 - x2*y1 - x3*y2 - y3*x1);
}

int bfs(int pos, double dis, vector<nodo>&aux){
	visited[pos] = 1;
	queue<int>Q;
	Q.push(pos);// metemos a la cola el nodo actual
	int cont = 1;
	aux.push_back(v[pos]);

	while (!Q.empty()){
		int id = Q.front();
		Q.pop();
		for (int i = 0; i < v.size(); i++){
			// si hay algun nodo no visitado que se encuentra a una distancia dada entonces
			// lo agregamos a la cola
			if (!visited[i] && hypot(v[i].cx - v[id].cx, v[i].cy - v[id].cy) <= dis){
				aux.push_back(v[i]);
				cont++;
				Q.push(i);
				visited[i] = 1;
			}
		}
	}

	return cont;
}


double computeReprojectionErrors(const vector<vector<Point3f> >& objectPoints,
	const vector<vector<Point2f> >& imagePoints,
	const vector<Mat>& rvecs, const vector<Mat>& tvecs,
	const Mat& cameraMatrix, const Mat& distCoeffs,
	vector<float>& perViewErrors){
	vector<Point2f> imagePoints2;
	int i, totalPoints = 0;
	double totalErr = 0, err;
	perViewErrors.resize(objectPoints.size());

	for (i = 0; i < (int)objectPoints.size(); ++i){
		projectPoints(Mat(objectPoints[i]), rvecs[i], tvecs[i], cameraMatrix,
			distCoeffs, imagePoints2);
		err = norm(Mat(imagePoints[i]), Mat(imagePoints2), CV_L2);

		int n = (int)objectPoints[i].size();
		perViewErrors[i] = (float)std::sqrt(err*err / n);
		totalErr += err*err;
		totalPoints += n;
	}

	return std::sqrt(totalErr / totalPoints);
}


void calcBoardCornerPositions(Size boardSize, float squareSize, vector<Point3f>& corners){
	corners.clear();
	//  esto servira para el padron des anillos
	for (int i = 0; i < boardSize.height; ++i)
		for (int j = 0; j < boardSize.width; ++j)
			corners.push_back(Point3f(float(j*squareSize), float(i*squareSize), 0));
	
}

bool runCalibration(Mat &view, Size& imageSize, Mat& cameraMatrix, Mat& distCoeffs,
	vector<vector<Point2f> > imagePoints, vector<Mat>& rvecs, vector<Mat>& tvecs,
	vector<float>& reprojErrs, double& totalAvgErr){

	cameraMatrix = Mat::eye(3, 3, CV_64F);
	//if (s.flag & CV_CALIB_FIX_ASPECT_RATIO)cameraMatrix.at<double>(0, 0) = 1.0;

	distCoeffs = Mat::zeros(8, 1, CV_64F);

	vector<vector<Point3f> > objectPoints(1);
	calcBoardCornerPositions(boardSize, squareSize, objectPoints[0]);

	objectPoints.resize(imagePoints.size(), objectPoints[0]);

	//Find intrinsic and extrinsic camera parameters
	double rms = calibrateCamera(objectPoints, imagePoints, imageSize, cameraMatrix,
		distCoeffs, rvecs, tvecs, CV_CALIB_ZERO_TANGENT_DIST);

	cout << "Re-projection error reported by calibrateCamera: " << rms << endl;

	bool ok = checkRange(cameraMatrix) && checkRange(distCoeffs);

	totalAvgErr = computeReprojectionErrors(objectPoints, imagePoints,
		rvecs, tvecs, cameraMatrix, distCoeffs, reprojErrs);

	return ok;
}


// Print camera parameters to the output file
void saveCameraParams(Mat view, Size& imageSize, Mat& cameraMatrix, Mat& distCoeffs,
	const vector<Mat>& rvecs, const vector<Mat>& tvecs,
	const vector<float>& reprojErrs, const vector<vector<Point2f> >& imagePoints,
	double totalAvgErr){
	//FileStorage fs(s.outputFileName, FileStorage::WRITE);
	ofstream myfile;
	myfile.open("salida.txt");
	//myfile << "Writing this to a file.\n";

	if (!rvecs.empty() || !reprojErrs.empty())
		cout << "nrOfFrames" << (int)std::max(rvecs.size(), reprojErrs.size());
	myfile << "image_Width" << imageSize.width << endl;
	myfile << "image_Height" << imageSize.height << endl;
	myfile << "board_Width" << boardSize.width << endl;
	myfile << "board_Height" << boardSize.height << endl;
	myfile << "square_Size" << squareSize << endl;

	//if (s.flag & CV_CALIB_FIX_ASPECT_RATIO)
	//	cout << "FixAspectRatio" << s.aspectRatio;
	/*
	if (s.flag){
	sprintf(buf, "flags: %s%s%s%s",
	s.flag & CV_CALIB_USE_INTRINSIC_GUESS ? " +use_intrinsic_guess" : "",
	s.flag & CV_CALIB_FIX_ASPECT_RATIO ? " +fix_aspectRatio" : "",
	s.flag & CV_CALIB_FIX_PRINCIPAL_POINT ? " +fix_principal_point" : "",
	s.flag & CV_CALIB_ZERO_TANGENT_DIST ? " +zero_tangent_dist" : "");
	//cvWriteComment(*fs, buf, 0);

	}*/

	//cout << "flagValue" << s.flag;


	myfile << "Camera_Matrix" << cameraMatrix << endl;
	myfile << "Distortion_Coefficients" << distCoeffs << endl;
	myfile << "Avg_Reprojection_Error" << totalAvgErr << endl;

	if (!reprojErrs.empty())
		myfile << "Per_View_Reprojection_Errors" << Mat(reprojErrs) << endl;

	if (!rvecs.empty() && !tvecs.empty()){
		CV_Assert(rvecs[0].type() == tvecs[0].type());
		Mat bigmat((int)rvecs.size(), 6, rvecs[0].type());
		for (int i = 0; i < (int)rvecs.size(); i++)
		{
			Mat r = bigmat(Range(i, i + 1), Range(0, 3));
			Mat t = bigmat(Range(i, i + 1), Range(3, 6));

			CV_Assert(rvecs[i].rows == 3 && rvecs[i].cols == 1);
			CV_Assert(tvecs[i].rows == 3 && tvecs[i].cols == 1);
			//*.t() is MatExpr (not Mat) so we can use assignment operator
			r = rvecs[i].t();
			t = tvecs[i].t();
		}
		//cvWriteComment(*fs, "a set of 6-tuples (rotation vector + translation vector) for each view", 0);
		myfile << "Extrinsic_Parameters" << bigmat << endl;
	}

	if (!imagePoints.empty()){
		Mat imagePtMat((int)imagePoints.size(), imagePoints[0].size(), CV_32FC2);
		for (int i = 0; i < (int)imagePoints.size(); i++){
			Mat r = imagePtMat.row(i).reshape(2, imagePtMat.cols);
			Mat imgpti(imagePoints[i]);
			imgpti.copyTo(r);
		}
		myfile << "Image_points" << imagePtMat << endl;
	}

	myfile.close();
}


double runCalibrationAndSave(Mat view, Size imageSize, Mat&  cameraMatrix, Mat& distCoeffs, vector<vector<Point2f> > imagePoints){
	vector<Mat> rvecs, tvecs;
	vector<float> reprojErrs;
	double totalAvgErr = 0;

	bool ok = runCalibration(view, imageSize, cameraMatrix, distCoeffs, imagePoints, rvecs, tvecs,
		reprojErrs, totalAvgErr);
	cout << (ok ? "Calibration succeeded" : "Calibration failed")
		<< ". avg re projection error = " << totalAvgErr<<endl;

	if (ok)
		saveCameraParams(view, imageSize, cameraMatrix, distCoeffs, rvecs, tvecs, reprojErrs,
		imagePoints, totalAvgErr);
	return totalAvgErr;
}


string fdist(double d){
	stringstream st;
	st << d;
	return st.str();
}

void distancia(Mat &original, vector<Point2f> pointBuf){
	
	cv::Mat rvec(3, 1, CV_64F);
	cv::Mat tvec(3, 1, CV_64F);
	vector<Point3f> corners;

	for (int i = 0; i < boardSize.height; i++)
		for (int j = 0; j < boardSize.width; j++)
			corners.push_back(Point3f(float(j*squareSize), float(i*squareSize), 0));

	Mat cameraMatrix = Mat::eye(3, 3, CV_64F);
	Mat distCoeffs = Mat::zeros(8, 1, CV_64F);

	if (nombrevideo == "medir_640x360_anillos.wmv"){
		cameraMatrix = (Mat_<double>(3, 3) << 488.0321954583645, 0, 330.2765238775828,
		0, 485.5352079137646, 172.0592124328771,
		0, 0, 1);
		distCoeffs = (Mat_<double>(8, 1) << -0.004211367806464171,
			-0.00825232247307175,
		0,
		0,
		-0.04104284377324537, 0, 0, 0);
		
	}

	if (nombrevideo == "medir_1280x720_anillos.wmv"){
		cameraMatrix = (Mat_<double>(3, 3) << 8931.729453468537, 0, 845.740330746426,
		0, 14857.81881131294, 377.2998192120343,
		0, 0, 1);
		distCoeffs = (Mat_<double>(8, 1) << -89.31293800581575,
		13812.14756155819,
		0,
		0,
		51.4529672655805, 0, 0, 0);
	}

	//cout << "camera Matrix: " << cameraMatrix<< endl;
	//cout << "distCoeffs: " << distCoeffs<< endl;

	if (solvePnP(corners, pointBuf, cameraMatrix, distCoeffs, rvec, tvec)){
		Mat dst = Mat::zeros(3, 3, CV_64F);
		Rodrigues(rvec, dst, noArray());
		double daux = norm((-dst).t()*tvec);
		dist.push_back(daux);
		//cout << "distancia " << daux << endl;
		string valuedist = "";
		valuedist = fdist(daux);
		valuedist = "distancia " + valuedist;

		putText(original, valuedist, Point2f(15, 20), FONT_HERSHEY_PLAIN, 1, Scalar(0, 0, 255, 255));
		putText(original, fdist(numframe), Point2f(50,50), FONT_HERSHEY_PLAIN, 1, Scalar(0, 0, 255, 255));
	}
}


vector<Point2f>nuevosPuntos(Mat image,double inix,double iniy){
	Mat filter;
	cvtColor(image, filter, CV_BGR2GRAY);
	//GaussianBlur(filter, filter, Size(5, 5), 10, 10);
	//if (getthreshold)adaptiveThreshold(filter, filter, 255.0, CV_THRESH_BINARY, CV_THRESH_BINARY, 101, -5);
	Canny(filter, filter, thresh, 2 * thresh);

	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;
	//Encontramos todos los contornos de la imagen y lo guardamos en el vector Contours
	findContours(filter, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
	v = vector<nodo>(contours.size(), nodo(1, 1, 1, 1));
	for (int ii = 0; ii < (int)contours.size(); ii++){
		vector<Point>P = contours[ii];
		Mat pointsf;
		RotatedRect box;
		if (contourArea(contours[ii])>1 && contourArea(contours[ii])<10000 && contours[ii].size()>5){
			Mat pointsf;RotatedRect box;
			try{
				Mat(contours[ii]).convertTo(pointsf, CV_32F);
				box = fitEllipse(pointsf);
			}
			catch (cv::Exception& e){
				continue;
			}

			float f1 = box.size.height;
			float f2 = box.size.width;
			if (fabs(f2 - f1) <= 12 && f1*f2 <= 1500){
				Point2f centro = box.center;
				v[ii] = nodo(ii, 0, centro.x, centro.y);
			}
		}
	}

	v = limpiar(v);
	sort(v.begin(), v.end());
	set<nodo>S;
	memset(visited, 0, sizeof(visited));

	for (int i = 0; i < v.size(); i++){
		if (visited[i])continue;
		double sumx = 0; double sumy = 0;
		int cont = 0;
		for (int j = i + 1; j < v.size(); j++){
			if (visited[j])continue;
			if (abs((v[i].cx - v[j].cx)) <= 1.5  && abs(v[i].cy - v[j].cy) <= 1.5){
				visited[j] = 1;
				sumx += v[i].cx;
				sumy += v[i].cy;
				cont++;
			}
			if (v[j].cx > v[i].cx + 10)break;
		}

		if (cont >= 1){
			v[i].cx = sumx / cont;
			v[i].cy = sumy / cont;
			S.insert(v[i]);
		}
	}

	v = vector<nodo>(S.begin(), S.end());
	//luego de eliminar duplicados con Set , si hay 2 centros cercanos solo se considerara 1 de ellos
	for (int i = 0; i < v.size(); i++)
		for (int j = i + 1; j < v.size(); j++)
			if (abs((v[i].cx - v[j].cx)) <= 1.5 && abs(v[i].cy - v[j].cy) <= 1.5){
				v[j].visit = 1;
			}


	//eliminaremos todos los elementos que tengan el atributo visit diferente de 0
	v = limpiar(v);
	double lo = 1; double hi = 65;
	if (imageSize.width > 700)hi = 120;


	vector<nodo>aux;
	for (int it = 0; it < 20; it++){
		double me = (lo + hi) / 2;
		memset(visited, 0, sizeof(visited));
		int maxi = 0;
		for (int i = 0; i < v.size(); i++){
			if (!visited[i]){
				aux.clear();
				// el bfs me retorna el tamanho de la componente conexa. la busqueda se hace en anchura 
				// la busqueda en profundidad tambien es valido
				maxi = max(maxi, bfs(i, me, aux));
			}
		}
		if (maxi >= 30){
			hi = me;
		}
		else{
			lo = me;
		}
	}


	memset(visited, 0, sizeof(visited));
	// aqui solo recuperamos el vector con los 30 puntos deseados
	for (int i = 0; i < v.size(); i++){
		if (!visited[i]){
			aux.clear();
			int val = bfs(i, hi, aux);
			if (val == 30){
				break;
			}
		}
	}

	v = aux;
	//cout << "sera pe v.size() " << v.size()<<endl;

	if (v.size() == 30){
		vector<vector<nodo> >vsegmento;
		for (int caso = 0; caso < 5; caso++){
			double minarea = 1e+10;
			vector<pair<double, int> >puntos6;
			for (int i = 0; i < v.size(); i++){
				if (v[i].visit != 0)continue;
				for (int j = 0; j < v.size(); j++){
					if (v[j].visit != 0)continue;
					if (i != j && v[i].cx <= v[j].cx){
						vector<pair<double, int> >dis;
						for (int k = 0; k < v.size(); k++){
							if (v[k].visit != 0)continue;
							if (k == i || k == j)continue;
							double aux = area(make_pair(v[i].cx, v[i].cy), make_pair(v[j].cx, v[j].cy), make_pair(v[k].cx, v[k].cy)) / hypot(v[i].cx - v[j].cx, v[i].cy - v[j].cy);
							dis.push_back(make_pair(aux, k));
						}

						sort(dis.begin(), dis.end());
						double area4 = dis[0].first + dis[1].first + dis[2].first + dis[3].first;
						if (minarea > area4){
							minarea = area4;
							puntos6.clear();
							puntos6.push_back(make_pair(v[i].cx, i));
							puntos6.push_back(make_pair(v[dis[0].second].cx, dis[0].second));
							puntos6.push_back(make_pair(v[dis[1].second].cx, dis[1].second));
							puntos6.push_back(make_pair(v[dis[2].second].cx, dis[2].second));
							puntos6.push_back(make_pair(v[dis[3].second].cx, dis[3].second));
							puntos6.push_back(make_pair(v[j].cx, j));
						}
					}
				}
			}

			vector<nodo>aux;
			sort(puntos6.begin(), puntos6.end());
			for (int i = 0; i < puntos6.size(); i++){
				v[puntos6[i].second].visit = caso + 1;
				aux.push_back(v[puntos6[i].second]);
			}

			vsegmento.push_back(aux);
		}

		for (int i = 0; i < vsegmento.size(); i++){
			for (int j = i + 1; j < vsegmento.size(); j++){
				pair<double, double>vectora = make_pair(vsegmento[i][5].cx - vsegmento[i][0].cx, vsegmento[i][5].cy - vsegmento[i][0].cy);
				pair<double, double>vectorb = make_pair(vsegmento[j][5].cx - vsegmento[i][0].cx, vsegmento[j][5].cy - vsegmento[i][0].cy);
				if (vectora.first*vectorb.second - vectora.second*vectorb.first < 0){
					swap(vsegmento[i], vsegmento[j]);
				}
			}
		}

		vector<Point2f> pointBuf, pointBuf2(30, Point2f(0, 0));

		// el metodo anterior ordena los segmentos de arriba hacia abajo usando ordenamiento burbuja
		for (int i = 0; i < 5; i++){  // 5 segmentos ordenados de abajo hacia arriba
			for (int j = 0; j < vsegmento[i].size(); j++){ // cada segmento de tamanho 6 esta ordenado por x
				pointBuf.push_back(Point2f(vsegmento[i][j].cx, vsegmento[i][j].cy));
			}
		}

		int cont = 0;
		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 5; j++)
				pointBuf2[cont++] = pointBuf[i + j * 6];

		if (pointBuf2.size() == 30)
			if (pointBuf2[0].x < pointBuf2[29].x){
				reverse(pointBuf2.begin(), pointBuf2.end());
			}

		//cout << "Image Size() " << image.size() << endl;

		drawChessboardCorners(image, boardSize, Mat(pointBuf2), 1);
		/*
		cout << "Inicio frontal" << endl;
		for (int i = 0; i < 30; i++)
			cout << pointBuf2[i].x << " " << pointBuf2[i].y << endl;
		cout << "Fin frontal" << endl;
		*/
		imshow("Frontal Calibracion", image);

		Mat nuevote;
		Mat homoinv = Hinv;

		//cout << "homoinv " << homoinv.size() << endl;
		vector<Point2f> pointBuf3(30, Point2f(0, 0));

		/*
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++)
				cout << H.at<double>(i, j) << " ";
			cout << endl;
		}

		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++)
				cout << homoinv.at<double>(i, j) << " ";
			cout << endl;
		}
		*/
		for (int i = 0; i < 30; i++){
			double sum = (double)(homoinv.at<double>(0, 0))*(pointBuf2[i].x) + (double)(homoinv.at<double>(0, 1))*(pointBuf2[i].y) + homoinv.at<double>(0, 2);
			sum /= ((double)(homoinv.at<double>(2, 0))*(pointBuf2[i].x) + (double)(homoinv.at<double>(2, 1))*(pointBuf2[i].y) + (double)(homoinv.at<double>(2, 2)));

			pointBuf3[i].x = sum;
			sum = (double)(homoinv.at<double>(1, 0))*(pointBuf2[i].x) + (double)(homoinv.at<double>(1, 1))*(pointBuf2[i].y) + homoinv.at<double>(1, 2);
			sum /= ((double)(homoinv.at<double>(2, 0))*(pointBuf2[i].x) + (double)(homoinv.at<double>(2, 1))*(pointBuf2[i].y) + (double)(homoinv.at<double>(2, 2)));
			pointBuf3[i].y = sum;
		}

		/*
		cout << "Inicio Reproyeccion" << endl;
		for (int i = 0; i < 30; i++){
			cout << pointBuf3[i].x << " " << pointBuf3[i].y << endl;
		}
		cout << "Fin Nueva projection" << endl;
		*/
		warpPerspective(image, nuevote, Hinv, imageSize);
		imshow("Reproyeccion Imagen", nuevote);

		drawChessboardCorners(finalrectify, boardSize, Mat(pointBuf3), 1);
		
		imshow("Reproyeccion puntos Final", finalrectify);
		return pointBuf3;
	}



	vector<Point2f>ans;
	return ans;
}


void ReDistortPoints(vector<Point2f> & src, vector<Point2f> & dst,Mat & cameraMatrix, const cv::Mat & distorsionMatrix){
	dst.clear();
	double fx = cameraMatrix.at<double>(0, 0);
	double fy = cameraMatrix.at<double>(1, 1);
	double ux = cameraMatrix.at<double>(0, 2);
	double uy = cameraMatrix.at<double>(1, 2);

	double k1 = distorsionMatrix.at<double>(0);
	double k2 = distorsionMatrix.at<double>(1);
	double p1 = distorsionMatrix.at<double>(2);
	double p2 = distorsionMatrix.at<double>(3);
	double k3 = distorsionMatrix.at<double>(4);

	for(int i = 0; i < src.size(); i++){
		Point2f pt = src[i];
		double x = (pt.x-ux)/fx;
		double y = (pt.y-uy)/fy;
		double xCorrected, yCorrected;
		double r2 = x*x + y*y;
		//radial distorsion
		xCorrected = x * (1. + k1 * r2 + k2 * r2 * r2 + k3 * r2 * r2*r2);
		yCorrected = y * (1. + k1 * r2 + k2 * r2 * r2 + k3 * r2 * r2*r2);

		//tangential distorsion
		xCorrected = xCorrected + (2. * p1 * x * y + p2 * (r2 + 2. * x * x));
		yCorrected = yCorrected + (p1 * (r2 + 2. * y * y) + 2. * p2 * x * y);
		xCorrected = xCorrected * fx + ux;
		yCorrected = yCorrected * fy + uy;
		dst.push_back(Point2f(xCorrected, yCorrected));
	}
}

void MyDistortPoints(vector<Point2f> & src, vector<Point2f> & dst,
	const cv::Matx33d & cameraMatrix, const cv::Matx<double, 1, 8> & distorsionMatrix){
	cv::Mat cameraMatrix2(cameraMatrix);
	cv::Mat distorsionMatrix2(distorsionMatrix);
	return ReDistortPoints(src, dst, cameraMatrix2, distorsionMatrix2);
}

int main(){
	// colocar PadronAnillos_01.avi , PadronAnillos_02.avi , PadronAnillos_03.avi,calibrar_anillo_nuevo_1280x720.wmv,calibrar_anillo_nuevo_640x360.wmv,medir_1280x720_anillos,medir_640x360_anillos.wmv que son los nombres de los videos
	nombrevideo = "medir_1280x720_anillos.wmv";
	int cantframes = 87;

	//vector<int>frames = takeframes(nombrevideo, cantframes);// los frames a considerar para la calibracion
	//int p1[91] = { 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 28, 30, 32, 35, 38, 40, 42, 44, 46, 48, 50, 52, 54, 58, 102, 191, 193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223, 230, 238, 239, 241, 243, 245, 247, 249, 250, 253, 255, 257, 259, 262, 264, 318, 320, 321, 391, 393, 394, 396, 416, 419, 431, 433, 435, 450, 452, 454, 456, 474, 477, 479, 481, 482, 484, 486, 506, 508, 509, 519, 540, 543, 564, 567 };
	int p1[87] = { 1, 7, 12, 19, 31, 36, 51, 59, 76, 87, 95, 107, 129, 130, 140, 172, 178, 181, 203, 208, 211, 241, 260, 294, 302, 308, 310, 348, 354, 373, 377, 383, 405, 435, 453, 475, 481, 487, 496, 513, 516, 524, 527, 533, 541, 547, 551, 562, 574, 586, 587, 595, 596, 606, 615, 619, 634, 636, 653, 656, 676, 792, 824, 830, 847, 859, 861, 865, 876, 989, 998, 1011, 1036, 1215, 1226, 1343, 1377, 1379, 1380, 1391, 1427, 1434, 1478, 1494, 1499, 1501, 1504 };

	vector<int>frames;
	for (int i = 0; i < cantframes; i++)frames.push_back(p1[i]);

	//vector<int>frames = primeroscantframes(cantframes);
	cout << "cantidad de frames a procesar" << cantframes << " " << frames.size() << endl;
	VideoCapture inputCapture(nombrevideo);
	numframe = 0; // indica el numero de frame que esta siendo procesado en el video.
	int posframe = 0; // sirve para posicionar en el arreglo
	setNumThreads(8);
	/*    ---------------------------------------------------- 1ERA PARTE  ---------------------------------------------*/
	//1ERA PARTE del trabajo corresponde en obtener el tablero de anillos (total 30) sin ruidos
	vector<vector<Point2f> > imagePoints;
	ofstream myfile;
	myfile.open("vector.txt");
	vector<pair<double, int> >errorframe;

	int contx = 0;
	Mat cameraMatrix = Mat::eye(3, 3, CV_64F);
	Mat distCoeffs = Mat::zeros(8, 1, CV_64F);


	Mat undi, rview, map1, map2;

	while (true){
		Mat original;
		// Mat original contiene los frames del video a colores
		if (!inputCapture.read(original))break;
		imageSize = original.size();

		//convertimos a escala de grises
		cvtColor(original, src, CV_BGR2GRAY);
		//eliminamos el ruido con desenfoque gaussiano
		GaussianBlur(src, src, Size(5, 5), 10, 10);
		//Obtenemos los bordes de la imagen con Canny
		if (imageSize.width < 700){
			if (nombrevideo != "PadronAnillos_03.avi"){
				getthreshold = 1;
				adaptiveThreshold(src, src, 255.0, CV_THRESH_BINARY, CV_THRESH_BINARY, 101, -5);
				//imshow("Adaptative", src);
			}
		}
		Canny(src, src, thresh, 2 * thresh);
		//imshow("Canny", src);

		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;
		//Encontramos todos los contornos de la imagen y lo guardamos en el vector Contours
		findContours(src, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

		/// Draw contours
		Mat drawing2 = Mat::zeros(src.size(), CV_8UC3);
		Mat drawing = Mat::zeros(src.size(), CV_8U);
		Mat image2 = Mat::zeros(src.size(), CV_8U);
		Mat image3 = Mat::zeros(src.size(), CV_8U);
		Mat image4 = Mat::zeros(src.size(), CV_8U);

		// vector v contendra informacion de los centros de los contornos validos
		v = vector<nodo>(contours.size(), nodo(1, 1, 1, 1));

		//paralelizamos esta parte ya que la funcion fitellipse es lenta
		#pragma omp parallel for 
		for (int ii = 0; ii < (int)contours.size(); ii++){
			vector<Point>P = contours[ii];
			Mat pointsf;
			RotatedRect box;
			// Un contorno es valido si el area es positiva y no sea tan grande
			if (contourArea(contours[ii])>1 && contourArea(contours[ii])<10000 && contours[ii].size()>5){
				Mat pointsf;
				RotatedRect box;
				//utilizamos fitellipse para obtener un rectangulo que lo contenga
				try{
					Mat(contours[ii]).convertTo(pointsf, CV_32F);
					box = fitEllipse(pointsf);
				}
				catch (cv::Exception& e){
					// en caso se produzca excepcion, si no se pueda formar la elipse con los puntos del contorno
					continue;
				}

				float f1 = box.size.height;
				float f2 = box.size.width;
				// si el largo y el ancho difieren minimamente en tamaño y el area no es tan grande sera considerado
				// como posible contorno
				if (fabs(f2 - f1) <= 12 && f1*f2 <= 1500){
					Point2f centro = box.center;
					v[ii] = nodo(ii, 0, centro.x, centro.y);
				}
			}
		}

		
		v = limpiar(v);
		
		sort(v.begin(), v.end());
		set<nodo>S;
		// ahora trabajaremos con los centros de cada contorno
		// tenemos circulos concentricos si existen 2 puntos iguales o muy cercanos  solo consideraremos 1 de ellos
		// Set nos permite eliminar duplicados
		memset(visited, 0, sizeof(visited));
		
		for (int i = 0; i < v.size(); i++){
			if (visited[i])continue;
			double sumx = 0; double sumy = 0;
			int cont = 0;
			for (int j = i + 1; j < v.size(); j++){
				if (visited[j])continue;
				if (abs((v[i].cx - v[j].cx)) <= 1.5  && abs(v[i].cy - v[j].cy) <= 1.5){
					visited[j] = 1;
					sumx += v[i].cx;
					sumy += v[i].cy;
					cont++;
				}
				if (v[j].cx > v[i].cx + 10)break;
			}

			if (  cont>=1){
				v[i].cx = sumx / cont;
				v[i].cy = sumy / cont;
				S.insert(v[i]);
			}
		}

		v = vector<nodo>(S.begin(), S.end());
		//luego de eliminar duplicados con Set , si hay 2 centros cercanos solo se considerara 1 de ellos
		for (int i = 0; i < v.size(); i++)
			for (int j = i + 1; j < v.size(); j++)
				if (abs((v[i].cx - v[j].cx)) <= 1.5 && abs(v[i].cy - v[j].cy) <= 1.5){
					v[j].visit = 1;
				}
		

		//eliminaremos todos los elementos que tengan el atributo visit diferente de 0
		v = limpiar(v);
		/*
		cout << "circle2 " << v.size() << endl;
		for (int i = 0; i < v.size(); i++)
			circle(image3, Point2f(v[i].cx, v[i].cy), 1, Scalar(255, 255, 255));
		imshow("circle2 ", image3);
		*/

		// Tenemos los centros de varios contornos y debemos quedarnos con los 30 del patron
		// La observacion esta en que la distancia entre cada centro del patron es pequeña
		// Por lo tanto dado una distancia debemos encontrar el tamaño maximo de las componentes conexas formadas (este debe ser 30)
		// Si una distancia es pequenha entonces la componente tendra tamanho 1 si es muy grande el tamanho sera el total de nodos.
		// Por lo tanto aplicamos binary search para encontrar la minima distancia de tal modo que la componente conexa mas grande sea 44 - tablero
		double lo = 1; double hi = 65;
		if (imageSize.width > 700)hi = 120;

		
		vector<nodo>aux;
		for (int it = 0; it < 20; it++){
			double me = (lo + hi) / 2;
			memset(visited, 0, sizeof(visited));
			int maxi = 0;
			for (int i = 0; i < v.size(); i++){
				if (!visited[i]){
					aux.clear();
					// el bfs me retorna el tamanho de la componente conexa. la busqueda se hace en anchura 
					// la busqueda en profundidad tambien es valido
					maxi = max(maxi, bfs(i, me, aux));
				}
			}
			if (maxi >= 30){
				hi = me;
			}
			else{
				lo = me;
			}
		}


		memset(visited, 0, sizeof(visited));
		// aqui solo recuperamos el vector con los 30 puntos deseados
		for (int i = 0; i < v.size(); i++){
			if (!visited[i]){
				aux.clear();
				int val = bfs(i, hi, aux);
				if (val == 30){
					break;
				}
			}
		}

		v = aux;

		// dibujamos los contornos que solo corresponden al tablero de anillos (total 30)
		for (int i = 0; i < v.size(); i++)
			drawContours(drawing, contours, v[i].id, Scalar(255, 255, 255), CV_FILLED);

		/*    ---------------------------------------------------- 2DA PARTE  ---------------------------------------------*/
		//En la 2DA PARTE obtendremos el orden de los centros de los anillos estos para realizar el ZIGZAG
		namedWindow("Tablero", CV_WINDOW_AUTOSIZE);
		imshow("Tablero", drawing);

		//El tablero es de 5X6
		Size boardSize(5, 6);
		//cout << "elementos= " << v.size() << endl;

		// en caso de que se haya filtrado los 30 elementos correctamente hacer lo siguiente
		if (v.size() == 30){
			vector<vector<nodo> >vsegmento;
			// vamos a obtener 5 segmentos los cuales contienen 6 puntos cada uno
			for (int caso = 0; caso < 5; caso++){
				double minarea = 1e+10;
				vector<pair<double, int> >puntos6;

				for (int i = 0; i < v.size(); i++){
					if (v[i].visit != 0)continue;
					for (int j = 0; j < v.size(); j++){
						if (v[j].visit != 0)continue;
						if (i != j && v[i].cx <= v[j].cx){
							vector<pair<double, int> >dis;
							// fijamos los puntos i , j los cuales vienen hacer los extremos del segmento que debe incluir 4 puntos
							for (int k = 0; k < v.size(); k++){
								if (v[k].visit != 0)continue;
								if (k == i || k == j)continue;
								// hallamos la distancia de un punto al segmento formado por los puntos i , j   
								// distancia= area/base
								double aux = area(make_pair(v[i].cx, v[i].cy), make_pair(v[j].cx, v[j].cy), make_pair(v[k].cx, v[k].cy)) / hypot(v[i].cx - v[j].cx, v[i].cy - v[j].cy);
								dis.push_back(make_pair(aux, k));
							}

							sort(dis.begin(), dis.end());
							// ordenamos las menores distancias, los 4 primeros puntos vienen hacer los puntos dentros del segmento i,j
							// cabe mencionar de que la distancia deberia ser 0 ya que los puntos buscados deben pertenecer a la recta
							// pero por errores de precision y/o distorsion esto no se da por ellos la busqueda de las 4 menores distancias
							double area4 = dis[0].first + dis[1].first + dis[2].first + dis[3].first;
							if (minarea > area4){
								minarea = area4;
								puntos6.clear();
								puntos6.push_back(make_pair(v[i].cx, i));
								puntos6.push_back(make_pair(v[dis[0].second].cx, dis[0].second));
								puntos6.push_back(make_pair(v[dis[1].second].cx, dis[1].second));
								puntos6.push_back(make_pair(v[dis[2].second].cx, dis[2].second));
								puntos6.push_back(make_pair(v[dis[3].second].cx, dis[3].second));
								puntos6.push_back(make_pair(v[j].cx, j));
							}
						}
					}
				}

				vector<nodo>aux;
				// ordenamiento por x
				sort(puntos6.begin(), puntos6.end());
				// luego de obtener la recta con 6 puntos lo marcamos como visitados
				// Para que no obtener la misma recta en la siguiente iteracion
				for (int i = 0; i < puntos6.size(); i++){
					v[puntos6[i].second].visit = caso + 1;
					aux.push_back(v[puntos6[i].second]);
				}

				vsegmento.push_back(aux);
			}

			// hemos obtenido 6 segmentos de 5 puntos pero no existe un orden entre ellos
			// Un metodo robusto para definir el orden es utilizar el producto vectorial

			for (int i = 0; i < vsegmento.size(); i++){
				for (int j = i + 1; j < vsegmento.size(); j++){
					pair<double, double>vectora = make_pair(vsegmento[i][5].cx - vsegmento[i][0].cx, vsegmento[i][5].cy - vsegmento[i][0].cy);
					pair<double, double>vectorb = make_pair(vsegmento[j][5].cx - vsegmento[i][0].cx, vsegmento[j][5].cy - vsegmento[i][0].cy);

					// dado 2 vectores a y b (cola menor x , cabeza mayor x)  
					// el vector a es un segmento i
					// el vector b es la cola del segmento i y la cabeza del segmento j
					// si el producto vectorial es positivo entonces segmento i esta debajo del segmento j
					if (vectora.first*vectorb.second - vectora.second*vectorb.first < 0){
						swap(vsegmento[i], vsegmento[j]);
					}
				}
			}

			vector<Point2f> pointBuf, pointBuf2(30, Point2f(0, 0));

			// el metodo anterior ordena los segmentos de arriba hacia abajo usando ordenamiento burbuja
			for (int i = 0; i < 5; i++){  // 5 segmentos ordenados de abajo hacia arriba
				for (int j = 0; j < vsegmento[i].size(); j++){ // cada segmento de tamanho 6 esta ordenado por x
					pointBuf.push_back(Point2f(vsegmento[i][j].cx, vsegmento[i][j].cy));
				}
			}

			// con pointBuf tenemos el orden en horizontal, para colocarlo en vertical realizamos lo siguiente
			int cont = 0;
			for (int i = 0; i < 6; i++)
				for (int j = 0; j < 5; j++)
					pointBuf2[cont++] = pointBuf[i + j * 6];

			if (pointBuf2.size() == 30)
				if (pointBuf2[0].x < pointBuf2[29].x){
					reverse(pointBuf2.begin(), pointBuf2.end());
				}

			/*
			cout << "Inicio PointBuf2" << endl;
			for (int i = 0; i < 30; i++)
				cout << pointBuf2[i] << endl;
			cout << "Fin PointBuf2" << endl;
			*/

			//Ya no quiero hallar las distancias<< endl;
			if (nombrevideo == "medir_1280x720_anillos.wmv" || nombrevideo == "medir_640x360_anillos.wmv"){
				distancia(original, pointBuf2);
				numframe++;
			}
			else{
				string namesave = "calibrar_anillo_nuevo_640x360/file_";
				// esta seccion comentada me permite calibrar segun frames que se encuentren en el arreglo frames se crea el archivo salida con los parametros de camara
				/*
				if (numframe == frames[posframe]){
					string ww = namesave + fdist(contx++);
					ww += ".jpg";
					imwrite(ww,original);
					imagePoints.push_back(pointBuf2);
					posframe++;
				}

				if (posframe == frames.size()){
					cout << "inicia calibracion" << endl;
					runCalibrationAndSave(original, imageSize, cameraMatrix, distCoeffs, imagePoints);
					break;
				}
				*/

				// esta seccion te genera el archivo vector con los RMS - independiente por cada frame para la seleccion de frames - se pone en el arreglo int p1[75] de la funcion takeframes()
				/*
				double error = runCalibrationAndSave(original, imageSize, cameraMatrix, distCoeffs, vector<vector<Point2f> >(1, pointBuf2));
				cout << "error " << error << endl;
				errorframe.push_back(make_pair(error, numframe));
				*/
				numframe++;

			}

			/*
			cout << "Inicio Sin distorsion" << endl;
			for (int i = 0; i < pointBuf2.size(); i++)
				cout << pointBuf2[i].x << " " << pointBuf2[i].y << endl;
			cout << "Fin sin distorsion" << endl;
			*/

			namedWindow("Original", WINDOW_AUTOSIZE);
			imshow("Original", original);
			//Generamos los frames sin distorsion
			undistort(original, undi, cameraMatrix, distCoeffs);
			imshow("Sin distorsion", undi);
			
			original = undi.clone();


			vector< Point2f> P1; 
			vector< Point2f> P2(4);

			double dx = (pointBuf2[0].x - pointBuf2[25].x) / 20;
			double dy = (pointBuf2[25].y - pointBuf2[29].y) / 20;
			dx *= 4;
			dy *= 4;
			//cout << dx << " " << dy << endl;

			P1.push_back(Point2f(pointBuf2[29].x - dx, pointBuf2[29].y - dy));
			P1.push_back(Point2f(pointBuf2[4].x + dx, pointBuf2[4].y - dy));
			P1.push_back(Point2f(pointBuf2[0].x + dx, pointBuf2[0].y + dy));
			P1.push_back(Point2f(pointBuf2[25].x - dx, pointBuf2[25].y + dy));


			P2[0].x = 0; P2[0].y = 0;
			P2[1].x = squareSize*10; P2[1].y = 0;
			P2[2].x = squareSize*10; P2[2].y = squareSize*8;
			P2[3].x = 0; P2[3].y = squareSize*8;
			//Mat H = findHomography(P1, P2);
			H = findHomography(P1, P2);
			//warping
			Mat warped_image;
			//Mat original2;
			Mat original3;
			finalrectify = original.clone();
			warpPerspective(original, warped_image, H, cv::Size(squareSize * 10, squareSize * 8));
			
			imshow("Frontal", warped_image);
			Hinv = findHomography(P2, P1);

			//cout << "Size Nuevos puntos warped_image " << warped_image.size() << endl;
			vector<Point2f> pointBufx = nuevosPuntos(warped_image, pointBuf2[29].x - dx, pointBuf2[29].y - dy);
			vector<Point2f> redistortedPoints;
			//cout << "nuevos puntos " << pointBufx.size() << endl;
			ReDistortPoints(pointBufx, redistortedPoints, cameraMatrix, distCoeffs);
			/*
			cout << "inicio reproy" << endl;
			for (int i = 0; i < redistortedPoints.size(); i++)
				cout << pointBufx[i].x << " " << pointBufx[i].y << endl;
			cout << "fin " << endl;
			*/

			//warpPerspective(warped_image, original2, H.inv(), imageSize);
			
			drawChessboardCorners(original, boardSize, Mat(pointBuf2), 1);
			imshow("Original Calibracion", original);

			//imshow("View", original2);

		}

		// dibujo del tablero sin ruido
		waitKey(1);
	}

	sort(errorframe.begin(), errorframe.end());
	myfile << "errores" << endl;
	for (int i = 0; i < errorframe.size(); i++)
		myfile << errorframe[i].first << "," << errorframe[i].second << endl;

	ofstream myfiledist;
	myfiledist.open("distancia.txt");
	for (int i = 0; i < dist.size();i++)
		myfiledist << dist[i] << endl;


	return 0;
}