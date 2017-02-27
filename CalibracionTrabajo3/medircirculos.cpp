// PatternGrid.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include <iostream>
#include <queue>
#include <cmath>
#include <time.h>
#include <fstream>
#include <omp.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;
int acc[5000];
Mat src, grey;
int thresh = 100;
vector<Point>P;
bool visited[1001];
Size imageSize;
Size boardSize(4, 11);
double squareSize = 12.4;
string nombrevideo = "";

struct nodo{
	int id;// id del contorno
	int visit;// si es 1 entonces es un contorno fallido, 0 contorno valido
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

	if (a.visit != b.visit)return a.visit < b.visit;//no importante, elimina ambiguedad
	return a.id < b.id;
}

vector<int> get(int c[]){
	vector<int>dev;
	for (int i = 0; i < 75; i++)dev.push_back(c[i]);
	sort(dev.begin(), dev.end());
	return dev;
}

vector<int>takeframes(string namefile,int cant){
	vector<int>values;
	if (namefile == "PadronCirculos_01.avi"){
		int p1[75] = { 0, 5, 9, 12, 14, 20, 23, 24, 28, 33, 36, 40, 41, 54, 58, 78, 89, 90, 93, 94, 122, 125, 155, 159, 180, 182, 189, 195, 202, 243, 249, 312, 352, 407, 408, 415, 418, 420, 424, 485, 491, 494, 501, 509, 513, 516, 533, 536, 547, 550, 551, 557, 568, 596, 598, 610, 669, 691, 696, 697, 698, 700, 701, 703, 706, 722, 727, 799, 802, 805, 807, 809, 814, 815, 821 };
		values = get(p1);
	}

	if (namefile == "PadronCirculos_02.avi"){
		int p1[75] = { 19, 26, 36, 43, 54, 66, 70, 78, 106, 116, 120, 155, 174, 190, 214, 227, 252, 257, 267, 268, 271, 279, 303, 328, 329, 330, 340, 347, 351, 358, 377, 381, 387, 414, 428, 441, 445, 464, 474, 488, 490, 500, 506, 511, 539, 563, 578, 609, 660, 670, 685, 695, 708, 728, 812, 819, 830, 837, 852, 857, 863, 868, 869, 877, 884, 891, 904, 919, 927, 938, 998, 1009, 1021, 1034, 1043};
		//int p1[150] = { 19, 20, 22, 25, 26, 36, 37, 38, 42, 49, 66, 70, 78, 106, 107, 109, 116, 118, 120, 159, 174, 176, 190, 191, 196, 197, 203, 227, 252, 257, 258, 267, 271, 277, 279, 280, 282, 329, 330, 340, 347, 350, 351, 354, 358, 377, 381, 382, 387, 414, 418, 420, 428, 444, 445, 471, 473, 488, 490, 498, 500, 501, 504, 511, 515, 539, 577, 578, 644, 668, 670, 672, 682, 685, 688, 692, 693, 694, 695, 702, 708, 715, 716, 727, 728, 785, 806, 808, 809, 811, 812, 816, 819, 830, 831, 832, 833, 834, 835, 837, 838, 839, 849, 851, 852, 854, 856, 857, 858, 867, 868, 869, 873, 875, 878, 879, 884, 885, 886, 887, 888, 889, 890, 891, 892, 894, 904, 918, 919, 925, 938, 990, 998, 999, 1000, 1008, 1009, 1016, 1020, 1021, 1022, 1023, 1024, 1026, 1032, 1034, 1039, 1040, 1041, 1043 };
		values = get(p1);
	}

	if (namefile == "PadronCirculos_03.avi"){
		int p1[75] = { 2, 4, 11, 12, 13, 50, 60, 65, 71, 93, 94, 95, 187, 188, 200, 202, 205, 224, 228, 246, 252, 253, 265, 289, 291, 294, 302, 322, 323, 356, 365, 399, 405, 406, 410, 430, 432, 447, 449, 458, 491, 498, 530, 548, 562, 563, 568, 574, 599, 607, 614, 647, 661, 698, 715, 731, 749, 767, 769, 770, 793, 796, 798, 800, 805, 810, 815, 817, 821, 824, 859, 860, 905, 907};
		values = get(p1);
	}

	if (namefile == "calibrar_circulo_nuevo_1280x720.wmv"){
		int p1[75] = { 1, 3, 5, 8, 10, 13, 15, 175, 176, 177, 178, 179, 180, 181, 182, 183, 185, 188, 273, 274, 275, 279, 280, 283, 285, 287, 288, 290, 291, 293, 294, 296, 298, 300, 301, 303, 305, 307, 308, 310, 313, 314, 315, 316, 368, 369, 370, 372, 373, 374, 377, 381, 384, 389, 392, 399, 400, 404, 405, 406, 442, 443, 447, 454, 455, 456, 457, 458, 630, 631, 632, 639, 640, 641, 642};
		values = get(p1);
	}

	if (namefile == "calibrar_circulo_nuevo_640x360.wmv"){
		int p1[75] = { 21, 32, 36, 49, 50, 85, 122, 126, 129, 194, 198, 207, 218, 262, 281, 285, 308, 310, 312, 323, 324, 327, 391, 393, 456, 461, 468, 470, 476, 483, 490, 491, 499, 500, 506, 509, 514, 523, 532, 539, 542, 545, 546, 571, 580, 581, 596, 607, 608, 609, 611, 619, 622, 626, 627, 644, 646, 669, 670, 671, 682, 690, 692, 701, 708, 710, 711, 724, 727, 730, 749, 755, 761, 774, 777};
		values = get(p1);
	}

	if (namefile == "medir_1280x720_circulos.wmv"){
		int p1[75] = { 0, 3, 6, 9, 12, 14, 16, 18, 20, 22, 24, 26, 28, 29, 32, 34, 36, 40, 45, 46, 93, 95, 97, 99, 101, 103, 107, 145, 146, 148, 150, 152, 155, 158, 161, 163, 164, 168, 172, 176, 179, 182, 185, 188, 191, 194, 196, 197, 199, 200, 203, 204, 207, 209, 302, 303, 309, 310, 313, 316, 319, 324, 325, 326, 329, 333, 336, 340, 344, 345, 347, 350, 353, 354, 357 };
		//int p1[75] = { 0, 1, 4, 9, 10, 11, 13, 16, 18, 20, 21, 22, 23, 24, 26, 27, 29, 30, 31, 33, 34, 35, 36, 38, 40, 45, 95, 98, 99, 100, 101, 103, 151, 155, 157, 159, 161, 165, 167, 168, 170, 173, 174, 176, 179, 182, 183, 184, 187, 191, 192, 195, 197, 200, 202, 203, 313, 315, 316, 321, 322, 326, 327, 331, 334, 337, 340, 341, 344, 346, 348, 349, 351, 354, 357};
		values = get(p1);
	}

	if (namefile == "medir_640x360_circulos.wmv"){
		int p1[75] = { 0, 2, 3, 5, 7, 9, 10, 12, 15, 17, 18, 20, 21, 22, 25, 26, 28, 33, 36, 46, 52, 53, 54, 56, 61, 62, 64, 66, 69, 71, 74, 75, 76, 78, 79, 81, 83, 86, 88, 90, 91, 94, 96, 102, 104, 106, 113, 116, 122, 126, 130, 131, 133, 136, 137, 139, 141, 142, 144, 146, 147, 150, 152, 154, 157, 158, 160, 161, 163, 165, 171, 172, 173, 174, 176};
		values = get(p1);
	}


	//medir_640x360_circulos.wmv

	vector<int>dev;
	if (cant == 75){
		dev = values;
	}else if(cant==50){
		for (int i = 0; i < 75; i++)if (i % 3 == 0 || i % 3 == 2)dev.push_back(values[i]);
	}else if (cant == 25){
		for (int i = 0; i < 75; i++)if (i % 3 == 1 )dev.push_back(values[i]);
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

// area del triangulo
double area(pair<double, double>a, pair<double, double>b, pair<double, double>c){
	double x1 = a.first, y1 = a.second;
	double x2 = b.first, y2 = b.second;
	double x3 = c.first, y3 = c.second;
	return x2*y3 + x1*y2 + y1*x3 - x2*y1 - x3*y2 - y3*x1;
}

vector<nodo>v;

// consideramos solo los nodos considerados positivos
vector<nodo>limpiar(vector<nodo>v){
	vector<nodo>ans;
	for (int i = 0; i < v.size(); i++)
		if (v[i].visit == 0)
			ans.push_back(v[i]);
	return ans;
}

bool ordenpairy(pair<int, int>a, pair<int, int>b){
	if (a.second != b.second)return a.second < b.second;
	return a.first < b.first;
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
	vector<float>& perViewErrors)
{
	vector<Point2f> imagePoints2;
	int i, totalPoints = 0;
	double totalErr = 0, err;
	perViewErrors.resize(objectPoints.size());

	for (i = 0; i < (int)objectPoints.size(); ++i)
	{
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
	/*  esto servira para el padron des anillos
	for (int i = 0; i < boardSize.height; ++i)
	for (int j = 0; j < boardSize.width; ++j)
	corners.push_back(Point3f(float(j*squareSize), float(i*squareSize), 0));
	*/

	for (int i = 0; i < boardSize.height; i++)
		for (int j = 0; j < boardSize.width; j++)
			corners.push_back(Point3f(float((2 * j + i % 2)*squareSize), float(i*squareSize), 0));
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
		<< ". avg re projection error = " << totalAvgErr;

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
			corners.push_back(Point3f(float((2 * j + i % 2)*squareSize), float(i*squareSize), 0));

	Mat cameraMatrix = Mat::eye(3, 3, CV_64F);
	Mat distCoeffs = Mat::zeros(8, 1, CV_64F);

	if (nombrevideo == "medir_640x360_circulos.wmv"){
		cameraMatrix = (Mat_<double>(3, 3) << 474.7362331173333, 0.0, 336.9159183412804,
		0, 471.6468943643084, 170.6371361334853,
		0, 0, 1.0);
		distCoeffs = (Mat_<double>(8, 1) << -0.01566076683042372, 0.0357104158888431, 0
		, 0, -0.08833872766225771, 0, 0,0);
	}

	if (nombrevideo == "medir_1280x720_circulos.wmv"){
		cameraMatrix = (Mat_<double>(3, 3) << 1007.855168120671, 0, 652.9744034479452,
		0, 1012.871800701923, 348.6053600183893,
		0, 0, 1);
		distCoeffs = (Mat_<double>(8, 1) << 0.01446192806267651,
		0.3477133256230234,
		0,
		0,
		-1.339508273671581, 0, 0, 0);
	}

	//cout << "camera Matrix: " << cameraMatrix<< endl;
	//cout << "distCoeffs: " << distCoeffs<< endl;
	
	if (solvePnP(corners, pointBuf, cameraMatrix, distCoeffs, rvec, tvec)){
		Mat dst = Mat::zeros(3, 3, CV_64F);
		Rodrigues(rvec, dst, noArray());
		cout <<"distancia "<< norm((-dst).t()*tvec) << endl;
		string valuedist = "";
		valuedist = fdist(norm((-dst).t()*tvec));
		valuedist = "distancia " + valuedist;

		putText(original, valuedist, Point2f(15, 20), FONT_HERSHEY_PLAIN, 1, Scalar(0, 0, 255, 255));
	}
}

int main(){
	// colocar PadronCirculos_01.avi , PadronCirculos_02.avi , PadronCirculos_03.avi calibrar_circulo_nuevo_1280x720.wmv calibrar_circulo_nuevo_640x360.wmv, medir_1280x720_circulos.wmv,medir_640x360_circulos.wmv  que son los nombres de los videos
    nombrevideo = "medir_1280x720_circulos.wmv";
	int cantframes = 178;
	vector<int>frames;
	//vector<int>frames = takeframes(nombrevideo,cantframes);// los frames a considerar para la calibracion
	//vector<int>frames = primeroscantframes(cantframes);
	int p1[178] = { 22, 23, 35, 37, 39, 51, 52, 53, 54, 55, 57, 95, 96, 99, 139, 140, 143, 144, 145, 147, 148, 149, 218, 223, 234, 235, 236, 237, 245, 246, 247, 301, 317, 320, 322, 323, 324, 326, 327, 351, 354, 355, 356, 357, 369, 370, 373, 374, 376, 435, 436, 447, 448, 449, 455, 508, 519, 524, 525, 526, 527, 528, 530, 531, 532, 533, 535, 536, 539, 547, 548, 555, 556, 565, 566, 568, 569, 572, 573, 576, 577, 580, 582, 583, 584, 585, 586, 591, 592, 600, 601, 603, 604, 605, 607, 608, 609, 610, 611, 612, 613, 614, 624, 642, 643, 650, 654, 655, 657, 659, 670, 671, 672, 684, 685, 686, 688, 689, 690, 693, 694, 696, 697, 698, 699, 701, 702, 703, 706, 707, 708, 726, 728, 729, 754, 755, 756, 757, 759, 768, 777, 778, 779, 780, 790, 791, 792, 793, 796, 797, 799, 801, 803, 805, 818, 819, 821, 822, 826, 833, 845, 846, 848, 850, 851, 852, 853, 855, 857, 858, 859, 862, 863, 865, 866, 872, 873, 876 };
	for (int i = 0; i < cantframes; i++)
		frames.push_back(p1[i]);
	
	cout <<"cantidad de frames a procesar"<<cantframes<<" "<<frames.size() << endl;
	VideoCapture inputCapture(nombrevideo);
	int numframe = 0; // indica el numero de frame que esta siendo procesado en el video.
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

	while (true){
		Mat original;
		if (!inputCapture.read(original))break;
		imageSize = original.size();
		cout << imageSize << endl;

		//convertimos a escala de grises
		cvtColor(original, src, CV_BGR2GRAY);
		//eliminamos el ruido con desenfoque gaussiano
		GaussianBlur(src, src, Size(5, 5), 10, 10);
		//Obtenemos los bordes de la imagen con Canny
		Canny(src, src, thresh, 2 * thresh);

		Mat image2 = src;
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;
		//Encontramos todos los contornos de la imagen y lo guardamos en el vector Contours
		findContours(src, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

		/// Draw contours
		Mat drawing = Mat::zeros(src.size(), CV_8U);
		Mat dcircles = Mat::zeros(src.size(), CV_8U);
		// vector v contendra informacion de los *centros de los contornos validos
		v = vector<nodo>(contours.size(), nodo(1, 1, 1, 1));

		//paralelizamos esta parte ya que la funcion fitellipse es lenta
		#pragma omp parallel for 
		for (int ii = 0; ii < (int)contours.size(); ii++){
			vector<Point>P = contours[ii];
			Mat pointsf;
			RotatedRect box;
			// Un contorno es valido si el area es positiva y no sea tan grande
			if (contourArea(contours[ii])>1 && contourArea(contours[ii]) < 10000 && contours[ii].size()>5){
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
					//cout << box.center.x << " " << box.center.y << endl;
					Point2f centro = box.center;
					v[ii] = nodo(ii, 0, centro.x, centro.y);//valido
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
					v[j].visit = 1;
				}
			}
		}

		//nos quedamos con los validos
		v = limpiar(v);

		// Tenemos los centros de varios contornos y debemos quedarnos con los 44 del patron
		// La observacion esta en que la distancia entre cada centro del patron es pequeña
		// Por lo tanto dado una distancia debemos encontrar el tamaño maximo de las componentes conexas formadas (este debe ser 44)
		// Si una distancia es pequenha entonces la componente tendra tamanho 1 si es muy grande el tamanho sera el total de nodos.
		// Por lo tanto aplicamos binary search para encontrar la minima distancia de tal modo que la componente conexa mas grande sea 44 - tablero
		double lo = 1; double hi = 2000;
		vector<nodo>aux;
		for (int it = 0; it < 30; it++){
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
			if (maxi >= 44){
				hi = me;
			}
			else{
				lo = me;
			}
		}

		memset(visited, 0, sizeof(visited));
		// aqui solo recueperamos el vector con los 44 puntos deseados
		for (int i = 0; i < v.size(); i++){
			if (!visited[i]){
				aux.clear();
				int val = bfs(i, hi, aux);
				if (val == 44){
					break;
				}
			}
		}
		
		v = aux;

		/// Draw contours	
		for (int i = 0; i < v.size(); i++)
			drawContours(dcircles, contours, v[i].id, Scalar(255, 255, 255), CV_FILLED);
		// mostramos los contornos , el metodo descrito anteriormente es robusto ya que no se pierden los 44 puntos en la imagen
		imshow("circles", dcircles);

		/*    ---------------------------------------------------- 2DA PARTE  ---------------------------------------------*/
		//En la 2DA PARTE obtendremos el orden de los centros de los circulos ya encontrados esto para realizar el ZIGZAG

		vector<Point2f> pointBuf;
		if (v.size() != 44)continue;

		vector<Point2f> originalPoints;
		vector<Point2f> ch;  // Convex hull points

		for (int i = 0; i < v.size(); i++)
			originalPoints.push_back(Point2f(v[i].cx, v[i].cy));
		// Realizamos el convex hull de los 44 puntos
		convexHull(Mat(originalPoints), ch, false);
		vector<nodo>nododir;

		// luego de obtener el convex hull, queremos obtener la "horizontal" de la parte superior del padron
		// Debido a la forma del poligono convexo habran 2 vertices que forman angulo de 90 grados
		// Para encontrar los vertices obtenemos el seno del angulo utilizando el producto vectorial
		for (int i = 0; i < ch.size(); i++){
			Point p1 = Point(ch[(i + 1) % ch.size()].x - ch[i].x, ch[(i + 1) % ch.size()].y - ch[i].y);
			Point p2 = Point(ch[i].x - ch[(i - 1 + ch.size()) % ch.size()].x, ch[i].y - ch[(i - 1 + ch.size()) % ch.size()].y);
			double prod = abs(p1.x*p2.y - p1.y*p2.x);
			prod /= hypot(p1.x, p1.y);
			prod /= hypot(p2.x, p2.y);

			for (int j = 0; j < v.size(); j++){
				// sen(X)= |AxB|/(|A||B|)
				// si sen(X)=1 el angulo es de 90 grados
				if (prod >= 0.9){
					if (hypot(v[j].cx - ch[i].x, v[j].cy - ch[i].y) <= 0.1){
						nododir.push_back(v[j]);
					}
				}
			}
		}

		// luego de obtener el segmento "horizontal" encontraremos su pendiente
		// Despues debemos encontrar los 8 segmentos paralelos
		if (nododir.size() != 2)continue;
		double arcotan = -1;
		float mindx = 1 << 20; float mindy = 1 << 20;

		if (nododir[0].cx < nododir[1].cx){
			arcotan = atan((nododir[1].cy - nododir[0].cy) / (0.0 + nododir[1].cx - nododir[0].cx));
			mindx = nododir[0].cx;
			mindy = nododir[0].cy;
		}
		else{
			arcotan = atan((nododir[0].cy - nododir[1].cy) / (0.0 + nododir[0].cx - nododir[1].cx));
			mindx = nododir[1].cx;
			mindy = nododir[1].cy;
		}

		int sizetam[] = { 8, 8, 7, 6, 5, 4, 3, 2, 1 };

		vector<vector<nodo> >vnodo;
		for (int caso = 0; caso < 7; caso++){
			vector<pair<double, int> >puntos4;
			double mindis = 1 << 30;
			for (int i = 0; i < v.size(); i++){
				if (v[i].visit != 0)continue;
				if (caso == 0)
					if (v[i].cx != mindx || v[i].cy != mindy)continue;

				for (int j = 0; j < v.size(); j++){
					if (v[j].visit != 0)continue;
					// si la recta tiene una pendiente cercana a lo buscado
					if (i != j && v[i].cx < v[j].cx){

						vector<pair<double, int> >dis;
						for (int k = 0; k < v.size(); k++){
							if (v[k].visit != 0)continue;
							if (k == i || k == j)continue;
							double aux = area(make_pair(v[i].cx, v[i].cy), make_pair(v[j].cx, v[j].cy), make_pair(v[k].cx, v[k].cy)) / hypot(v[i].cx - v[j].cx, v[i].cy - v[j].cy);
							aux = abs(aux);
							// si los puntos son colineales
							dis.push_back(make_pair(aux, k));
						}

						sort(dis.begin(), dis.end());
						double sum = 0;
						for (int k = 0; k + 2 < sizetam[caso]; k++)
							sum += dis[k].first;

						if (mindis > sum){
							mindis = sum;
							puntos4.clear();
							puntos4.push_back(make_pair(v[i].cx, i));
							puntos4.push_back(make_pair(v[j].cx, j));
							for (int k = 0; k + 2 < sizetam[caso]; k++)
								puntos4.push_back(make_pair(v[dis[k].second].cx, dis[k].second));
						}

					}
				}

			}

			vector<nodo>aux;
			sort(puntos4.begin(), puntos4.end());
			// en caso de encontrar el segmento, lo marcamos como visitados
			for (int i = 0; i < puntos4.size(); i++){
				v[puntos4[i].second].visit = caso + 1;
				aux.push_back(v[puntos4[i].second]);
			}

			if (aux.size() > 0)
				vnodo.push_back(aux);
		}

		// hay 3 puntos faltantes, considerar la diagonal de tamanho 2 la que tenga la menor distancia
		double mindiag = 1 << 30;
		vector<nodo>aux2;
		int pid1 = -1, pid2 = -1;

		for (int i = 0; i < v.size(); i++)
			if (v[i].visit == 0)
				for (int j = 0; j < v.size(); j++)
					if (v[j].visit == 0 && v[i].cx<v[j].cx){
						if (hypot(v[i].cx - v[j].cx, v[i].cy - v[j].cy) < mindiag){
							mindiag = hypot(v[i].cx - v[j].cx, v[i].cy - v[j].cy);
							pid1 = i; pid2 = j;
						}
					}

		v[pid1].visit = 1;
		v[pid2].visit = 1;
		aux2.push_back(v[pid1]);
		aux2.push_back(v[pid2]);
		vnodo.push_back(aux2);

		// el ultimo punto no es considerado
		aux2.clear();
		for (int i = 0; i < v.size(); i++)
			if (v[i].visit == 0)
				aux2.push_back(v[i]);
		if (aux2.size() > 0)
			vnodo.push_back(aux2);

		// ordenamos los segmentos de abajo hacia arriba mediante el producto vectorial
		for (int i = 0; i < vnodo.size(); i++){
			for (int j = i + 1; j < vnodo.size(); j++){
				pair<double, double>vectora = make_pair(vnodo[i][vnodo[i].size() - 1].cx - vnodo[i][0].cx, vnodo[i][vnodo[i].size() - 1].cy - vnodo[i][0].cy);
				pair<double, double>vectorb = make_pair(vnodo[j][vnodo[j].size() - 1].cx - vnodo[i][0].cx, vnodo[j][vnodo[j].size() - 1].cy - vnodo[i][0].cy);
				if (vectora.first*vectorb.second - vectora.second*vectorb.first >= 0){
					swap(vnodo[i], vnodo[j]);
				}
			}
		}

		// hay 9 segmentos que corresponden a 9 diagonales de tamanho 8 , 8 , 7 ,6 ,5,4,3,2,1
		// los cuales estan ordenados de arriba hacia abajo por producto vectorial
		// lo que deseamos es tener 11 verticales de tamanho 4, por lo que realizamos el siguiente procedimiento
		// obtener el menor x de cada segmento de abajo hacia arriba 1 por cada uno hasta llegar a 4. marcarlo como visitado
		if (vnodo.size() == 9){
			bool visited2[9][100];
			memset(visited2, 0, sizeof(visited2));
			// si es una diagonal de tamanho 2 entonces ordenar de menor a mayor x
			if (vnodo[0].size() == 2){
				for (int caso = 0; caso < 11; caso++){
					int count = 0;
					for (int i = 0; i < 9 && count < 4; i++){
						for (int j = 0; j < vnodo[i].size() && count < 4; j++){
							if (!visited2[i][j]){
								pointBuf.push_back(Point2f(vnodo[i][j].cx, vnodo[i][j].cy));
								count++;
								visited2[i][j] = 1;
								break;
							}
						}
					}
				}
			}// si es una diagonal de tamanho 1 entonces ordenar de mayor a menor x (es decir leer viceversa)
			else{
				for (int caso = 0; caso < 11; caso++){
					int count = 0;
					for (int i = 0; i < 9 && count < 4; i++){
						for (int j = vnodo[i].size() - 1; j >= 0 && count < 4; j--){
							if (!visited2[i][j]){
								pointBuf.push_back(Point2f(vnodo[i][j].cx, vnodo[i][j].cy));
								count++;
								visited2[i][j] = 1;
								break;
							}
						}
					}
				}
			}
		}

		cout << "siiii" << pointBuf.size() << endl;
		bool ok = 0;

		if (pointBuf.size() == 44){
			if (pointBuf[0].x < pointBuf[43].x)
				reverse(pointBuf.begin(), pointBuf.end());
			//reverse(pointBuf.begin(), pointBuf.end());
			//drawChessboardCorners(original, boardSize, Mat(pointBuf), 1);
			ok = 1;
			/*
			if (imagePoints.size() <= 1){
			cout << "empieza" << endl;
			for (int i = 0; i < pointBuf.size(); i++)
			cout << pointBuf[i].x << " " << pointBuf[i].y << endl;
			cout << "fin " << endl;
			}
			*/

			////
			cout << "numFrame: " << numframe <<" "<<posframe<<" "<<frames.size()<< endl;
			if (nombrevideo == "medir_1280x720_circulos.wmv" || nombrevideo == "medir_640x360_circulos.wmv"){
				distancia(original, pointBuf);
				
			}
			else{
				string namesave = "calibrar_circulo_nuevo_640x360/file_";
				/**/
				if (numframe == frames[posframe]){
					string ww = namesave + fdist(contx++);
					ww += ".jpg";
					imwrite(ww,original);
					imagePoints.push_back(pointBuf);
					drawChessboardCorners(original, boardSize, Mat(pointBuf), 1);
					imshow("super original", original);
					posframe++;
				}


				if (posframe == frames.size()){
					cout << "inicia calibracion" << endl;
					runCalibrationAndSave(original, imageSize, cameraMatrix, distCoeffs, imagePoints);
					break;
				}
				
				/*
				double error = runCalibrationAndSave(original, imageSize, cameraMatrix, distCoeffs, vector<vector<Point2f> >(1, pointBuf));
				cout << "error " << error << endl;
				errorframe.push_back(make_pair(error, numframe));
				*/
				numframe++;

			}
			
		}


		namedWindow("Original", WINDOW_AUTOSIZE);
		imshow("Original", original);

		waitKey(1);
		//if (cont==1)break;
	}

	sort(errorframe.begin(), errorframe.end());
	myfile << "errores" << endl;
	for (int i = 0; i < errorframe.size();i++)
		myfile << errorframe[i].first << "," << errorframe[i].second<<endl;

	return 0;
}


