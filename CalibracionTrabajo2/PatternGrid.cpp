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
		int p1[75] = { 41, 44, 48, 96, 100, 161, 173, 221, 230, 240, 289, 292, 376, 388, 441, 492, 500, 539, 543, 587, 601, 611, 615, 657, 660, 670, 676, 680, 731, 752, 811, 851, 873, 896, 923, 964, 988, 1002, 1050, 1102, 1111, 1172, 1211, 1249, 1263, 1318, 1349, 1389, 1404, 1461, 1524, 1569, 1603, 1610, 1716, 1733, 1754, 1801, 1846, 1863, 1883, 1907, 1930, 1949, 2005, 2050, 2107, 2113, 2151, 2207, 2209, 2334, 2367, 2406, 2431};
		values = get(p1);
	}


	if (namefile == "calibrar_anillo_nuevo_1280x720.wmv"){
		int p1[75] = { 0, 1, 6, 8, 15, 21, 27, 28, 30, 33, 35, 38, 41, 45, 47, 51, 53, 55, 58, 191, 194, 197, 200, 205, 213, 215, 219, 221, 226, 229, 232, 234, 237, 240, 241, 245, 251, 254, 262, 265, 313, 314, 315, 317, 320, 321, 322, 372, 373, 380, 389, 391, 394, 397, 419, 420, 421, 429, 430, 431, 434, 448, 454, 476, 479, 486, 496, 507, 509, 517, 520, 532, 542, 563, 567 };
		values = get(p1);
	}

	if (namefile == "calibrar_anillo_nuevo_640x360.wmv"){
		int p1[75] = { 2, 7, 13, 19, 32, 36, 49, 76, 95, 107, 129, 130, 140, 165, 172, 181, 203, 227, 247, 261, 282, 299, 301, 302, 312, 320, 344, 347, 353, 371, 375, 385, 407, 429, 442, 456, 466, 485, 497, 513, 525, 541, 551, 562, 574, 585, 596, 618, 619, 636, 668, 670, 688, 733, 783, 791, 815, 835, 861, 868, 900, 994, 1002, 1016, 1114, 1118, 1211, 1330, 1368, 1382, 1414, 1461, 1481, 1488, 1491 };
		values = get(p1);
	}
	
	if (namefile == "medir_1280x720_anillos.wmv"){
		int p1[75] = { 0, 1, 3, 4, 7, 10, 11, 14, 16, 18, 19, 20, 22, 24, 25, 26, 30, 31, 32, 65, 66, 70, 74, 76, 77, 108, 109, 110, 122, 125, 126, 127, 128, 133, 137, 140, 142, 143, 148, 153, 154, 156, 160, 166, 168, 170, 173, 175, 178, 179, 186, 189, 190, 191, 193, 196, 199, 201, 206, 246, 248, 250, 251, 253, 257, 259, 260, 261, 263, 264, 266, 271, 275, 276, 278};
		values = get(p1);
	}

	if (namefile == "medir_640x360_anillos.wmv"){
		int p1[75] = { 1, 3, 16, 18, 23, 28, 35, 36, 40, 43, 45, 48, 56, 59, 65, 72, 86, 93, 95, 100, 103, 106, 111, 116, 118, 119, 124, 127, 129, 130, 134, 140, 143, 144, 149, 156, 191, 196, 197, 209, 210, 239, 241, 256, 258, 260, 269, 274, 281, 282, 290, 300, 304, 314, 318, 320, 340, 352, 365, 367, 385, 408, 426, 450, 452, 462, 477, 481, 489, 492, 497, 504, 511, 518, 531};
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

	if (nombrevideo == "medir_1280x720_circulos.wmv"){
		cameraMatrix = (Mat_<double>(3, 3) << 995.2879103202191, 0, 652.366423379633,
		0, 990.9977657087171, 342.6447609965957,
		0, 0, 1);
		distCoeffs = (Mat_<double>(8, 1) << 0.01894022119472604,
			-0.00378900400831297,
			0,
			0,
			-0.05370745733157262, 0, 0, 0);
	}

	//cout << "camera Matrix: " << cameraMatrix<< endl;
	//cout << "distCoeffs: " << distCoeffs<< endl;

	if (solvePnP(corners, pointBuf, cameraMatrix, distCoeffs, rvec, tvec)){
		Mat dst = Mat::zeros(3, 3, CV_64F);
		Rodrigues(rvec, dst, noArray());
		double daux = norm((-dst).t()*tvec);
		dist.push_back(daux);
		cout << "distancia " << daux << endl;
		string valuedist = "";
		valuedist = fdist(daux);
		valuedist = "distancia " + valuedist;

		putText(original, valuedist, Point2f(15, 20), FONT_HERSHEY_PLAIN, 1, Scalar(0, 0, 255, 255));
		putText(original, fdist(numframe), Point2f(50,50), FONT_HERSHEY_PLAIN, 1, Scalar(0, 0, 255, 255));
	}
}


int main(){
	// colocar PadronAnillos_01.avi , PadronAnillos_02.avi , PadronAnillos_03.avi,calibrar_anillo_nuevo_1280x720.wmv,calibrar_anillo_nuevo_640x360.wmv,medir_1280x720_anillos,medir_640x360_anillos.wmv que son los nombres de los videos
    nombrevideo = "medir_640x360_anillos.wmv";
	int cantframes = 75;

	vector<int>frames = takeframes(nombrevideo, cantframes);// los frames a considerar para la calibracion
	//vector<int>frames = primeroscantframes(cantframes);
	cout << "cantidad de frames a procesar" << cantframes << " " << frames.size() << endl;
	VideoCapture inputCapture(nombrevideo);
	numframe = 0; // indica el numero de frame que esta siendo procesado en el video.
	int posframe = 0; // sirve para posicionar en el arreglo
	setNumThreads(8);
	/*    ---------------------------------------------------- 1ERA PARTE  ---------------------------------------------*/
	//1ERA PARTE del trabajo corresponde en obtener el tablero de anillos (total 30) sin ruidos
	vector<vector<Point2f> > imagePoints;
	Mat cameraMatrix, distCoeffs;
	ofstream myfile;
	myfile.open("vector.txt");
	vector<pair<double, int> >errorframe;

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
			adaptiveThreshold(src, src, 255.0, CV_THRESH_BINARY, CV_THRESH_BINARY, 101, -5);
			imshow("Adaptative", src);
		}
		Canny(src, src, thresh, 2 * thresh);
		imshow("Canny", src);

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
		double lo = 1; double hi = 60;
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
		cout << "loooo " << lo <<" "<<v.size()<< endl;

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


			//cout << "numFrame: " << numframe << " " << posframe << " " << frames.size() << endl;
			if (nombrevideo == "medir_1280x720_anillos.wmv" || nombrevideo == "medir_640x360_anillos.wmv"){
				distancia(original, pointBuf2);
				numframe++;
			}
			else{
				/**/
				if (numframe == frames[posframe]){
					imagePoints.push_back(pointBuf2);
					posframe++;
				}

				if (posframe == frames.size()){
					cout << "inicia calibracion" << endl;
					runCalibrationAndSave(original, imageSize, cameraMatrix, distCoeffs, imagePoints);
					break;
				}
				
				/*
				double error = runCalibrationAndSave(original, imageSize, cameraMatrix, distCoeffs, vector<vector<Point2f> >(1, pointBuf2));
				cout << "error " << error << endl;
				errorframe.push_back(make_pair(error, numframe));
				*/
				numframe++;
			}
			//ahora simplemente dibujamos el zigzag con el orden ya definido previamente
			drawChessboardCorners(original, boardSize, Mat(pointBuf2), 1);
			namedWindow("Original", WINDOW_AUTOSIZE);
			imshow("Original", original);
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
