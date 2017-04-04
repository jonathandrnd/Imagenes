// PatternGrid.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include <iostream>
#include <sstream>
#include <omp.h>
#include <queue>
#include <time.h>
#include <fstream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
//#include <GL/glew.h>
#include "freeglut/include/GL/glut.h"
#include <GL/gl.h>
//#include "glm.h"
#include <time.h>
#include <fstream>
using namespace std;
using namespace cv;
ofstream myfile;
ifstream myinfile;
int readlines = 0;
Mat src, grey;
Mat viewgl;
int thresh = 100;
bool visited[301];
int menor12 = 0;
Size imageSize;
Size boardSize(3, 4);
double tetha = 0.0;
double squareSize = 67.5;
string nombrevideo;
int numframe = 0;
vector<double>dist;
bool getthreshold = 0;
bool lectura = 0;
int doce = 0;
Mat H;
Mat Hinv;
Mat finalrectify;
VideoCapture inputCapture;
vector<Point2f>imgpoints;
vector<int>numframes;
int posnumframes = 0;
bool calibra = 0;
vector<vector<Point2f> >imgp;
vector<vector<Point3f> >objpoints;
double fovx, fovy, fx, fy, px, py, aspect, k1, k2, p1, p2, k3;
// estructura que almacena la informacion de contorno
struct nodo {
	int id;// id del contorno
	int visit;// si es 1 entonces es un contorno fallido y se debe borrar con el metodo limpiar(), 0 contorno valido
	float cx;// centro x
	float cy;// centro y
	float area;
	nodo(int _id, int _visit, float _cx, float _cy, float _area) {
		id = _id; visit = _visit;
		cx = _cx; cy = _cy;
		area = _area;
	}
};
// define el metodo de ordenamiento sort en un nodo
bool operator<(nodo a, nodo b) {
	if (a.area != b.area)return a.area > b.area;
	if (a.cx != b.cx)return a.cx < b.cx;//los menores centros en x
	if (a.cy != b.cy)return a.cy < b.cy;//sino los menores centros en y
	if (a.visit != b.visit)return a.visit < b.visit;//prioridad los .visit=0 es decir los validos
	return a.id < b.id;
}
// consideramos solo los nodos considerados positivos
vector<nodo>limpiar(vector<nodo>v) {
	vector<nodo>ans;
	for (int i = 0; i < v.size(); i++)
		if (v[i].visit == 0)
			ans.push_back(v[i]);
	return ans;
}
vector<nodo>v;
string f(string s) {
	if (s == "calibrar_anillo_nuevo_640x360.wmv")return "salida1.txt";
	if (s == "calibrar_anillo_nuevo_1280x720.wmv")return "salida2.txt";
	if (s == "PadronAnillos_03.avi")return "salida3.txt";
}
//hallar el area de un triangulo
double area(pair<double, double>a, pair<double, double>b, pair<double, double>c) {
	double x1 = a.first, y1 = a.second;
	double x2 = b.first, y2 = b.second;
	double x3 = c.first, y3 = c.second;
	return abs(x2*y3 + x1*y2 + y1*x3 - x2*y1 - x3*y2 - y3*x1);
}
int bfs(int pos, double dis, vector<nodo>&aux) {
	visited[pos] = 1;
	queue<int>Q;
	Q.push(pos);// metemos a la cola el nodo actual
	int cont = 1;
	aux.push_back(v[pos]);
	while (!Q.empty()) {
		int id = Q.front();
		Q.pop();
		for (int i = 0; i < v.size(); i++) {
			// si hay algun nodo no visitado que se encuentra a una distancia dada entonces
			// lo agregamos a la cola
			if (!visited[i] && (v[i].cx - v[id].cx)*(v[i].cx - v[id].cx) + (v[i].cy - v[id].cy)*(v[i].cy - v[id].cy) <= dis*dis) {
				aux.push_back(v[i]);
				cont++;
				Q.push(i);
				visited[i] = 1;
			}
		}
	}
	return cont;
}
vector<int> get(int c[]) {
	vector<int>dev;
	for (int i = 0; i < 65; i++)dev.push_back(c[i]);
	sort(dev.begin(), dev.end());
	return dev;
}
vector<int>takeframes(string namefile, int cant) {
	vector<int>values;
	if (namefile == "Prueba.wmv") {
		int p1[65] = { 92, 95, 122, 153, 160, 201, 207, 277, 283, 338, 341, 350, 360, 375, 397, 398, 438, 449, 470, 474, 493, 500, 503, 522, 573, 575, 607, 613, 629, 632, 661, 669, 677, 682, 694, 700, 710, 734, 737, 860, 881, 916, 925, 944, 953, 969, 971, 1262, 1276, 1290, 1299, 1305, 1354, 1361, 1370, 1383, 1390, 1397, 1402, 1411, 1506, 1507, 1645, 1656, 1661 };
		//int p1[75] = { 92, 95, 122, 153, 160, 201, 207, 277, 283, 338, 341, 350, 360, 371, 375, 389, 397, 398, 438, 449, 465, 470, 474, 483, 493, 500, 503, 522, 541, 544, 573, 575, 607, 613, 629, 632, 661, 669, 677, 682, 694, 700, 710, 734, 737, 751, 767, 771, 774, 860, 881, 916, 925, 944, 953, 969, 971, 1262, 1276, 1290, 1299, 1305, 1354, 1361, 1370, 1383, 1390, 1397, 1402, 1411, 1506, 1507, 1645, 1656, 1661 };
		values = get(p1);
	}
	cout << "values " << values.size() << endl;
	return values;
}
void calcBoardCornerPositions(Size boardSize, float squareSize, vector<Point3f>& corners) {
	corners.clear();
	for (int i = 0; i < boardSize.width; ++i)
		for (int j = 0; j < boardSize.height; ++j)
			corners.push_back(Point3f(float(j*squareSize), float(i*squareSize), 0));
}
string fdist(double d) {
	stringstream st;
	st << d;
	return st.str();
}
bool teapot = true;
bool sphere = false;
// a useful function for displaying your coordinate system
void drawAxes(float length) {
	glPushAttrib(GL_POLYGON_BIT | GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT);
	//glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3f(1, 0, 0);
	glVertex3f(0, 0, 0);
	glVertex3f(length, 0, 0);
	glColor3f(0, 1, 0);
	glVertex3f(0, 0, 0);
	glVertex3f(0, length, 0);
	glColor3f(0, 0, 1);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, length);
	glEnd();
	glPopAttrib();
}
GLuint texture[2];
const GLfloat light_ambient[] = { 1.0f, 0.0f, 0.0f, 1.0f };
const GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 3.0f, 5.0f, 5.0f, 0.0f };
const GLfloat light_ambient2[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_ambient3[] = { 1.0f, 1.0f, 0.0f, 1.0f };
const GLfloat light_diffuse2[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular2[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position2[] = { 3.0f, 5.0f, 5.0f, 0.0f };
const GLfloat mat_ambient[] = { 0.2f, 0.2f, 0.2f, 0.5f };
const GLfloat mat_diffuse[] = { 0.5f, 0.5f, 0.5f, 0.5f };
const GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat high_shininess[] = { 100.0f };
void drawtea(int col, int row) {

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glClear(GL_DEPTH_BUFFER_BIT);

	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
	//glPushAttrib(GL_POLYGON_BIT | GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT);
	for (int i = 0; i < boardSize.width; ++i)
		for (int j = 0; j < boardSize.height; ++j) {
			if (i == col && j == row)continue;
			glPushMatrix();
			glTranslatef(j*squareSize, i*squareSize, 0.0);
			if (j % 2 == 0) {

				glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
				glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
				glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
				glLightfv(GL_LIGHT0, GL_POSITION, light_position);
				glutSolidSphere(squareSize / 4, 20, 200);

				glEnable(GL_BLEND);
				glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient2);
				glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse2);
				glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular2);
				glLightfv(GL_LIGHT0, GL_POSITION, light_position2);
				glutSolidCube(squareSize / 3 * 2.1);
				glDisable(GL_BLEND);

			}
			if (j % 2 == 1) {
				glRotatef(90, 1, 0, 0);
				glRotated(180, 0, 1, 0);
				glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient3);
				glutSolidTeapot(squareSize / 3);
			}
			glPopMatrix();
		}

	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_NORMALIZE);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
	glFlush();

}
void drawsphere(int col, int row) {
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glClear(GL_DEPTH_BUFFER_BIT);

	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);

	glTranslatef(squareSize*1.5, squareSize, squareSize);

	glScaled(2, 2, 2);

	glPushMatrix();

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	tetha = tetha + 15;
	glRotatef(tetha, 1, 1, 1);
	glutSolidCube(squareSize / 2.0);

	glPopMatrix();

	glEnable(GL_BLEND);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient2);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse2);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular2);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position2);
	glutSolidCube(squareSize);
	glDisable(GL_BLEND);

	//glColor4f(0, 1, 1, 1.0);
	//glutSolidSphere(squareSize / 3 * 1.2, 20, 20);

	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_NORMALIZE);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
	glPopMatrix();
	glFlush();
}
void getextrinsicparameters(Mat &cameraMatrix, Mat &distCoeffs) {
	if (nombrevideo == "Prueba.wmv") {
		cameraMatrix = (Mat_<double>(3, 3) << 1030.661648718432, 0, 661.231015316684,
			0, 1035.250924387226, 351.9960472684449,
			0, 0, 1);
		distCoeffs = (Mat_<double>(8, 1) << 0.003534174932036642,
			0.1974014017386418,
			0,
			0,
			-0.5168381889426796, 0, 0, 0);
	}
}
vector<Point2f> readimg(int &id3,int &idfaltante) {
	string line;
	getline(myinfile, line);
	vector<Point2f> ans;
	readlines++;
	if (line.size() == 0) {
		cout << "not found" << endl;
		return ans;
	}
	float x1, y1;
	istringstream is(line);
	while (is >> x1 >> y1) {
		ans.push_back(Point2f(x1, y1));
		id3 = x1; idfaltante = y1;
	}

	ans.erase(ans.begin() + (ans.size() - 1));
	return ans;
}
double computeReprojectionErrors(const vector<vector<Point3f> >& objectPoints,
	const vector<vector<Point2f> >& imagePoints,
	const vector<Mat>& rvecs, const vector<Mat>& tvecs,
	const Mat& cameraMatrix, const Mat& distCoeffs,
	vector<float>& perViewErrors) {
	vector<Point2f> imagePoints2;
	int i, totalPoints = 0;
	double totalErr = 0, err;
	perViewErrors.resize(objectPoints.size());
	for (i = 0; i < (int)objectPoints.size(); ++i) {
		projectPoints(Mat(objectPoints[i]), rvecs[i], tvecs[i], cameraMatrix,
			distCoeffs, imagePoints2);
		err = norm(Mat(imagePoints[i]), Mat(imagePoints2), CV_L2);
		int n = (int)objectPoints[i].size();
		perViewErrors[i] = (float)std::sqrt(err*err / n);
		cout << perViewErrors[i] << endl;
		totalErr += err*err;
		totalPoints += n;
	}
	return std::sqrt(totalErr / totalPoints);
}
bool runCalibration(Size& imageSize, Mat& cameraMatrix, Mat& distCoeffs,
	vector<vector<Point2f> > imagePoints, vector<Mat>& rvecs, vector<Mat>& tvecs, vector<vector<Point3f> > objectPoints) {
	cameraMatrix = Mat::eye(3, 3, CV_64F);
	distCoeffs = Mat::zeros(8, 1, CV_64F);
	//Find intrinsic and extrinsic camera parameters
	double rms = calibrateCamera(objectPoints, imagePoints, imageSize, cameraMatrix,
		distCoeffs, rvecs, tvecs, CV_CALIB_ZERO_TANGENT_DIST);
	vector<float>reprojErrs; double totalAvgErr;
	cout << "Re-projection error reported by calibrateCamera: " << rms << endl;
	totalAvgErr = computeReprojectionErrors(objectPoints, imagePoints,
		rvecs, tvecs, cameraMatrix, distCoeffs, reprojErrs);
	cout << cameraMatrix << endl;
	cout << distCoeffs << endl;
	return 1;
}
double runCalibrationAndSave(Size imageSize, Mat& cameraMatrix, Mat& distCoeffs, vector<vector<Point2f> > imagePoints, vector<vector<Point3f> > objectPoints) {
	vector<Mat> rvecs, tvecs;
	bool ok = runCalibration(imageSize, cameraMatrix, distCoeffs, imagePoints, rvecs, tvecs,
		objectPoints);
	return 0;
}
void display() {
	/* ---------------------------------------------------- 1ERA PARTE ---------------------------------------------*/
	//1ERA PARTE del trabajo corresponde en obtener el tablero de anillos (total 30) sin ruidos
	//vector<vector<Point2f> > imagePoints;
	//vector<pair<double, int> >errorframe;
	numframe++;
	int contx = 0;
	Mat cameraMatrix = Mat::eye(3, 3, CV_64F);
	Mat distCoeffs = Mat::zeros(8, 1, CV_64F);
	//Mat original = viewgl.clone();
	// Mat original contiene los frames del video a colores
	imageSize = viewgl.size();
	cout << numframe << " " << readlines << endl;
	//imshow("viewgl ", viewgl);
	int idfaltante = -1;
	int id3 = -1;
	vector<Point2f>puntos;
	if (!lectura) {
		//convertimos a escala de grises
		cv::cvtColor(viewgl, src, CV_BGR2GRAY);
		//eliminamos el ruido con desenfoque gaussiano
		GaussianBlur(src, src, Size(5, 5), 17, 17);
		//Obtenemos los bordes de la imagen con Canny
		if (imageSize.width > 700) {
			if (nombrevideo != "PadronAnillos_03.avi") {
				getthreshold = 1;
				//adaptiveThreshold(src, src, 255.0, CV_THRESH_BINARY, CV_THRESH_BINARY, 101, -5);
				//imshow("Adaptative", src);
			}
		}
		cv::Canny(src, src, 50, 2 * 50);
		//imshow("Canny", src);
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;
		//Encontramos todos los contornos de la imagen y lo guardamos en el vector Contours
		cv::findContours(src, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
		/// Draw contours
		//Mat drawing2 = Mat::zeros(src.size(), CV_8UC3);
		Mat drawing = Mat::zeros(src.size(), CV_8U);
		vector<Point2f>xxx;
		// vector v contendra informacion de los centros de los contornos validos
		v = vector<nodo>(contours.size(), nodo(1, 1, 1, 1, 1));
		for (int ii = 0; ii < (int)contours.size(); ii++) {
			vector<Point>P = contours[ii];
			Mat pointsf;
			RotatedRect box;
			// Un contorno es valido si el area es positiva y no sea tan grande
			if (contourArea(contours[ii])>150 && contourArea(contours[ii]) < 10000 && contours[ii].size() > 5) {
				Mat pointsf;
				RotatedRect box;
				//utilizamos fitellipse para obtener un rectangulo que lo contenga
				try {
					Mat(contours[ii]).convertTo(pointsf, CV_32F);
					box = fitEllipse(pointsf);
				}
				catch (cv::Exception& e) {
					// en caso se produzca excepcion, si no se pueda formar la elipse con los puntos del contorno
					continue;
				}
				float f1 = box.size.height;
				float f2 = box.size.width;
				// si el largo y el ancho difieren minimamente en tamaño y el area no es tan grande sera considerado
				// como posible contorno
				if (fabs(f2 - f1) <= 35 && f1*f2 <= 8000) {
					Point2f centro = box.center;
					v[ii] = nodo(ii, 0, centro.x, centro.y, f1*f2);
				}
			}
		}
		//sort(prueba.begin(), prueba.end());
		//for (int i = 0; i < prueba.size(); i++)
		// cout << prueba[i].first.first << " " << prueba[i].first.second << " " << prueba[i].second << endl;
		v = limpiar(v);
		//bool found = findCirclesGrid(drawing, Size(3,3), xxx, CALIB_CB_ASYMMETRIC_GRID);
		//cout << "found " << found <<" "<<xxx.size()<<endl;
		std::sort(v.begin(), v.end());
		set<nodo>S;
		cout << "inicia v.size() " << v.size() << endl;
		// ahora trabajaremos con los centros de cada contorno
		// tenemos circulos concentricos si existen 2 puntos iguales o muy cercanos solo consideraremos 1 de ellos
		// Set nos permite eliminar duplicados
		std::memset(visited, 0, sizeof(visited));
		for (int i = 0; i < v.size(); i++) {
			if (visited[i])continue;
			double sumx = 0; double sumy = 0;
			int cont = 0;
			for (int j = 0; j < v.size(); j++) {
				if (i == j)continue;
				if (visited[j])continue;
				if (abs((v[i].cx - v[j].cx)) <= 3.5 && abs(v[i].cy - v[j].cy) <= 3.5 && v[i].area>v[j].area) {
					visited[j] = 1;
					sumx += v[i].cx;
					sumy += v[i].cy;
					cont++;
				}
				//if (v[j].cx > v[i].cx + 10)break;
			}
			if (cont >= 1) {
				v[i].cx = sumx / cont;
				v[i].cy = sumy / cont;
				S.insert(v[i]);
			}
		}
		v = vector<nodo>(S.begin(), S.end());
		//luego de eliminar duplicados con Set , si hay 2 centros cercanos solo se considerara 1 de ellos
		for (int i = 0; i < v.size(); i++)
			for (int j = 0; j < v.size(); j++)
				if (i != j && abs((v[i].cx - v[j].cx)) <= 3.5 && abs(v[i].cy - v[j].cy) <= 3.5 && v[i].area>v[j].area) {
					v[j].visit = 1;
				}
		//eliminaremos todos los elementos que tengan el atributo visit diferente de 0
		v = limpiar(v);
		/*
		cout << "yohohohoho " << v.size() << endl;
		for (int i = 0; i < v.size(); i++)
		cout << v[i].cx << " " << v[i].cy << " " << v[i].area << endl;
		*/
		if (v.size() > 12) {
			vector<double>distancia(v.size());
			for (int i = 0; i < v.size(); i++) {
				double minimo1 = 1e+10;
				double minimo2 = 1e+9;
				for (int j = 0; j < v.size(); j++) {
					if (i == j)continue;
					double disx = (double)hypot(v[i].cx - v[j].cx, v[i].cy - v[j].cy);
					if (disx < minimo2) {
						minimo1 = minimo2;
						minimo2 = disx;
					}
					else {
						if (disx < minimo1) {
							minimo1 = disx;
						}
					}
				}
				//cout <<"minimo "<< minimo2 << endl;
				distancia[i] = minimo2;
			}
			vector<double>distancia2 = distancia;
			sort(distancia2.begin(), distancia2.end());
			double dx = distancia2[distancia2.size() / 2];
			for (int i = 0; i < v.size(); i++) {
				if (distancia[i] > dx*1.3)
					v[i].visit = 1;
			}
			v = limpiar(v);
		}
		cout << "v.size() " << v.size() << endl;
		for (int i = 0; i < v.size(); i++)
			drawContours(drawing, contours, v[i].id, Scalar(255, 255, 255), CV_FILLED);
		//imshow("drawing", drawing);
		//El tablero es de 5X6
		Size boardSize(3, 4);
		// en caso de que se haya filtrado los 30 elementos correctamente hacer lo siguiente
		if (v.size() < 12)
			menor12++;
		if (v.size() == 12) {
			doce++;
			cout << "doce " << doce << " " << menor12 << endl;
			vector<Point3f> corners;
			for (int i = 0; i < boardSize.width; ++i)
				for (int j = 0; j < boardSize.height; ++j)
					corners.push_back(Point3f(float(j*squareSize), float(i*squareSize), 0));
			//for (int i = 0; i < v.size(); i++)
			// puntos.push_back(Point2f((float)v[i].cx, (float)v[i].cy));
			vector<vector<nodo> >vsegmento;
			for (int caso = 0; caso < 3; caso++) {
				double minarea = 1e+10;
				vector<pair<double, int> >puntos4;
				for (int i = 0; i < v.size(); i++) {
					if (v[i].visit != 0)continue;
					for (int j = 0; j < v.size(); j++) {
						if (v[j].visit != 0)continue;
						if (i != j && v[i].cx <= v[j].cx) {
							vector<pair<double, int> >dis;
							// fijamos los puntos i , j los cuales vienen hacer los extremos del segmento que debe incluir 4 puntos
							for (int k = 0; k < v.size(); k++) {
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
							double area2 = dis[0].first + dis[1].first;
							if (minarea > area2) {
								minarea = area2;
								puntos4.clear();
								puntos4.push_back(make_pair(v[i].cx, i));
								puntos4.push_back(make_pair(v[dis[0].second].cx, dis[0].second));
								puntos4.push_back(make_pair(v[dis[1].second].cx, dis[1].second));
								puntos4.push_back(make_pair(v[j].cx, j));
							}
						}
					}
				}
				vector<nodo>aux;
				// ordenamiento por x
				sort(puntos4.begin(), puntos4.end());
				// luego de obtener la recta con 4 puntos lo marcamos como visitados
				// Para que no obtener la misma recta en la siguiente iteracion
				for (int i = 0; i < puntos4.size(); i++) {
					v[puntos4[i].second].visit = caso + 1;
					aux.push_back(v[puntos4[i].second]);
				}
				vsegmento.push_back(aux);
			}
			for (int i = 0; i < vsegmento.size(); i++) {
				for (int j = i + 1; j < vsegmento.size(); j++) {
					pair<double, double>vectora = make_pair(vsegmento[i][3].cx - vsegmento[i][0].cx, vsegmento[i][3].cy - vsegmento[i][0].cy);
					pair<double, double>vectorb = make_pair(vsegmento[j][3].cx - vsegmento[i][0].cx, vsegmento[j][3].cy - vsegmento[i][0].cy);
					if (vectora.first*vectorb.second - vectora.second*vectorb.first < 0) {
						swap(vsegmento[i], vsegmento[j]);
					}
				}
			}
			// el metodo anterior ordena los segmentos de arriba hacia abajo usando ordenamiento burbuja
			for (int i = 0; i < 3; i++) { // 5 segmentos ordenados de abajo hacia arriba
				for (int j = 0; j < vsegmento[2 - i].size(); j++) { // cada segmento de tamanho 6 esta ordenado por x
					puntos.push_back(Point2f((float)vsegmento[2 - i][j].cx, (float)vsegmento[2 - i][j].cy));
					//pointBuf.push_back(Point2f(vsegmento[i][j].cx, vsegmento[i][j].cy));
				}
			}
		}
		if (v.size() == 11) {
			vector<vector<nodo> >vsegmento;
			for (int caso = 0; caso < 3; caso++) {
				double minarea = 1e+10;
				vector<pair<double, int> >puntos4;
				for (int i = 0; i < v.size(); i++) {
					if (v[i].visit != 0)continue;
					for (int j = 0; j < v.size(); j++) {
						if (v[j].visit != 0)continue;
						if (i != j && v[i].cx <= v[j].cx) {
							vector<pair<double, int> >dis;
							// fijamos los puntos i , j los cuales vienen hacer los extremos del segmento que debe incluir 4 puntos
							for (int k = 0; k < v.size(); k++) {
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
							double area2 = dis[0].first;
							if (caso != 2)area2 = dis[0].first + dis[1].first;
							if (minarea > area2) {
								minarea = area2;
								puntos4.clear();
								if (caso != 2) {
									puntos4.push_back(make_pair(v[i].cx, i));
									puntos4.push_back(make_pair(v[dis[0].second].cx, dis[0].second));
									puntos4.push_back(make_pair(v[dis[1].second].cx, dis[1].second));
									puntos4.push_back(make_pair(v[j].cx, j));
								}
								else {
									puntos4.push_back(make_pair(v[i].cx, i));
									puntos4.push_back(make_pair(v[dis[0].second].cx, dis[0].second));
									puntos4.push_back(make_pair(v[j].cx, j));
								}
							}
						}
					}
				}
				vector<nodo>aux;
				// ordenamiento por x
				sort(puntos4.begin(), puntos4.end());
				// luego de obtener la recta con 4 puntos lo marcamos como visitados
				// Para que no obtener la misma recta en la siguiente iteracion
				for (int i = 0; i < puntos4.size(); i++) {
					v[puntos4[i].second].visit = caso + 1;
					aux.push_back(v[puntos4[i].second]);
				}
				vsegmento.push_back(aux);
			}
			for (int i = 0; i < vsegmento.size(); i++) {
				for (int j = i + 1; j < vsegmento.size(); j++) {
					pair<double, double>vectora = make_pair(vsegmento[i][vsegmento[i].size() - 1].cx - vsegmento[i][0].cx, vsegmento[i][vsegmento[i].size() - 1].cy - vsegmento[i][0].cy);
					pair<double, double>vectorb = make_pair(vsegmento[j][vsegmento[j].size() - 1].cx - vsegmento[i][0].cx, vsegmento[j][vsegmento[j].size() - 1].cy - vsegmento[i][0].cy);
					if (vectora.first*vectorb.second - vectora.second*vectorb.first < 0) {
						swap(vsegmento[i], vsegmento[j]);
					}
				}
			}
			for (int i = 0; i < 3; i++)
				if (vsegmento[i].size() == 3)id3 = i;
			if (id3 == -1)return;
			for (int i = 0; i < 3; i++)if (i != id3 && vsegmento[i].size() != 4)return;
			double maxarea = 0;
			for (int j = 0; j < 4; j++) {
				vector<pair<double, double> >x(3, make_pair(0, 0));
				int cont = 0;
				for (int i = 0; i < 3; i++) {
					if (i != id3)
						x[cont++] = make_pair(vsegmento[i][j].cx, vsegmento[i][j].cy);
				}
				double minarea = 1e+10;
				for (int i = 0; i < 3; i++)
					minarea = min(minarea, area(x[0], x[1], make_pair(vsegmento[id3][i].cx, vsegmento[id3][i].cy)));
				if (maxarea < minarea) {
					maxarea = minarea;
					idfaltante = j;
				}
			}
			if (id3 == 0) {
				pair<double, double>pseg = make_pair(2 * vsegmento[1][idfaltante].cx - vsegmento[2][idfaltante].cx,
					2 * vsegmento[1][idfaltante].cy - vsegmento[2][idfaltante].cy);
				vsegmento[id3].push_back(nodo(0, 0, pseg.first, pseg.second, 0));
			}
			if (id3 == 1) {
				pair<double, double>pseg = make_pair((vsegmento[2][idfaltante].cx + vsegmento[0][idfaltante].cx) / 2,
					(vsegmento[2][idfaltante].cy + vsegmento[0][idfaltante].cy) / 2);
				vsegmento[id3].push_back(nodo(0, 0, pseg.first, pseg.second, 0));
			}
			if (id3 == 2) {
				pair<double, double>pseg = make_pair(2 * vsegmento[1][idfaltante].cx - vsegmento[0][idfaltante].cx,
					2 * vsegmento[1][idfaltante].cy - vsegmento[0][idfaltante].cy);
				vsegmento[id3].push_back(nodo(0, 0, pseg.first, pseg.second, 0));
			}
			for (int i = 0; i < vsegmento.size(); i++)
				for (int j = 0; j < 4; j++)
					for (int k = j + 1; k < 4; k++)
						if (vsegmento[i][j].cx >vsegmento[i][k].cx)
							swap(vsegmento[i][j], vsegmento[i][k]);
			// el metodo anterior ordena los segmentos de arriba hacia abajo usando ordenamiento burbuja
			for (int i = 0; i < 3; i++) { // 5 segmentos ordenados de abajo hacia arriba
				for (int j = 0; j < vsegmento[2 - i].size(); j++) { // cada segmento de tamanho 6 esta ordenado por x
					puntos.push_back(Point2f((float)vsegmento[2 - i][j].cx, (float)vsegmento[2 - i][j].cy));
					//pointBuf.push_back(Point2f(vsegmento[i][j].cx, vsegmento[i][j].cy));
				}
			}
		}
	}else{
		puntos = readimg(id3, idfaltante);
		cout << "puntos " << puntos.size() << endl;
	}

	if (puntos.size() != 12)return;
	// clear the window
	glClear(GL_COLOR_BUFFER_BIT);
	// show the current camera frame
	//imshow("display ", viewgl);
	//based on the way cv::Mat stores data, you need to flip it before displaying it
	//Mat tempimage;
	flip(viewgl, viewgl, 0);
	// OPENCV utiliza canales BGR y OpenGL RGB por lo tanto invertimos el canal B y R
	for (int i = 0; i < viewgl.size().width; i++)
		for (int j = 0; j < viewgl.size().height; j++)
			swap(viewgl.ptr<Vec3b>(j)[i][0], viewgl.ptr<Vec3b>(j)[i][2]);
	glDrawPixels(viewgl.size().width, viewgl.size().height, GL_RGB, GL_UNSIGNED_BYTE, viewgl.ptr());
	//////////////////////Read Camera Matrix y Object points/////////////////////////////////
	vector<Point3f> corners;
	calcBoardCornerPositions(boardSize, squareSize, corners);
	getextrinsicparameters(cameraMatrix, distCoeffs);
	////////////////// Projection Matrix ///////////////////////////////////////
	Mat projMat = Mat::zeros(4, 4, CV_64FC1);
	float zfar = 10000.0f, znear = 0.1f;
	projMat.at<double>(0, 0) = 2 * cameraMatrix.at<double>(0, 0) / viewgl.size().width;
	projMat.at<double>(0, 1) = 2 * cameraMatrix.at<double>(0, 1) / viewgl.size().width;
	projMat.at<double>(0, 2) = 1 - (2 * cameraMatrix.at<double>(0, 2) / viewgl.size().width); // en la diapo del profe es su negativo se equivoco
	projMat.at<double>(1, 1) = 2 * cameraMatrix.at<double>(1, 1) / viewgl.size().height;
	projMat.at<double>(1, 2) = -1 + ((2 * cameraMatrix.at<double>(1, 2)) / viewgl.size().height);// en la diapo del profe es su negativo se equivoco
	projMat.at<double>(2, 2) = -(zfar + znear) / (zfar - znear);
	projMat.at<double>(2, 3) = -2 * zfar*znear / (zfar - znear);
	projMat.at<double>(3, 2) = -1;
	double projectionMatrix[16];
	//La matriz projMat sacamos su transpuesta y lo metemos todo a un arreglo
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			projectionMatrix[4 * j + i] = projMat.at<double>(i, j);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, viewgl.size().width, viewgl.size().height);
	glLoadMatrixd(projectionMatrix);
	Mat rvec, tvec;
	vector<Mat> rvecs; vector<Mat> tvecs;
	int chessflags = 0;
	cv::Mat gimage;
	cvtColor(viewgl, gimage, CV_BGR2GRAY);
	Mat rotation;
	Mat viewMatrix = cv::Mat::zeros(4, 4, CV_64FC1);
	solvePnP(corners, puntos, cameraMatrix, distCoeffs, rvec, tvec);
	Rodrigues(rvec, rotation);
	for (int row = 0; row < 3; ++row) {
		for (int col = 0; col < 3; ++col)
			viewMatrix.at<double>(row, col) = rotation.at<double>(row, col);
		viewMatrix.at<double>(row, 3) = tvec.at<double>(row, 0);
	}
	viewMatrix.at<double>(3, 3) = 1.0f;
	//cout<<"viewmatrix"<<viewMatrix<<endl;
	//cameraMatrix*viewMatrix
	Mat auxcameraMatrix = cv::Mat::zeros(3, 3, CV_64FC1);
	auxcameraMatrix.at<double>(0, 0) = cameraMatrix.at<double>(0, 0);
	auxcameraMatrix.at<double>(1, 1) = cameraMatrix.at<double>(1, 1);
	auxcameraMatrix.at<double>(0, 2) = cameraMatrix.at<double>(0, 2);
	auxcameraMatrix.at<double>(1, 2) = cameraMatrix.at<double>(1, 2);
	auxcameraMatrix.at<double>(2, 2) = cameraMatrix.at<double>(2, 2);
	Mat auxview = cv::Mat::zeros(3, 4, CV_64FC1);
	for (int row = 0; row < 3; ++row) {
		for (int col = 0; col < 3; ++col)
			auxview.at<double>(row, col) = viewMatrix.at<double>(row, col);
		auxview.at<double>(row, 3) = viewMatrix.at<double>(row, 3);
	}
	cv::Mat cvToGl = cv::Mat::zeros(4, 4, CV_64F);
	cvToGl.at<double>(0, 0) = 1.0f;
	cvToGl.at<double>(1, 1) = -1.0f;
	cvToGl.at<double>(2, 2) = -1.0f;
	cvToGl.at<double>(3, 3) = 1.0f;
	viewMatrix = cvToGl * viewMatrix;
	Mat glViewMatrix;
	transpose(viewMatrix, glViewMatrix);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glLoadMatrixd(&glViewMatrix.at<double>(0, 0));
	glPushMatrix();
	drawAxes(500.0);
	//glScalef(1.0,-1.0,-1.0);
	//rook();
	if (teapot) {
		drawtea(2 - id3, idfaltante);
	}
	if (sphere) {
		drawsphere(2 - id3, idfaltante);
	}
	glPopMatrix();
	// show the rendering on the screen
	glutSwapBuffers();
	// post the next redisplay
	waitKey(1);
}
void reshape(int w, int h) {
	w = imageSize.width;
	h = imageSize.height;
	// set OpenGL viewport (drawable area)
	glViewport(0, 0, w, h);
}
void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {}
}
void keyboard(unsigned char key, int x, int y) {
	if (key == 't') {
		if (teapot == true) {
			teapot = false;
			sphere = true;
		}
		else {
			teapot = true;
			sphere = false;
		}
	}
	if (key == 's') {
		if (sphere == false) {
			sphere = true;
			teapot = false;
		}
		else {
			sphere = false;
			teapot = true;
		}
	}
}
void idle() {
	// grab a frame from the camera
	(inputCapture) >> viewgl;
	glutPostRedisplay();
}
int main(int argc, char **argv) {
	// colocar PadronAnillos_01.avi , PadronAnillos_02.avi , PadronAnillos_03.avi,calibrar_anillo_nuevo_1280x720.wmv,calibrar_anillo_nuevo_640x360.wmv,medir_1280x720_anillos,medir_640x360_anillos.wmv que son los nombres de los videos
	nombrevideo = "Prueba.wmv";
	int cantframes = 87;
	myinfile.open("salida.txt");
	//myinfile.open("vector.txt");
	numframes = takeframes("Prueba.wmv", 65);
	//vector<int>frames = takeframes(nombrevideo, cantframes);// los frames a considerar para la calibracion
	//int p1[91] = { 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 28, 30, 32, 35, 38, 40, 42, 44, 46, 48, 50, 52, 54, 58, 102, 191, 193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223, 230, 238, 239, 241, 243, 245, 247, 249, 250, 253, 255, 257, 259, 262, 264, 318, 320, 321, 391, 393, 394, 396, 416, 419, 431, 433, 435, 450, 452, 454, 456, 474, 477, 479, 481, 482, 484, 486, 506, 508, 509, 519, 540, 543, 564, 567 };
	vector<int>frames;
	//for (int i = 0; i < cantframes; i++)frames.push_back(p1[i]);
	//vector<int>frames = primeroscantframes(cantframes);
	cout << "cantidad de frames a procesar" << cantframes << " " << frames.size() << endl;
	inputCapture = VideoCapture(nombrevideo);
	numframe = 0; // indica el numero de frame que esta siendo procesado en el video.
	int posframe = 0; // sirve para posicionar en el arreglo
	/* ---------------------------------------------------- 1ERA PARTE ---------------------------------------------*/
	//1ERA PARTE del trabajo corresponde en obtener el tablero de anillos (total 30) sin ruidos
	vector<vector<Point2f> > imagePoints;
	vector<pair<double, int> >errorframe;
	(inputCapture) >> viewgl;
	glutInit(&argc, argv);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(viewgl.cols, viewgl.rows);
	glutCreateWindow("OpenGL");
	Mat undi, rview, map1, map2;
	// set up GUI callback functions
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard);
	glutIdleFunc(idle);
	/**/
	// start GUI loop
	glutMainLoop();
	return 0;
}