// PatternGrid.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include <iostream>
#include <queue>
#include <cmath>
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

struct nodo{
	int id;// id del contorno
	int visit;// si es 1 entonces es un contorno fallido, 0 contorno valido
	int cx;// centro x
	int cy;// centro y
	nodo(int _id, int _visit, int _cx, int _cy){
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

// area del triangulo
int area(pair<int, int>a, pair<int, int>b, pair<int, int>c){
	int x1 = a.first, y1 = a.second;
	int x2 = b.first, y2 = b.second;
	int x3 = c.first, y3 = c.second;

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


int main(){
	// colocar PadronCirculos_01.avi , PadronCirculos_02.avi , PadronCirculos_03.avi que son los nombres de los videos
	VideoCapture inputCapture("PadronCirculos_02.avi");
	setNumThreads(8);
	/*    ---------------------------------------------------- 1ERA PARTE  ---------------------------------------------*/
	//1ERA PARTE del trabajo corresponde en obtener el tablero de anillos (total 30) sin ruidos

	while (true){
		Mat original;
		if (!inputCapture.read(original))break;

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
		// vector v contendra informacion de los centros de los contornos validos
		v = vector<nodo>(contours.size(), nodo(1, 1, 1, 1));

		//paralelizamos esta parte ya que la funcion fitellipse es lenta
		for (int ii = 0; ii < (int)contours.size(); ii++){
			vector<Point>P = contours[ii];
			Mat pointsf;
			RotatedRect box;
			// Un contorno es valido si el area es positiva y no sea tan grande
			if (contourArea(contours[ii])>1 && contourArea(contours[ii]) < 10000){
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
					v[ii] = nodo(ii, 0, centro.x, centro.y);//valido
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
		Size boardSize(4, 11);
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
		if (nododir[0].cx < nododir[1].cx){
			arcotan = atan((nododir[1].cy - nododir[0].cy) / (0.0 + nododir[1].cx - nododir[0].cx));
		}
		else{
			arcotan = atan((nododir[0].cy - nododir[1].cy) / (0.0 + nododir[0].cx - nododir[1].cx));
		}


		vector<vector<nodo> >vnodo;
		for (int caso = 0; caso < 8; caso++){
			int maxipoints = 0;
			vector<pair<int, int> >puntos4;
			for (int i = 0; i < v.size(); i++){
				if (v[i].visit != 0)continue;
				for (int j = 0; j < v.size(); j++){
					if (v[j].visit != 0)continue;
					// si la recta tiene una pendiente cercana a lo buscado
					if (i != j && v[i].cx < v[j].cx &&  abs(atan((v[j].cy - v[i].cy + 0.0) / (v[j].cx - v[i].cx)) - arcotan) < 0.05){
						int cant = 2;
						vector<pair<int, int> >dis;
						for (int k = 0; k < v.size(); k++){
							if (v[k].visit != 0)continue;
							if (k == i || k == j)continue;
							if (min(v[i].cx, v[j].cx) <= v[k].cx &&  v[k].cx <= max(v[i].cx, v[j].cx) &&
								min(v[i].cy, v[j].cy) <= v[k].cy &&  v[k].cy <= max(v[i].cy, v[j].cy)){
								double aux = area(make_pair(v[i].cx, v[i].cy), make_pair(v[j].cx, v[j].cy), make_pair(v[k].cx, v[k].cy)) / hypot(v[i].cx - v[j].cx, v[i].cy - v[j].cy);
								// si los puntos son colineales
								if (aux <= 4)
									dis.push_back(make_pair(v[k].cx, k));
							}
						}

						cant += dis.size();
						if (maxipoints < cant){
							maxipoints = cant;
							puntos4.clear();
							puntos4.push_back(make_pair(v[i].cx, i));
							puntos4.push_back(make_pair(v[j].cx, j));
							for (int k = 0; k < dis.size(); k++)
								puntos4.push_back(dis[k]);
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

		// ordenamos los segmentos de abajo hacia arriba mediante el producto vectorial
		for (int i = 0; i < vnodo.size(); i++){
			for (int j = i + 1; j < vnodo.size(); j++){
				pair<int, int>vectora = make_pair(vnodo[i][vnodo[i].size() - 1].cx - vnodo[i][0].cx, vnodo[i][vnodo[i].size() - 1].cy - vnodo[i][0].cy);
				pair<int, int>vectorb = make_pair(vnodo[j][vnodo[j].size() - 1].cx - vnodo[i][0].cx, vnodo[j][vnodo[j].size() - 1].cy - vnodo[i][0].cy);
				if (vectora.first*vectorb.second - vectora.second*vectorb.first <= 0){
					swap(vnodo[i], vnodo[j]);
				}
			}
		}

		vector<nodo>auxv;

		if (vnodo.size() == 8)
			for (int i = 0; i < 8; i++){
				for (int j = 0; j < vnodo[i].size(); j++){
					pointBuf.push_back(Point2f(vnodo[i][j].cx, vnodo[i][j].cy));
				}
			}

		/*   // ordenamiento por X o Y --- otro metodo de ordenamiento
		pointBuf = originalPoints;
		for (int i = 0; i < pointBuf.size(); i++)
			for (int j = i + 1; j < pointBuf.size(); j++)
				if (pointBuf[i].x>pointBuf[j].x){
					swap(pointBuf[i].x, pointBuf[j].x);
					swap(pointBuf[i].y, pointBuf[j].y);
				}
				else{
					if (pointBuf[i].x == pointBuf[j].x && pointBuf[i].y > pointBuf[j].y)
						swap(pointBuf[i].y, pointBuf[j].y);
				}
		*/
	
		if (pointBuf.size() == 44){
			drawChessboardCorners(original, boardSize, Mat(pointBuf), 1);
		}

		namedWindow("Original", WINDOW_AUTOSIZE);
		imshow("Original", original);

		waitKey(1);
		//if (cont==1)break;
	}
	return 0;
}

