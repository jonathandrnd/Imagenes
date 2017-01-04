// PatternGrid.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include <iostream>
#include <omp.h>
#include <queue>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;
Mat src, grey;
int thresh = 100;
bool visited[1001];

// estructura que almacena la informacion de contorno
struct nodo{
	int id;// id del contorno
	int visit;// si es 1 entonces es un contorno fallido y se debe borrar con el metodo limpiar(), 0 contorno valido
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
	if (a.visit != b.visit)return a.visit < b.visit;//prioridad los .visit=0 es decir los validos
	return a.id < b.id;
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
int area(pair<int, int>a, pair<int, int>b, pair<int, int>c){
	int x1 = a.first, y1 = a.second;
	int x2 = b.first, y2 = b.second;
	int x3 = c.first, y3 = c.second;
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




int main(){
	// colocar PadronAnillos_01.avi , PadronAnillos_02.avi , PadronAnillos_03.avi que son los nombres de los videos
	VideoCapture inputCapture("PadronAnillos_03.avi");
	setNumThreads(8);
	/*    ---------------------------------------------------- 1ERA PARTE  ---------------------------------------------*/
	//1ERA PARTE del trabajo corresponde en obtener el tablero de anillos (total 30) sin ruidos


	while (true){
		Mat original;
		// Mat original contiene los frames del video a colores
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
		findContours(src, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

		/// Draw contours
		Mat drawing = Mat::zeros(src.size(), CV_8U);
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
		set<nodo>S;
		// ahora trabajaremos con los centros de cada contorno
		// tenemos circulos concentricos si existen 2 puntos iguales o muy cercanos  solo consideraremos 1 de ellos
		// Set nos permite eliminar duplicados
		for (int i = 0; i < v.size(); i++){
			for (int j = i + 1; j < v.size(); j++)
				if (hypot((v[i].cx - v[j].cx), (v[i].cy - v[j].cy)) <= 1.5){
					S.insert(v[i]);// luego elimino los que tengan el valor 1
				}
		}

		v = vector<nodo>(S.begin(), S.end());
		//luego de eliminar duplicados con Set , si hay 2 centros cercanos solo se considerara 1 de ellos
		for (int i = 0; i < v.size(); i++)
			for (int j = i + 1; j < v.size(); j++)
				if (hypot((v[i].cx - v[j].cx), (v[i].cy - v[j].cy)) <= 1.5){
					v[j].visit = 1;
				}

		//eliminaremos todos los elementos que tengan el atributo visit diferente de 0
		v = limpiar(v);

		// Tenemos los centros de varios contornos y debemos quedarnos con los 30 del patron
		// La observacion esta en que la distancia entre cada centro del patron es pequeña
		// Por lo tanto dado una distancia debemos encontrar el tamaño maximo de las componentes conexas formadas (este debe ser 30)
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

		// en caso de que se haya filtrado los 30 elementos correctamente hacer lo siguiente
		if (v.size() == 30){
			vector<vector<nodo> >vsegmento;
			// vamos a obtener 5 segmentos los cuales contienen 6 puntos cada uno
			for (int caso = 0; caso < 5; caso++){
				double minarea = 1e+10;
				vector<pair<int, int> >puntos6;

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
					pair<int, int>vectora = make_pair(vsegmento[i][5].cx - vsegmento[i][0].cx, vsegmento[i][5].cy - vsegmento[i][0].cy);
					pair<int, int>vectorb = make_pair(vsegmento[j][5].cx - vsegmento[i][0].cx, vsegmento[j][5].cy - vsegmento[i][0].cy);

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

			//ahora simplemente dibujamos el zigzag con el orden ya definido previamente
			drawChessboardCorners(original, boardSize, Mat(pointBuf2), 1);
			namedWindow("Original", WINDOW_AUTOSIZE);
			imshow("Original", original);
		}

		// dibujo del tablero sin ruido
		waitKey(1);
	}
	return 0;
}
