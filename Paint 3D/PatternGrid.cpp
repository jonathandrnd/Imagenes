// PatternGrid.cpp : Defines the entry point for the console application.
#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <queue>
#include <time.h>
#include <fstream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <GL/glut.h>
//#include "freeglut-2.8.1/include/GL/glut.h"
#include <GL/gl.h>
#include <SOIL/SOIL.h>
#include <time.h>
#include <fstream>


#define MAXLINES 5000      // Numero maximo de lineas en 3D
using namespace std;
using namespace cv;
ofstream myfile;
ifstream myinfile;
int readlines = 0;int cz=0;
double anglerotatex=0,anglerotatey=0,anglerotatez=0;
int nlinesdraw = 0;
int lastframe=-1<<20;
struct lines {
	double pos[3];
	int color;
	int width;
	bool drawon;
};
GLuint texture[3];
int numframe=0;
double xxxxxx = 0;
lines vlines[MAXLINES + 1];
int color = 1;
int width = 4;
double basecolor[10][3];
double delta=0;
bool drawon = 1;// modo de dibujo activo
double iniciox = -1;
double inicioy = -1;
double inicioz = -1;
double purplex=-1;
double purpley=-1;
double verdex=-1;
double verdey=-1;
double angryx=-1;
double angryy=-1;
double patroninix=-1<<30;double patronx=-1<<30;
double patroniniy=-1<<30;double patrony=-1<<30;
double inidirx=-1<<30;double inidiry=-1<<30;
double dirx=-1<<30;double diry=-1<<30;
Mat src, grey;
Mat viewgl;
int thresh = 100;
bool visited[301];
int framecorrecto = 0;
Size imageSize;
Size boardSize(3, 4);
double tetha = 0.0; double rotatepaint = 0.0;
double squareSize = 30;
double squareSize2 = 40;
string nombrevideo;
int idangry=0;
vector<double>dist;
bool getthreshold = 0;
bool lectura = 0;
int cte = 90;
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

int LoadGLTextures(){
	texture[0]= SOIL_load_OGL_texture(
        "angry.png",
        SOIL_LOAD_AUTO,
        SOIL_CREATE_NEW_ID,
        SOIL_FLAG_INVERT_Y
        );
    if (texture[0] == 0)
        return false;
    // Typical Texture Generation Using Data From The Bitmap
    glBindTexture(GL_TEXTURE_2D, texture[0]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	texture[1]= SOIL_load_OGL_texture(
        "angry2.png",
        SOIL_LOAD_AUTO,
        SOIL_CREATE_NEW_ID,
        SOIL_FLAG_INVERT_Y
        );
    if (texture[1] == 0)
        return false;
    // Typical Texture Generation Using Data From The Bitmap
    glBindTexture(GL_TEXTURE_2D, texture[1]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	texture[2]= SOIL_load_OGL_texture(
        "angry3.png",
        SOIL_LOAD_AUTO,
        SOIL_CREATE_NEW_ID,
        SOIL_FLAG_INVERT_Y
        );
    if (texture[2] == 0)
        return false;
    // Typical Texture Generation Using Data From The Bitmap
    glBindTexture(GL_TEXTURE_2D, texture[2]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);



	return true;
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
	for (int i = 0; i < 75; i++)dev.push_back(c[i]);
	sort(dev.begin(), dev.end());
	return dev;
}

vector<int>takeframes(string namefile, int cant) {
	vector<int>values;
	if (namefile == "video2.mp4") {
		int p1[75] = { 4, 8, 16, 20, 26, 31, 35, 37, 43, 46, 50, 55, 58, 63, 65, 70, 77, 83, 89, 93, 98, 102, 117, 125, 131, 140, 148, 160, 173, 185, 316, 330, 340, 355, 362, 386, 396, 400, 409, 411, 419, 423, 430, 439, 449, 452, 457, 461, 466, 470, 477, 481, 488, 493, 498, 508, 513, 518, 528, 534, 574, 579, 583, 587, 591, 598, 604, 610, 707, 712, 719, 722, 726, 729, 736 };
		//int p1[75] = { 92, 95, 122, 153, 160, 201, 207, 277, 283, 338, 341, 350, 360, 371, 375, 389, 397, 398, 438, 449, 465, 470, 474, 483, 493, 500, 503, 522, 541, 544, 573, 575, 607, 613, 629, 632, 661, 669, 677, 682, 694, 700, 710, 734, 737, 751, 767, 771, 774, 860, 881, 916, 925, 944, 953, 969, 971, 1262, 1276, 1290, 1299, 1305, 1354, 1361, 1370, 1383, 1390, 1397, 1402, 1411, 1506, 1507, 1645, 1656, 1661 };
		values = get(p1);
	}
	cout << "values " << values.size() << endl;
	return values;
}

void readpoints(vector<Point2f>&img30, vector<Point2f>&img12){
	string readline;
	getline(myinfile, readline);
	
	istringstream is(readline);
	if (readline.size()!=0){
		for (int i = 0; i < 30; i++){
			float f1, f2;
			is >> f1 >> f2;
			img30.push_back(Point2f(f1, f2));
		}
	}

	getline(myinfile, readline);

	istringstream is2(readline);
	if (readline.size() != 0){
		for (int i = 0; i < 12; i++){
			float f1, f2;
			is2 >> f1 >> f2;
			img12.push_back(Point2f(f1, f2));
		}
	}
	
	return;
}

void calcBoardCornerPositions(Size boardSize, float squareSize, vector<Point3f>& corners) {
	corners.clear();
	for (int i = 0; i < boardSize.height; ++i)
		for (int j = 0; j < boardSize.width; ++j)
			corners.push_back(Point3f(float(j*squareSize), float(i*squareSize), 0));
}

string fdist(double d) {
	stringstream st;
	st << d;
	return st.str();
}
bool teapot = true;
bool sphere = false;

 GLfloat light_ambient[] = { 1.0f, 0.0f, 0.0f, 1.0f };
const GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 3.0f, 5.0f, 5.0f, 0.0f };
const GLfloat light_ambient2[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_ambient3[] = { 1.0f, 1.0f, 0.0f, 1.0f };
const GLfloat light_diffuse2[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular2[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position2[] = { 3.0f, 5.0f, 5.0f, 0.0f };
GLfloat mat_ambient[] = { 0.2f, 0.2f, 0.2f, 0.5f };
const GLfloat mat_diffuse[] = { 0.5f, 0.5f, 0.5f, 0.5f };
const GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat high_shininess[] = { 100.0f };

void drawPaint(){

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glClear(GL_DEPTH_BUFFER_BIT);

    glPushMatrix();

	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_ambient[0]=1.0f;light_ambient[1]= 1.0f;light_ambient[2]= 1.0f;light_ambient[3]=1.0f;
	

	angryx=940+60*cos(rotatepaint);angryy=viewgl.size().height-570+60*sin(rotatepaint);
	glTranslatef(940, 570,0);
	//mat_ambient[0]=0.0f;mat_ambient[1]= 0.0f;mat_ambient[2]=1.0f;mat_ambient[3]=0.0f;
	//glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glTranslatef(60*cos(rotatepaint),60*sin(rotatepaint),300*sin(rotatepaint));
	glRotatef(270,1.0f,0.0f,0.0f);
	glRotatef(30+cz*2,0.0f,cz+0.0f,0.0f);	
	cz++;
	//glTranslatef(rotatepaint,rotatepaint,rotatepaint);
	//if(rotatepaint>100)rotatepaint-=10;
	//glTranslatef(rotatepaint>10?rotatepaint:0,rotatepaint>10?rotatepaint:0,rotatepaint>10?rotatepaint:0);

	GLUquadricObj *qObj = gluNewQuadric();
	gluQuadricNormals(qObj, GLU_SMOOTH);
	gluQuadricTexture(qObj, GL_TRUE);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture[idangry]);    // texID is the texture ID of a
	gluSphere(qObj, 30.0f, 360, 180);
	glDisable(GL_TEXTURE_2D);
	glPopMatrix();

	//mat_ambient[0]=0.2f;mat_ambient[1]= 0.2f;mat_ambient[2]= 0.2f;mat_ambient[3]=0.5f;
	mat_ambient[0]=1.0f;mat_ambient[1]= 0.0f;mat_ambient[2]= 0.0f;mat_ambient[3]=0.0f;
	light_ambient[0]=1.0f;light_ambient[1]= 0.0f;light_ambient[2]= 0.0f;light_ambient[3]=1.0f;
	

	glPushMatrix();
	purplex=520;purpley=100;
	glTranslatef(purplex, viewgl.size().height-purpley, 800);
	//rotatepaint += 3;
	//glRotatef(rotatepaint, rotatepaint / 1, rotatepaint / 3, rotatepaint/2);
	//glutSolidSphere(30, 20, 20);
	glRotatef(94, 5, 12,12);	
	mat_ambient[3]=0.5f;mat_ambient[0]= 204.0f/ 255;mat_ambient[1]= 0.0f;mat_ambient[2]=102.0f/255;
	light_ambient[0]=1.0f;light_ambient[1]= 1.0f;light_ambient[2]= 1.0f;light_ambient[3]=1.0f;	
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glutSolidCube(35);	
	glPopMatrix();

	glPushMatrix();
	//glRotatef(-rotatepaint,- rotatepaint / 1, -rotatepaint / 3, -rotatepaint / 2);
	verdex=800;verdey=400;	
	glTranslatef(verdex, viewgl.size().height-400,800);
	//rotatepaint += 3;
	mat_ambient[0]=0.1f;mat_ambient[1]= 1.0f;mat_ambient[2]= 0.0f;mat_ambient[3]=0.0f;
	light_ambient[0]=1.0f;light_ambient[1]= 0.0f;light_ambient[2]= 0.0f;light_ambient[3]=1.0f;	
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);	
	//glRotatef(rotatepaint, rotatepaint / 1, rotatepaint / 3, rotatepaint / 2);
	rotatepaint +=0.1;
	//glRotatef(rotatepaint, rotatepaint , rotatepaint, rotatepaint);
	
	glRotatef(-24+rotatepaint*9, 60, 80,32);	
	
	glutSolidTorus(10,20,10,10); 


	glPopMatrix();
	mat_ambient[0]=0.2f;mat_ambient[1]= 0.2f;mat_ambient[2]= 0.2f;mat_ambient[3]=0.5f;	
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
}

void drawline(float w, float r, float g, float b,float x1, float y1, float z1,float x2, float y2, float z2) {
	glLineWidth(w);  glColor3f(r, g, b);
	glBegin(GL_LINES);
		glVertex3f(x1, y1, z1);
		glVertex3f(x2, y2, z2);
	glEnd();
}

void drawreference() {
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_QUADS);
	glVertex3f(0, 0, -0.002);
	glVertex3f(297, 0, -0.002);
	glVertex3f(297, -210, -0.002);
	glVertex3f(0, -210, -0.002);
	glEnd();
	drawline(5, 1, 0, 0, 0, 0, 0, 297, 0, 0);
	drawline(5, 0, 1, 0, 0, 0, 0, 0, -210, 0);
	drawline(5, 0, 0, 1, 0, 0, 0, 0, 0, 210);
}


double distance(double x1, double y1, double z1,
	double x2, double y2, double z2) {
	return (sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)));
}

// a useful function for displaying your coordinate system
void drawAxes(float length) {
	glLineWidth(1);
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

void drawAxes2(float length,float x1,float y1,float z1) {
	glLineWidth(1);
	glPushAttrib(GL_POLYGON_BIT | GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT);
	//glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3f(1, 0, 0);
	glVertex3f(x1, y1, z1);
	glVertex3f(x1+length, y1, z1);
	glColor3f(0, 1, 0);
	glVertex3f(x1, y1, z1);
	glVertex3f(x1, y1+length, z1);
	glColor3f(0, 0, 1);
	glVertex3f(x1, y1, z1);
	glVertex3f(x1, y1, z1+length);
	glEnd();
	glPopAttrib();
}


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
	boardSize=Size(5, 6);

	for (int i = 0; i < boardSize.height; ++i)
		for (int j = 0; j < boardSize.width; ++j) {
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
	//glFlush();

}

void drawtea2(int col, int row) {

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
	boardSize = Size(4, 3);

	for (int i = 0; i < 1;/* boardSize.height;*/ ++i)
		for (int j = 0; j < 1;/*boardSize.width;*/ ++j) {
			if (i == col && j == row)continue;
			glPushMatrix();
			glTranslatef(j*squareSize2 - 3 * squareSize2*teapot, i*squareSize2 - 2 * squareSize2*teapot, 0.0);
			if (j % 2 == 0) {

				glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
				glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
				glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
				glLightfv(GL_LIGHT0, GL_POSITION, light_position);
				glutSolidSphere(squareSize2 / 3, 20, 200);

				glEnable(GL_BLEND);
				glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient2);
				glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse2);
				glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular2);
				glLightfv(GL_LIGHT0, GL_POSITION, light_position2);
				glutSolidCube(squareSize2 / 3 * 2.1);
				glDisable(GL_BLEND);

			}
			if (j % 2 == 1) {
				glRotatef(90, 1, 0, 0);
				glRotated(180, 0, 1, 0);
				glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient3);
				glutSolidTeapot(squareSize2 / 3);
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
	//glFlush();
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
}

void getextrinsicparameters(Mat &cameraMatrix, Mat &distCoeffs) {
	if (nombrevideo == "video2.mp4") {
		//cte 110 0.30
		/*
		cameraMatrix = (Mat_<double>(3, 3) << 1112.845179694531, 0, 857.1091763643549,
			0, 1109.979207704567, 555.251348969829,
			0, 0, 1);
		distCoeffs = (Mat_<double>(8, 1) << -0.09657156041303895,
			0.8647468107514414,
			0,
			0,
			-1.307334638548276, 0, 0, 0);
		*/
		//cte 90 0.25

		cameraMatrix = (Mat_<double>(3, 3) << 909.6455568576023, 0, 701.4340981783473,
			0, 907.3177208413616, 454.3474487800865,
			0, 0, 1);
		distCoeffs = (Mat_<double>(8, 1) << -0.09864424280708323,
			0.8778651983833405,
			0,
			0,
			-1.338511223223186, 0, 0, 0);


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
	cout << "sisisisize " << objectPoints.size() << " " << imagePoints.size() << endl;
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

vector<Point2f> leer12(){
	double lo = 1; double hi = 150;
	vector<Point2f>imgpoints;

	vector<nodo>aux;
	for (int it = 0; it < 10; it++){
		double me = (lo + hi) / 2;
		memset(visited, 0, sizeof(visited));
		int maxi = 0;
		for (int i = 0; i < v.size(); i++){
			if (!visited[i]){
				aux.clear();
				maxi = max(maxi, bfs(i, me, aux));
			}
		}

		if (maxi >= 12){
			hi = me;
		}
		else{
			lo = me;
		}
	}


	memset(visited, 0, sizeof(visited));
	for (int i = 0; i < v.size(); i++){
		if (!visited[i]){
			aux.clear();
			int val = bfs(i, hi, aux);
			if (val == 12){
				break;
			}
		}
	}

	v = aux;
	std::sort(v.begin(), v.end());

	if (v.size() != 12)return imgpoints;

	Size boardSize(3, 4);
	vector<vector<nodo> >vsegmento;

	if(v.size() == 12){
		vector<Point3f> corners;
		for (int i = 0; i < boardSize.width; ++i)
			for (int j = 0; j < boardSize.height; ++j)
				corners.push_back(Point3f(float(j*squareSize2), float(i*squareSize2), 0));
		
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
			std::sort(puntos4.begin(), puntos4.end());
			// luego de obtener la recta con 4 puntos lo marcamos como visitados
			// Para que no obtener la misma recta en la siguiente iteracion
			for (int i = 0; i < puntos4.size(); i++) {
				v[puntos4[i].second].visit = caso + 1;
				aux.push_back(v[puntos4[i].second]);
			}
			vsegmento.push_back(aux);
		}
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
			imgpoints.push_back(Point2f((float)vsegmento[2 - i][j].cx, (float)vsegmento[2 - i][j].cy));
		}
	}

	//drawChessboardCorners(viewgl, boardSize, Mat(imgpoints), 1);
	//cv::imshow("Crocante2", viewgl);

	return imgpoints;
}

vector<Point2f>leer30(){
	vector<Point2f>imgpoints;

	double lo = 1; double hi = 150;

	vector<nodo>aux;
	for (int it = 0; it < 10; it++){
		double me = (lo + hi) / 2;
		memset(visited, 0, sizeof(visited));
		int maxi = 0;
		for (int i = 0; i < v.size(); i++){
			if (!visited[i]){
				aux.clear();
				maxi = max(maxi, bfs(i, me, aux));
			}
		}
		
		if (maxi >= 30){
			hi = me;
		}else{
			lo = me;
		}
	}


	memset(visited, 0, sizeof(visited));
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

	//El tablero es de 5X6
	Size boardSize(5, 6);
	//cout << "elementos= " << v.size() << endl;
	std::sort(v.begin(), v.end());

	if (v.size() != 30)return imgpoints;

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

						std::sort(dis.begin(), dis.end());
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
			std::sort(puntos6.begin(), puntos6.end());
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

		//drawChessboardCorners(viewgl, boardSize, Mat(pointBuf2), 1);
		//cv::imshow("Crocante", viewgl);
		
		return pointBuf2;
	}

	return imgpoints;
}

void display() {
	/* ---------------------------------------------------- 1ERA PARTE ---------------------------------------------*/
	//1ERA PARTE del trabajo corresponde en obtener el tablero de anillos (total 30) sin ruidos
	//vector<vector<Point2f> > imagePoints;
	//vector<pair<double, int> >errorframe;
	numframe++;
	if(numframe==1)LoadGLTextures();		
	
	int contx = 0;
	Mat cameraMatrix = Mat::eye(3, 3, CV_64F);
	Mat distCoeffs = Mat::zeros(8, 1, CV_64F);
	//Mat original = viewgl.clone();
	// Mat original contiene los frames del video a colores
	resize(viewgl, viewgl, Size( (1920*cte)/120, (1080*cte)/120));

	imageSize = viewgl.size();

	//cv::imshow("images1", viewgl);

	//cout << imageSize << endl;
	//cout << numframe << " " << readlines << endl;
	//imshow("viewgl ", viewgl);
	int idfaltante = -1;
	int id3 = -1;
	vector<Point2f>img30;
	vector<Point2f>img12;

	if (!lectura) {
		//convertimos a escala de grises
		//cv::imshow("images", viewgl);
		Mat destviewgl;
		
		cv::cvtColor(viewgl, src, CV_BGR2GRAY);
		//eliminamos el ruido con desenfoque gaussiano
		GaussianBlur(src, src, Size(3, 3), 11, 11);
		//Obtenemos los bordes de la imagen con Canny
		if (imageSize.width > 700) {
			if (nombrevideo != "PadronAnillos_03.avi") {
				getthreshold = 1;
				//adaptiveThreshold(src, src, 255.0, CV_THRESH_BINARY, CV_THRESH_BINARY, 81, 3);
				//imshow("Adaptative", src);
			}
		}
		
		cv::Canny(src, src, 80, 2 * 80);
		//cv::imshow("Canny", src);
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;
		//Encontramos todos los contornos de la imagen y lo guardamos en el vector Contours
		cv::findContours(src, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
		/// Draw contours
		//Mat drawing2 = Mat::zeros(src.size(), CV_8UC3);
		Mat drawing = Mat::zeros(src.size(), CV_8U);
		Mat drawing2 = Mat::zeros(src.size(), CV_8U);

		// vector v contendra informacion de los centros de los contornos validos
		v = vector<nodo>(contours.size(), nodo(1, 1, 1, 1, 1));
		for (int ii = 0; ii < (int)contours.size(); ii++) {
			vector<Point>P = contours[ii];
			Mat pointsf;
			RotatedRect box;
			// Un contorno es valido si el area es positiva y no sea tan grande
			if (contourArea(contours[ii])>30 && contourArea(contours[ii]) < 7000 && contours[ii].size() >= 5) {
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
				if (fabs(f2 - f1) <= 100 && f1*f2 <= 8000) {
					Point2f centro = box.center;
					v[ii] = nodo(ii, 0, centro.x, centro.y, f1*f2);
				}
			}
		}
		
		v = limpiar(v);
		for (int i = 0; i < v.size(); i++)
			drawContours(drawing2, contours, v[i].id, Scalar(255, 255, 255), CV_FILLED);
		//cv::imshow("drawing2: ", drawing2);
		
		std::sort(v.begin(), v.end());
		set<nodo>S;
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
				if (abs((v[i].cx - v[j].cx)) <= 7.5 && abs(v[i].cy - v[j].cy) <= 7.5 && v[i].area>v[j].area) {
					visited[j] = 1;
					sumx += v[i].cx;
					sumy += v[i].cy;
					cont++;
				}
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
				if (i != j && abs((v[i].cx - v[j].cx)) <= 7.5 && abs(v[i].cy - v[j].cy) <= 7.5 && v[i].area>v[j].area)
					v[j].visit = 1;
				
		//eliminaremos todos los elementos que tengan el atributo visit diferente de 0
		v = limpiar(v);
		
		for (int i = 0; i < v.size(); i++)
			drawContours(drawing, contours, v[i].id, Scalar(255, 255, 255), CV_FILLED);
		
		//cv::imshow("drawing", drawing);

		vector<nodo>totalv = v;
		img30 = leer30();
		if (img30.size() == 30){
			/*
			vector<Point3f> objcorners;
			calcBoardCornerPositions(Size(5, 6), squareSize, objcorners);
			framecorrecto++;
			vector<Mat>rvecs;
			vector<Mat>tvecs;
			*/
			//vector<vector<Point2f> > auximgp(1, img30);
			//vector<vector<Point3f> > auxobjcorners(1, objcorners);
			//double rms = calibrateCamera(auxobjcorners, auximgp, imageSize, cameraMatrix,
			//	distCoeffs, rvecs, tvecs, CV_CALIB_ZERO_TANGENT_DIST);
			//myfile << numframe << " " << rms << endl;
			//printf("-------------RMS     %d %.5lf\n", numframe, rms);
			/*
			cout << "posnumframe<<<< " << posnumframes << endl;
			if (numframes[posnumframes] == numframe){
				imgp.push_back(img30);
				objpoints.push_back(objcorners);
				posnumframes++;
			}

			if (posnumframes == 75){
				
				runCalibrationAndSave(imageSize, cameraMatrix, distCoeffs, imgp,objpoints);
				return;
			}
			*/
		}

		
		v = totalv;
		
		for (int i = 0; i < v.size(); i++)
			for (int j = 0; j < img30.size(); j++)
				if (v[i].cx == img30[j].x && v[i].cy == img30[j].y)
					v[i].visit = 1;

		v = limpiar(v);
		img12 = leer12();


	}else{
		img30.clear();
		img12.clear();
		readpoints(img30,img12);	
	}

	if (img30.size() != 30){ /*myfile << "" << endl; myfile << "" << endl;*/ return; }
	if (img12.size() != 12){ /*myfile << "" << endl; myfile << "" << endl;*/ return; }

	if(patroninix==-1<<30){
		patroninix=img12[0].x;
		patroniniy=img12[0].y;
		patronx=patroninix;
		patrony=patroniniy;
		inidirx= img12[4].x-img12[0].x;
		inidiry=img12[4].y-img12[0].y;
		//for(int i=0;i<12;i++)
		//	cout<<img12[i].x<<" "<<img12[i].y<<endl;
	}else{
		patronx=img12[0].x;
		patrony=img12[0].y;
		dirx= img12[4].x-img12[0].x;
		diry=img12[4].y-img12[0].y;

		//cout<<"dir "<<atan(diry/dirx)*180/acos(-1)<<" "<<dirx<<" "<<diry<<endl;
		//cout<<"patron "<<patronx<<" "<<patrony<<endl;
	}

	for(int i=0;i<img30.size();i++){
		if(distance(img30[i].x,img30[i].y,0,purplex,purpley,0)<40 ){
			color=4;
		}
		else if(distance(img30[i].x,img30[i].y,0,verdex,verdey,0)<40 )
			color=1;
		else if(distance(img30[i].x,img30[i].y,0,angryx,angryy,0)<40 && numframe-lastframe>40){
			lastframe=numframe;			
			if(idangry==0){
				color=0;
				idangry=1;			
			}else{
				if(idangry==1){
					color=2;
					idangry=2;
				}else{
					if(idangry==2){
						color=3;
						idangry=0;					
					}				
				}
			}
		}
	}

	/*
	for (int i = 0; i < 30; i++)
		myfile << img30[i].x << " " << img30[i].y << " ";
	myfile << endl;
	for (int i = 0; i < 12; i++)
		myfile << img12[i].x << " " << img12[i].y << " ";
	myfile << endl;
	*/
	glEnable(GL_DEPTH_TEST);  
	glClear(GL_COLOR_BUFFER_BIT);
	flip(viewgl, viewgl, 0);
	
	// OPENCV utiliza canales BGR y OpenGL RGB por lo tanto invertimos el canal B y R
	for (int i = 0; i < viewgl.size().width; i++)
		for (int j = 0; j < viewgl.size().height; j++)
			swap(viewgl.ptr<Vec3b>(j)[i][0], viewgl.ptr<Vec3b>(j)[i][2]);

	glLineWidth(2);
	glDrawPixels(viewgl.size().width, viewgl.size().height, GL_RGB, GL_UNSIGNED_BYTE, viewgl.ptr());
	//////////////////////Read Camera Matrix y Object points/////////////////////////////////
	
	vector<Point3f> corners;
	calcBoardCornerPositions(Size(5,6), squareSize, corners);
	getextrinsicparameters(cameraMatrix, distCoeffs);
	////////////////// Projection Matrix  Padron Superior ///////////////////////////////////////
	Mat projMat = Mat::zeros(4, 4, CV_64FC1);
	float zfar = 100000.0f, znear = 0.1f;
	projMat.at<double>(0, 0) = 2 * cameraMatrix.at<double>(0, 0) / viewgl.size().width;
	projMat.at<double>(0, 1) = 2 * cameraMatrix.at<double>(0, 1) / viewgl.size().width;
	projMat.at<double>(0, 2) = 1 - (2 * cameraMatrix.at<double>(0, 2) / viewgl.size().width);
	projMat.at<double>(1, 1) = 2 * cameraMatrix.at<double>(1, 1) / viewgl.size().height;
	projMat.at<double>(1, 2) = -1 + ((2 * cameraMatrix.at<double>(1, 2)) / viewgl.size().height);
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
	solvePnP(corners, img30, cameraMatrix, distCoeffs, rvec, tvec);
	Rodrigues(rvec, rotation);

	double daux = norm((-rotation).t()*tvec);
	cout << "distancia " << daux << endl;
	double t1=tvec.at<double>(0, 0);
	double t2=tvec.at<double>(1, 0);
	double t3=tvec.at<double>(2, 0);
	

	cout<<"superior "<<tvec.at<double>(0, 0)<<" "<<tvec.at<double>(1, 0)<<" "<<tvec.at<double>(2, 0)<<endl;
	
	for (int row = 0; row < 3; ++row) {
		for (int col = 0; col < 3; ++col)
			viewMatrix.at<double>(row, col) =  rotation.at<double>(row, col);
		viewMatrix.at<double>(row, 3) = tvec.at<double>(row, 0);
	}
	viewMatrix.at<double>(3, 3) = 1.0f;

	cv::Mat cvToGl = cv::Mat::zeros(4, 4, CV_64F);
	cvToGl.at<double>(0, 0) = 1.0f;
	cvToGl.at<double>(1, 1) = -1.0f;
	cvToGl.at<double>(2, 2) = -1.0f;
	cvToGl.at<double>(3, 3) = 1.0f;
	Mat glViewMatrix;
	viewMatrix = cvToGl* viewMatrix;
	transpose(viewMatrix, glViewMatrix);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glLoadMatrixd(&glViewMatrix.at<double>(0, 0));

	//drawAxes(500.0);

	//glScalef(1.0,-1.0,-1.0);
	//rook();
	if (teapot) {
		//drawtea(2 - id3, idfaltante);
	}
	
	if (sphere) {
		drawsphere(2 - id3, idfaltante);
	}



	////////////----------Ortogonal----------////////////
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			projectionMatrix[4 * j + i] = projMat.at<double>(i, j);
	glViewport(0, 0, viewgl.size().width, viewgl.size().height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, viewgl.size().width, 0, viewgl.size().height, -10000, 10000.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	drawAxes(viewgl.size().width);
	//glEnable(GL_BLEND);

		if (nlinesdraw == 0){
			iniciox = img30[0].x;
			inicioy = img30[0].y;
			inicioz = 0;
		}

		if (nlinesdraw == 0 || vlines[nlinesdraw - 1].pos[0] == 1 << 30){
			vlines[nlinesdraw].pos[0] = img30[0].x;
			vlines[nlinesdraw].pos[1] = img30[0].y;
			vlines[nlinesdraw].pos[2] = daux;
			vlines[nlinesdraw].color = color;
			vlines[nlinesdraw].drawon = drawon;

			nlinesdraw++;
			vlines[nlinesdraw].pos[0] = 1 << 30;
		}
		else{
			if (distance(vlines[nlinesdraw - 1].pos[0], img30[0].x, vlines[nlinesdraw - 1].pos[1], img30[0].y, 0, 0) >= 4){
				vlines[nlinesdraw].pos[0] = img30[0].x;
				vlines[nlinesdraw].pos[1] = img30[0].y;
				vlines[nlinesdraw].pos[2] = daux;
				vlines[nlinesdraw].color = color;
				vlines[nlinesdraw].drawon = drawon;
				nlinesdraw++;
				vlines[nlinesdraw].pos[0] = 1 << 30;
			}
		}


	drawPaint();
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
		
	delta+=0.1;

	glRotatef(anglerotatex,1.0f,0.0f,0.0f);
	glRotatef(anglerotatey,0.0f,1.0f,0.0f);
	glRotatef(anglerotatez,0.0f,0.0f,1.0f);
	//glTranslatef(cos(delta)*20,sin(delta)*20,delta*10);
	drawAxes2(80.0, 120, 120, 0);
	// 1 up  -x  2 +x  izq  3 -y  der 4 +y izq  5 -z diag-  6 +z diag+
	
	if( abs(patroninix-patronx)>=2 ){
		anglerotatex+=(patronx-patroninix)/5;
		patroninix=patronx;
	}

	if( abs(patroniniy-patrony)>=2 ){
		anglerotatey+=(patrony-patroniniy)/5;
		patroniniy=patrony;
	}


/*		
	patroninix=img12[0].x;
	patroniniy=img12[0].y;
	inidirx= img12[4].x-img12[0].x;
	inidiry=img12[4].y-img12[0].y;
	patronx=img12[0].x;
	patrony=img12[0].y;
	dirx= img12[4].x-img12[0].x;
	diry=img12[4].y-img12[0].y;
*/

	for (int i = 1; i < nlinesdraw; i++){
		if (vlines[i].drawon)
			drawline(1.5+(900-vlines[i].pos[2])/300, basecolor[vlines[i].color][0], basecolor[vlines[i].color][1], basecolor[vlines[i].color][2], vlines[i].pos[0], viewgl.size().height - vlines[i].pos[1], vlines[i].pos[2],
				vlines[i - 1].pos[0], viewgl.size().height - vlines[i - 1].pos[1], vlines[i - 1].pos[2]);
		}

	glPushMatrix();	
	//glDisable(GL_BLEND);



	/////////////////////////////////////////////////////////////////////////////////////////////////////01
	/*        Padron inferior        	*/
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			projectionMatrix[4 * j + i] = projMat.at<double>(i, j);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glFlush();
	
	calcBoardCornerPositions(Size(4, 3), squareSize2, corners);
	glViewport(0, 0, viewgl.size().width, viewgl.size().height);
	glLoadMatrixd(projectionMatrix);
	cvtColor(viewgl, gimage, CV_BGR2GRAY);
	viewMatrix = cv::Mat::zeros(4, 4, CV_64FC1);
	solvePnP(corners, img12, cameraMatrix, distCoeffs, rvec, tvec);
	Rodrigues(rvec, rotation);
	for (int row = 0; row < 3; ++row) {
		for (int col = 0; col < 3; ++col)
			viewMatrix.at<double>(row, col) = rotation.at<double>(row, col);
		viewMatrix.at<double>(row, 3) = tvec.at<double>(row, 0);
	}
	viewMatrix.at<double>(3, 3) = 1.0f;

	cout<<"Inferior "<<tvec.at<double>(0, 0)<<" "<<tvec.at<double>(1, 0)<<" "<<tvec.at<double>(2, 0)<<endl;

	double t4=t1-tvec.at<double>(0, 0);
	double t5=viewgl.size().height-(t2-tvec.at<double>(1, 0));
	double t6=t3-tvec.at<double>(2, 0);
	

	cvToGl.at<double>(0, 0) = 1.0f;
	cvToGl.at<double>(1, 1) = -1.0f;
	cvToGl.at<double>(2, 2) = -1.0f;
	cvToGl.at<double>(3, 3) = 1.0f;
	viewMatrix = cvToGl * viewMatrix;
	transpose(viewMatrix, glViewMatrix);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glLoadMatrixd(&glViewMatrix.at<double>(0, 0));

	if (teapot){
		glBegin(GL_QUADS);
		glColor3d(1, 1, 1);
		glVertex3f(-3*squareSize2, -2*squareSize2, 0);
		glColor3d(1, 1, 1);
		glVertex3f(7.5*squareSize2, -2*squareSize2, 0);
		glColor3d(1, 1, 1);
		glVertex3f(7.5*squareSize2, 5*squareSize2, 0);
		glColor3d(1, 1, 1);
		glVertex3f(-3*squareSize2, 5*squareSize2, 0);
		glEnd();
	}


	for (int i = 1; i < nlinesdraw; i++){
		if (vlines[i].drawon)
			drawline(1, 0.8, 0.8,0.8, (vlines[i].pos[0] - iniciox)/2+1*squareSize,  (vlines[i].pos[1]-inicioy)/2+1*squareSize , 0,
			(vlines[i - 1].pos[0] - iniciox )/2+ 1 * squareSize, (vlines[i - 1].pos[1] - inicioy)/2 + 1 * squareSize, 0);
	}

	drawAxes2(80.0, -3 * squareSize*teapot, -2 * squareSize*teapot, 0);
	
	drawtea2(2 - id3, idfaltante);
	
	glFlush();
	glutSwapBuffers();
	waitKey(1);
}

void reshape(int w, int h) {
	w = imageSize.width;
	h = imageSize.height;
	// set OpenGL viewport (drawable area)
	glViewport(0, 0, (1920 * cte) / 120, (1080 * cte) / 120);
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

	if (key == 'r'){ color = 0; }
	if (key == 'g'){ color = 1; }
	if (key == 'b'){ color = 2; }
	if (key == 'y'){ color = 3; }
	if (key == 'p'){ color = 4; }
	if (key == 'n'){ color = 5; }
	if (key == 'x'){ drawon = !drawon; }
	if (key == '2'){ anglerotatex += 1;cout<<"anglerorateeeeeex "<<anglerotatex<<endl; }
	if (key == '1'){ anglerotatex -= 1;cout<<"anglerorateeeeeex "<<anglerotatex<<endl; }
	if (key == '4'){ anglerotatey += 1;cout<<"anglerorateeeeeey "<<anglerotatey<<endl; }
	if (key == '3'){ anglerotatey -= 1;cout<<"anglerorateeeeeey "<<anglerotatey<<endl; }
	if (key == '6'){ anglerotatez += 1;cout<<"anglerorateeeeeez "<<anglerotatez<<endl; }
	if (key == '5'){ anglerotatez -= 1;cout<<"anglerorateeeeeez "<<anglerotatez<<endl; }
	if (key == '0'){ anglerotatex = 0;anglerotatey = 0;anglerotatez = 0;}

}

void idle() {
	// grab a frame from the camera
	(inputCapture) >> viewgl;
	glutPostRedisplay();
}

int main(int argc, char **argv) {
	// colocar PadronAnillos_01.avi , PadronAnillos_02.avi , PadronAnillos_03.avi,calibrar_anillo_nuevo_1280x720.wmv,calibrar_anillo_nuevo_640x360.wmv,medir_1280x720_anillos,medir_640x360_anillos.wmv que son los nombres de los videos
	const GLubyte *Vstr;
	Vstr = glGetString(GL_VERSION);
	fprintf(stderr, "Your OpenGL version is %s\n", Vstr);
	
	nombrevideo = "video2.mp4";
	int cantframes = 87;
	//myfile.open("vector.txt");
	myinfile.open("salida.txt");

	//myinfile.open("vector.txt");
	numframes = takeframes("video2.mp4", 75);
	//vector<int>frames = takeframes(nombrevideo, cantframes);// los frames a considerar para la calibracion
	//int p1[91] = { 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 28, 30, 32, 35, 38, 40, 42, 44, 46, 48, 50, 52, 54, 58, 102, 191, 193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223, 230, 238, 239, 241, 243, 245, 247, 249, 250, 253, 255, 257, 259, 262, 264, 318, 320, 321, 391, 393, 394, 396, 416, 419, 431, 433, 435, 450, 452, 454, 456, 474, 477, 479, 481, 482, 484, 486, 506, 508, 509, 519, 540, 543, 564, 567 };
	vector<int>frames;
	//for (int i = 0; i < cantframes; i++)frames.push_back(p1[i]);
	//vector<int>frames = primeroscantframes(cantframes);
	cout << "cantidad de frames a procesar" << cantframes << " " << frames.size() << endl;
	inputCapture = VideoCapture(nombrevideo);
	numframe = 0; // indica el numero de frame que esta siendo procesado en el video.
	int posframe = 0; // sirve para posicionar en el arreglo

	basecolor[0][0] = 255.0/255; basecolor[0][1] = 0; basecolor[0][2] = 0;//Red
	basecolor[1][0] = 0; basecolor[1][1] = 200.4/255; basecolor[1][2] = 0;//Green
	basecolor[2][0] = 0; basecolor[2][1] = 0; basecolor[2][2] = 204.0/255;//Blue
	basecolor[3][0] = 255.0 / 255; basecolor[3][1] = 255.0/255; basecolor[3][2] = 51.0/255;//Yellow
	basecolor[4][0] = 204.0 / 255; basecolor[4][1] = 0; basecolor[4][2] = 102.0/255;//Purple
	basecolor[5][0] = 0; basecolor[5][1] = 0; basecolor[5][2] = 0;//Black /Negro

	vlines[0].pos[0] = 1<<30;
	vector<vector<Point2f> > imagePoints;
	vector<pair<double, int> >errorframe;
	(inputCapture) >> viewgl;
	glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize((1920 * cte) / 120, (1080 * cte) / 120);
	glutCreateWindow("OpenGL");
	Mat undi, rview, map1, map2;
	// set up GUI callback functions
	glutDisplayFunc(display);
	//glutReshapeFunc(reshape);
	glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard);
	glutIdleFunc(idle);
	/**/
	// start GUI loop
	glutMainLoop();
	return 0;
}
