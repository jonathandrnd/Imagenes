
#include "stdafx.h"
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;
using namespace std;

class Settings{
public:
	enum Pattern { NOT_EXISTING, CHESSBOARD, CIRCLES_GRID, ASYMMETRIC_CIRCLES_GRID };
	enum InputType { INVALID, CAMERA, VIDEO_FILE, IMAGE_LIST };

	Size boardSize;				// Tamaño del patron -> Numero de elemntos en ancho y alto
	Pattern calibrationPattern;	// Para seleccionar el tipo de patron que se usara Chessboard, circulos, o circulos asimetricos
	float squareSize;			// El tamaño de un cuadrado en las unidades definidas (puntos, milimetros, etc)
	int delay;					// En el caso de que la entrada sea un video
	string outputFileName;		// El nombre del archivo de salida
	string input;

	int cameraID;
	vector<string> imageList;
	int atImageList;
	VideoCapture inputCapture;
	InputType inputType;
	bool goodInput;
	int flag;

	// Constructor
	Settings() : goodInput(false) {}

	void read(const FileNode& node){
		node["BoardSize_Width"]	  >> boardSize.width;
		node["BoardSize_Height"]  >> boardSize.height;
		node["Calibrate_Pattern"] >> patternToUse;
		node["Square_Size"]		  >> squareSize;
		node["Input"]			  >> input;
		node["Input_Delay"]		  >> delay;
		interprate();
	}

	void interprate(){

		goodInput = true;
		if (boardSize.width <= 0 || boardSize.height <= 0){
			cerr << "Tamagno de tablero no valido: " << boardSize.width << " " << boardSize.height << endl;
			goodInput = false;
		}

		if (squareSize <= 10e-6){
			cerr << "Tamagno de cuadrado invalido." << squareSize << endl;
			goodInput = false;
		}

		if (input.empty())			// Se verifica que la entrada sea valida
			inputType = INVALID;
		else{
			if (input[0] >= '0' && input[0] <= '9'){
				stringstream ss(input);
				ss >> cameraID;
				inputType = CAMERA;
			}
			else if (input.substr(input.size() - 4, 4) != ".xml" && input.substr(input.size() - 4, 4)!=".yml") {
				inputType = VIDEO_FILE;
			}
			else if (readStringList(input, imageList)){
					inputType = IMAGE_LIST;
			}
			else
				inputType = VIDEO_FILE;
			
			if (inputType == CAMERA){
				inputCapture.open(0);
			}

			if (inputType == VIDEO_FILE)
				inputCapture.open(input);

			if (inputType != IMAGE_LIST && !inputCapture.isOpened())
				inputType = INVALID;
		}

		if (inputType == INVALID){
			cerr << " No existe la entrada: " << input;
			goodInput = false;
		}

		flag = 0;

		calibrationPattern = NOT_EXISTING;
		if (!patternToUse.compare("CHESSBOARD")) calibrationPattern = CHESSBOARD;
		if (!patternToUse.compare("CIRCLES_GRID")) calibrationPattern = CIRCLES_GRID;
		if (!patternToUse.compare("ASYMMETRIC_CIRCLES_GRID")) calibrationPattern = ASYMMETRIC_CIRCLES_GRID;
		if (calibrationPattern == NOT_EXISTING){
			cerr << " Modo de calibracion de camara inexistente: " << patternToUse << endl;
			goodInput = false;
		}
		atImageList = 0;
	}

	Mat nextImage(){
		Mat result;
		if (inputCapture.isOpened()){
			Mat view0;
			inputCapture >> view0;
			view0.copyTo(result);
		}
		else if (atImageList < (int)imageList.size())
			result = imread(imageList[atImageList++], IMREAD_COLOR);

		return result;
	}

	static bool readStringList(const string& filename, vector<string>& l){
		l.clear();
		FileStorage fs(filename, FileStorage::READ);
		if (!fs.isOpened())
			return false;
		FileNode n = fs.getFirstTopLevelNode();
		if (n.type() != FileNode::SEQ)
			return false;
		FileNodeIterator it = n.begin(), it_end = n.end();
		for (; it != it_end; ++it)
			l.push_back((string)*it);
		return true;
	}

private:
	string patternToUse;

};

static void read(const FileNode& node, Settings& x, const Settings& default_value = Settings()){
	if (node.empty())
		x = default_value;
	else
		x.read(node);
}

enum { DETECTION = 0, CAPTURING = 1, CALIBRATED = 2 };

int main(){

	const char ESC_KEY = 27;
	const Scalar RED(0, 0, 255), GREEN(0, 255, 0);

	Settings s;
	const string inputSettingsFile = "in_VID5.xml";
	FileStorage fs(inputSettingsFile, FileStorage::READ); // Se leen las configuraciones
	if (!fs.isOpened()){
		cout << "No se puede abrir el archivo de configuracion: \"" << inputSettingsFile << "\"" << endl;
		return -1;
	}

	fs["Settings"] >> s;
	
	fs.release();                                         // Se cierra el archivo de configuraciones
	if (!s.goodInput){
		cout << "Se ha detectado una entrada no valida. El programa a terminado. " << endl;
		return -1;
	}

	vector<vector<Point2f> > imagePoints;
	Mat cameraMatrix, distCoeffs;
	Size imageSize;
	int mode = s.inputType == Settings::IMAGE_LIST ? CAPTURING : DETECTION;
	clock_t prevTimestamp = 0;
	
	while(true){
		Mat view;
		bool blinkOutput = false;

		view = s.nextImage();

		imageSize = view.size();  // Formato de imagen de entrada

		vector<Point2f> pointBuf;

		bool found = false;
		switch (s.calibrationPattern){ // Se encuentran los puntos caracteristicos en el formato de entrada
			case Settings::CHESSBOARD:
				found = findChessboardCorners(view, s.boardSize, pointBuf, CALIB_CB_ADAPTIVE_THRESH | CALIB_CB_FAST_CHECK | CALIB_CB_NORMALIZE_IMAGE);
				break;
			case Settings::CIRCLES_GRID:
				found = findCirclesGrid(view, s.boardSize, pointBuf);
				break;
			case Settings::ASYMMETRIC_CIRCLES_GRID:
				found = findCirclesGrid(view, s.boardSize, pointBuf, CALIB_CB_ASYMMETRIC_GRID);
				break;
			default:
				found = false;
				break;
		}
		
		if (found){
			// Para el caso de un patron de chessboard
			if (s.calibrationPattern == Settings::CHESSBOARD){
				Mat viewGray;
				cvtColor(view, viewGray, COLOR_BGR2GRAY);
				cornerSubPix(viewGray, pointBuf, Size(11, 11),
					Size(-1, -1), TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 30, 0.1));
			}

			if (mode == CAPTURING &&  // Solo para entrada de camara
				(!s.inputCapture.isOpened() || clock() - prevTimestamp > s.delay*1e-3*CLOCKS_PER_SEC))
			{
				imagePoints.push_back(pointBuf);
				prevTimestamp = clock();
				blinkOutput = s.inputCapture.isOpened();
			}

			// Dibujado de las esquinas
			drawChessboardCorners(view, s.boardSize, Mat(pointBuf), found);
		}
		
		imshow("Patron circulos", view);
		waitKey(1);
	}

	return 0;
}
