Prueben con el video calibrar_anillo_nuevo_640x360.wmv y calibrar_anillo_nuevo_1280x720.wmv

Quiero poner fichas de ajedrez en el padron y que este rote de acuerdo al padron.
Todos los .obj estan en la carpeta para la muestra solo he usado con la torre pueden agregar los demas solo deben jugar con la funcion 
scalef(si no le pones sale muy chiquito casi no se ve ponganle 20) y translatef revisen el codigo

Para las fichas se lee de un .obj  glm.cpp y gl.h tienen la implementacion de como leer.
"Queria usar la libreria Assimp pero no se como instalarlo en windows"
link de como usar los .obj
https://www.d.umn.edu/~ddunham/cs5721f07/schedule/resources/lab_opengl07.html

Consideracion a tomar en OPENGL
Opencv usar el formato BGR y OPENGL RGB
por eso intercambio los canales antes de mostrar sino chichico aparecera como pitufo
for (int i = 0; i < tempimage.size().width; i++)
    for (int j = 0; j < tempimage.size().height; j++)
	swap(tempimage.ptr<Vec3b>(j)[i][0], tempimage.ptr<Vec3b>(j)[i][2]);


http://spottrlabs.blogspot.pe/
modelView = [1   0   0   0]   *   [r00   r01   r02   t1]
            [0  -1   0   0]       [r10   r11   r12   t2]
            [0   0  -1   0]       [r20   r21   r22   t3]
            [0   0   0   1]       [0     0     0     1]

link con muchos links: http://stackoverflow.com/questions/21997021/augmented-reality-openglopencv

funciones utilizadas

funcion solvePNP de opencv
permite obtener las matrices de rotacion y traslacion a partir de cameramatrix , imagepoints ,objectpoints (osea  parametros intrinsecos y extrinsecos) calculados previamente
solvePnP(corners, imgpoints, cameraMatrix, distCoeffs, rvec, tvec);

Luego La matriz view lo calculamos asi --- revisen el codigo

for (int row = 0; row < 3; ++row) {
	for (int col = 0; col < 3; ++col) {
		viewMatrix.at<double>(row, col) = rotation.at<double>(row, col);
	}
	viewMatrix.at<double>(row, 3) = tvec.at<double>(row, 0);
}


funcion gluperspective
http://stackoverflow.com/questions/16571981/gluperspective-parameters-what-do-they-mean

funcion flip
flip(viewgl, tempimage, 0);
Si no coloco esto la imagen estara de cabeza al procesar la imagen de opencv a opengl
Es decir al imprimir una imagen de tipo MAT de opencv a opengl

funcion glDrawPixels
permite dibujar en OPENGL 
glDrawPixels(tempimage.size().width, tempimage.size().height, GL_RGB, GL_UNSIGNED_BYTE, tempimage.ptr());

Revisen las fotos resultado.png y resultado2.png

GLUT : inicializacion
	glutInit(&argc, argv);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(viewgl.cols, viewgl.rows);
	glutCreateWindow("OpenGL");

GLUT:
// set up GUI callback functions
	glutDisplayFunc(display);  // muestra las ventanas a mostrar
	glutReshapeFunc(reshape);  // redimensionar el ancho y largo de la ventana
	glutMouseFunc(mouse);     // Eventos por mouse
	glutKeyboardFunc(keyboard); //Eventos por teclado
	glutIdleFunc(idle);  // Lectura de frame a frame




funciones utilizadas para dibujar lo adicional (tetera u objetos de ajedrez)
// dibuja la torre a partir del rook.obj  
void rook(){
	glPushMatrix();tarea
	glColor3f(0.2, 0.7, 0.4);
	glTranslatef(115.0f, 115.0f, 0.0);
	glRotatef(-90, 0, 0, 1);          
	glScalef(10.0f, 10.0f, 10.0f);
	GLMmodel *objmodel_ptr;
	objmodel_ptr = glmReadOBJ("rook.obj");
	glmUnitize(objmodel_ptr);
	glmFacetNormals(objmodel_ptr);
	glmVertexNormals(objmodel_ptr, 90.0);
	glmDraw(objmodel_ptr, GLM_SMOOTH | GLM_MATERIAL | GLM_COLOR | GLM_TEXTURE);

	glPopMatrix();
	//glFlush();
	//glutSwapBuffers();
}

// dibuja una tetera mediante una funcion propia de GLU
void drawtea(){
	glPushMatrix();
	glColor3f(0.5, 0.6, 0.8);
	glTranslatef(60.0f, 60.0f, 0.0);
	glRotatef(-90,0, 0, 1);           // top teapot

	glutWireTeapot(10.0);
	glRotatef(-90, 0., 1., 0.);
	glutSolidTeapot(1);
	
	glPopMatrix();
	glFlush();
	//glutSwapBuffers();
}



TAREAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
USAR TEXTURAS PARA QUE LOS OBJETOS SE VEAN EN 3D

