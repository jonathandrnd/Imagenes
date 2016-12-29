#Padron de Calibracion

Trabajo N°1 del curso de Imagenes de la UCSP
----------------------------------------------------------------
Dependencias
----------------------------------------------------------------
Windows 2010
Visual Studio 2013
Opencv 2.4.13
OpenMP

----------------------------------------------------------------
Especificaciones
----------------------------------------------------------------
El trabajo se encuentra en 3 carpetas que realizan lo siguiente:
* Padron Anillos
	Contiene la deteccion del padron de calibracion de anillos, mediante un algoritmo propio.
* Padron Circulos OpenCV
	Contiene la deteccion del padron de calibracion de anillos, mediante el algoritmo de OpenCV.
* Padron Circulos
	Contiene la deteccion del padron de calibracion de circulos, mediante un algoritmo propio.
El Padron es calculado automaticamente.
	
----------------------------------------------------------------
Instalacion
----------------------------------------------------------------
* Descargar el OpenCV 2.4.13 http://opencv.org/downloads.html
* El proyecto fue realizado en 32 bits, Core I7 8 nucleos.
* Configurar OpenMP en Visual Studio. (Properties-> C/C++ ->Language -> OpenMP Support -> Yes)

----------------------------------------------------------------
Equipo Victor.
----------------------------------------------------------------
Integrantes

*Jonathan Durand

*Victor Cornejo

*Javier Neira


----------------------------------------------------------------
MEJORA, RESPECTO A LO ENVIADO EL 22 DE DICIEMBRE
----------------------------------------------------------------
Tanto en el padron de circulo y padron de anillos obtenemos robustamente el padron  y el zigzag en el padron de anillos.
Para el caso del padron de circulos el zigzag presentaba muchos errores

La mejora esta en el zigzag del padron de circulos, ahora se obtiene de manera mas robusta.
El cambio consiste en obtener el punto superior izquierda, el cual se obtiene por medio del convex hull ya que este forma un angulo cercano a 90 grados.
Luego a partir de el trazaremos una diagonal de tamaño 8 en total son 9 diagonales de tamaños (8,8,7,6,5,4,3,2,1).
ver "dibujo_prueba_diagonal1"

Luego de obtener las diagonales los ordenamos de abajo hacia arriba mediante producto vectorial
y lo transformamos para que el zigzag sea en vertical. 
El metodo no es robusto con el video 3.
