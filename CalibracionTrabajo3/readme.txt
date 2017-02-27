Ambos .cpp imagenecesCalibracionAnillos.cpp como imagenesCalibracionCirculos.cpp
leen de carpetas donde se encuentran frames seleccionados para calibrar

imagenesCalibracionAnillos.cpp
funciona con los videos
medir_640x360_anillos.wmv  -> Carpeta ---  imagenes_medir_640x360
PadronAnillos_03.avi	-> Carpeta ---     imagenes_medir_video3_anillos
calibrar_anillo_nuevo_640x360.wmv -> carpeta --- calibrar_anillo_nuevo_640x360 
calibrar_anillo_nuevo_1280x720.wmv -> carpeta ---calibrar_anillo_nuevo_1280x720

imagenesCalibracionCirculos.cpp
2 videos de calibrar y el video 3

-------------------------COMO  CORRER ------------------
nombrevideo = "medir_640x360_anillos.wmv";
int cantframes = 25; // 25 50 Y 75
--------------------------------------------------------

------------------------- QUE HACE -------------------------------
Lo que realiza es dado un conjunto de imagenes los cuales se encuentran en las carpetas
se realiza el proceso iterativo para el afinamiento.
-------------------------------------------------------------------


---------------------------TAREA-----------------------
Ejecutar para 25,50,75 frames para cada video
y sacar los RMS - hacer una grafica.

¿COMO HAGO ESO?
Ejecutalo, y revisa todas las itereaciones quedate con la que tuvo menor RMS
El comportamiento es baja iteracion tras iteracion y luego  comienza a subir y finalmente baja de manera despreciable
ejemplo revisa  "ejemplo.jpg"
nos quedamos en la iteracion 13 ya que nos da el menor RMS en la iteracion 14 comienza a subir
ojo que tambien tenemos como varia el camera matrix
--------------------------------------
