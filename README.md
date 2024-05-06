# EBD - PopGenIntro

## Preparativos

Los scripts que vamos a correr para este curso necesitan muchas dependencias. La forma más sencilla de instalarlas todas es usando docker. Esto se puede conseguir bajando e instalando [Docker Desktop](https://www.docker.com/products/docker-desktop/).

Una vez instalado docker podemos proceder a descargar [esta carpeta](https://github.com/Enricobazzi/EBD-PopGenIntro/archive/refs/heads/main.zip) de github que contiene todo el material del curso.

Los datos que vamos a utilizar para las prácticas están disponibles para su descarga en [este enlace](https://saco.csic.es/index.php/s/Y9yw8XjJakc2gQY).

Una vez descargada y descomprimida la colocamos en un sitio de nuestro ordenador que nos guste. A continuación vamos a tener el programa de docker abierto a la vez que abrimos una terminal (o símbolo de sistema en Windows). Ahora nos podemos bajar la imagen del nuestro curso con `docker pull` y la corremos con `docker run`.

Copiando y pegando los comandos de abajo podemos tener todo listo. Solo tenemos que cambiar la ruta a la carpeta del curso por la de nuestro ordenador:


**para Windows** *si hemos guardado la carpeta del curso en "C:\ruta\a\la\carpeta\del\curso"*

[como extraer la ruta de una carpeta en windows?](https://www.sony-latin.com/es/electronics/support/personal-computers/articles/00015251)

```
docker pull enricobazzi/ebd-popgenintro:v2-amd64
docker run --security-opt seccomp=unconfined -it -v C:\ruta\a\la\carpeta\del\curso:/root enricobazzi/ebd-popgenintro:v2-amd64 /bin/bash
```

**para Linux y Mac** *si hemos guardado la carpeta del curso en "/ruta/a/carpeta/del/curso"*

[como extraer la ruta de una carpeta en mac?](https://iboysoft.com/es/como-hacer/copiar-ruta-de-archivo-en-mac.html)

[como extraer la ruta de una carpeta en linux?](https://www.hostinger.es/tutoriales/como-usar-comando-find-locate-en-linux/)

```
docker pull enricobazzi/ebd-popgenintro:v2-amd64
docker run --security-opt seccomp=unconfined -it -v /ruta/a/carpeta/del/curso:/root enricobazzi/ebd-popgenintro:v2-amd64 /bin/bash
```

Otra cosa que vamos a necesitar es un editor de código. Aconsejamos para este curso [Rstudio](https://posit.co/download/rstudio-desktop/), ya que nos va a permitir leer los markdown donde están escritos los comandos de las practicas y correr los scripts de R de forma interactiva.

Los paquetes de R necesarios para correr los scripts de R de esta practica están en el file [installation_packages.R](practica3/installation_packages.R) de la carpeta practica3. Ejecutando ese script o copiando su contenido en una terminal de R nos instalaremos todos los paquetes necesarios a analizar los datos y visualizar los resultados de las diferentes practicas.