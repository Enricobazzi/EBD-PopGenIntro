# EBD - PopGenIntro

## Preparativos

Los scripts que vamos a correr para este curso necesitan muchas dependencias. La forma más sencilla de instalarlas todas es usando docker. Esto se puede conseguir bajando e instalando [Docker Desktop](https://www.docker.com/products/docker-desktop/).

Una vez instalado docker podemos proceder a descargar [esta carpeta](https://github.com/Enricobazzi/EBD-PopGenIntro/archive/refs/heads/main.zip) de github que contiene todo el material del curso.

Una vez descargada y descomprimida la colocamos en un sitio de nuestro ordenador que nos guste. A continuación vamos a tener que abrir una terminal y navegamos a la carpeta que acabamos de descomprimir con `cd`, nos bajamos la imagen del nuestro curso con `docker pull` y la corremos con `docker run`.

Copiando y pegando los comandos de abajo podemos tener todo listo. Solo tenemos que cambiar la ruta a la carpeta por la de nuestro ordenador:

**para Windows**

```
cd C:\path\to\your\folder
docker pull enricobazzi/ebd-popgenintro:v2-amd64
docker run -it -v C:\path\to\your\folder:/root enricobazzi/ebd-popgenintro:v2-amd64 /bin/bash
```

**para Linux y Mac**

```
cd /path/to/your/folder
docker pull enricobazzi/ebd-popgenintro:v2-amd64
docker run -it -v /path/to/your/folder:/root enricobazzi/ebd-popgenintro:v2-amd64 /bin/bash
```

Otra cosa que vamos a necesitar es un editor de código. Aconsejamos para este curso [Rstudio](https://posit.co/download/rstudio-desktop/), ya que nos va a permitir leer los markdown donde están escritos los comandos de las practicas y correr los scripts de R de forma interactiva.

Los paquetes de R necesarios para correr los scripts de R de esta practica están en el file [](). Ejecutando ese script o copiando su contenido en una terminal de R nos instalaremos todos los paquetes necesarios a analizar los datos y visualizar los resultados de las diferentes practicas.