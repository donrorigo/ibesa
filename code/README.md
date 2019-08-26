# Multi-Objective ESA derived based on Indicators
## Información
Proyecto: Trabajo de Fin de Grado
Autor: Rodrigo Puerto Pedrera
Tutor: Miguel Ángel Vega Rodriguez

Grado en Ingeniería Informática en Ingeniería de los Computadores
Universidad de Extremadura
Escuela Politécnica

## Compilación: 
    - Si tenemos el archivo Makefile, únicamente invocaremos a la función: "make".
    - Si no la tenemos, invocaremos la siguiente línea: g++ ibesa.cpp lcs.cpp lrs.cpp mhd.cpp cai.cpp lrcs.cpp fitness.cpp mutations.cpp -o ibesa

## Ejecución: 
    - Para generar el fichero con toda la información de los elefantes: main.sh
    * Los resultados se tienen que convertir a formato de hypervolume: ./convert (japoneses) / ./IBESA convert to HYPERVOLUME --> file
    * Despues se debe normalizar: ./normalize file.txt proteins.txt (1=japoneses / 0=ibesa)
    * Convertir a decimal el valor: convert2decimal