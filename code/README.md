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
  - Los resultados se tienen que convertir a formato de hypervolume: ./IBESA convert to HYPERVOLUME [file]
  - Despues se debe normalizar: ./normalize RESULTADOS\:\ [code of protein] [code of protein]
  - Calcular hipervolumen: ./hyp_ind hyp_ind_param.txt NORMALIZADO\:\ [code of protein] notnecessaryfile.txt results.txt
  - Convertir a decimal el valor: ./convert2decimal results.txt
