#include "lrcs.h"

/*
    Comentario: Tercer objetivo a cumplir en la resolución de la recombinación homóloga. 
    Autor: Rodrigo Puerto Pedrera
    Fecha: lun abr 15, 2019
    Fuente: Art. del tutor Miguel Ángel Vega Rodriguez - Multi-Objective Artificial Bee Colony for designing multiple genes
    encoding the same protein.
*/

using namespace std;

aim mlrs(vector<string> cds)
/* Devuelve el substring más grande EN UNA MISMA CADENA de todos los CDSs pasados por parámetro */
{

    double max = -1, l;
    aim result;
    CDSmaxLRS cdsmaxlrs;
 
    for(int i=0; i < (int)cds.size(); ++i) {
        #pragma omp critical
        cdsmaxlrs = lrs(cds[i]); 
        if((l = cdsmaxlrs.maxsubstring.length()) > max){
            max = l;     
            result.cds1 = cdsmaxlrs.maxsubstring;
            result.cds2 = cds[i];
            result.cds3 = "";
            result.index = cdsmaxlrs.index;
        }
    }

    result.value = max;
   
    return result;
}


aim mlrcs(vector<string> cds)
/* Devuelve el tamaño máximo de todos los substrings comunes de la proteina */
{
    aim result;
    
    double l = cds[0].length();
    double max = -1, temp;
    
    
    for(int i=0; i < (int)cds.size(); ++i) // numero de CDSs que tiene el vector
        for(int j=0; j < (int)l; ++j) // numero de caracteres del CDS
            cds[i][j] += 32;

    
    CDSmaxLCS cdsmaxlcs;         
    for(int i=0; i < (int)cds.size(); ++i)
        for(int j=i+1; j<(int)cds.size(); ++j){
            cdsmaxlcs = lcs(cds[i], cds[j]); 
            if((temp = cdsmaxlcs.maxsubstring.size()) > max){
                max = temp;
                result.cds1 = cdsmaxlcs.maxsubstring;
                result.cds2 = cds[i];
                result.cds3 = cds[j];
                result.index = cdsmaxlcs.index;
            }
        }

    aim mlr = mlrs(cds);
    
    if(mlr.value > max){
        max = mlr.value;
        result.cds1 = mlr.cds1;
        cout << "MISMA CADENA!" << endl;
    } 

    result.value = max/l;


    return result;
}