#include "mhd.h"

/*  Comentario: Segundo objetivo a cumplir en la resolución de la recombinación homóloga. 
    Autor: Rodrigo Puerto Pedrera
    Fecha: sáb mar 9, 2019
    Fuente: Art. del tutor Miguel Ángel Vega Rodriguez - Multi-Objective Artificial Bee Colony for designing multiple genes
    encoding the same protein   
*/

using namespace std;


double HD(string cds_i, string cds_j)
/* Calcula el número total de discrepancias de carácteres entre ambos CDSs */
{
    
    double result = 0, l = cds_i.length();
    for(int k = 0; k < l; k++) if(cds_i[k] != cds_j[k]) result++;
    return result;
}


aim mHD(vector<string> cds_vector)
/*  Obtiene el mínimo de todos los CDSs pasados por copia a la función */
{
    aim result;
    double min = cds_vector[0].length(), L = min, hd;

    for(int i=0; i < (int)cds_vector.size(); ++i)
        for(int j=i+1; j < (int)cds_vector.size(); ++j)
            if((hd = HD(cds_vector[i], cds_vector[j])) < min){
                min = hd;
                result.cds1 = cds_vector[i];
                result.cds2 = cds_vector[j];
            } 
    result.value = min/L;
    
    return result;
}