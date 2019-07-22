#include "cai.h"

/*
    Comentario: Primer objetivo a cumplir en la resolución de la recombinación homóloga. 
    Autor: Rodrigo Puerto Pedrera
    Fecha: sáb mar 9, 2019
    Fuente: Art. del tutor Miguel Ángel Vega Rodriguez - Multi-Objective Artificial Bee Colony for designing multiple genes
    encoding the same protein
*/

using namespace std;

double cai(string cds)
/* Calcula el peso total del cds proporcionado por copia */
{
    double total = 1, length = cds.size(), cexp = 1.0L/(length/3.0L);
    for(int i = 0; i<length; i+=3) total = total * amino_weights[cds.substr(i,3)];  
    return (pow(total, cexp)); 
}

double cai2(string cds)
/* Calcula el peso total del cds proporcionado por copia */
{
    double total = 0, length = cds.size();
    for(int i = 0; i<length; i+=3) total = total + log(amino_weights[cds.substr(i,3)]);    
    return (exp(total / (length/3.0)));
}

aim mCAI(vector<string> cds_vector)
/* Obtiene el mínimo de todos los CDSs pasados por copia a la función */
{
    
    double min = 1, ncai;
    aim result;
    
    for(int i = 0; i < (int)cds_vector.size(); i++){
        if((ncai = cai2(cds_vector[i])) <= min){
            result.value = ncai;
            result.cds1 = cds_vector[i];
            min = ncai;
        } 
    }
   
    return result;
}