#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#include <map>
#include <vector>

using namespace std;

extern map<string, double> amino_weights;
extern map<char, vector<string> > amino_codons;
extern map <string, char> which_amino;

typedef struct single_st{
    bool gender;
    vector<string> cds;
    vector<double> objetives; /* [0] = HD, [1] = CAI, [2] = LRCS */
    double fitness;
    int age;
    string id;
    int lastmutation;
    int number; /* nº que lleva el elefante de la primera generación */
    bool operator==(const single_st& a) { return (id==a.id) ? true : false; }
    single_st& operator =(const single_st& a)
    {
        gender = a.gender;
        cds = a.cds;
        objetives[0] = a.objetives[0];
        objetives[1] = a.objetives[1];
        objetives[2] = a.objetives[2]; 
        fitness = a.fitness;
        age = a.age;
        id = a.id;
        lastmutation = a.lastmutation;
        number = a.number;
        return *this;
    }
} single; /* estructura de cada individuo */

typedef struct range_st{
    double min;
    double max;
} range;

typedef struct aim_st{
    string cds1; /* CAI = CDS con menos CAI :   mHD = CDSi : LRCS = substring */
    string cds2; /* CAI = no se utiliza :       mHD = CDSj : LRCS = CDS donde se encuentra el substring */
    string cds3;
    int index;
    double value;
} aim;

#endif // !DATASTRUCTURE_H