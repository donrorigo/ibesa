#ifndef MUTATIONS_H
#define MUTATIONS_H

#include "datastructure.h"
#include "mhd.h"
#include "cai.h"
#include "lrcs.h"
#include <iostream>
#include <ctype.h>
#include <set>

void random_mutation(const single & a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss); 
void mhd_mutation(const single & a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss);
void cai_mutation(const single & a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss);
void undue_cai_mutation(const single & a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss);
void lrcs_mutation(const single & a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss);
void greedy_mutation(const single & a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss);

/* módulos auxiliares */
string change_CDS(string cds, string new_codon, int index);
void update_CDSs(vector<string> CDSs, string new_cds, string curr_cds, vector<string> &new_vector);
void update_vector(vector<string> &src, vector<string> &dst);

#endif