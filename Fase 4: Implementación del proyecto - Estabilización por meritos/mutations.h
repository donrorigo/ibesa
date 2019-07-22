#ifndef MUTATIONS_H
#define MUTATIONS_H

#include "datastructure.h"
#include "mhd.h"
#include "cai.h"
#include "lrcs.h"
#include <iostream>
#include <ctype.h>
#include <set>

void random_mutation(single a, single &result, double Pm, int th, vector<unsigned int> &random_vector); 
void mhd_mutation(single a, single &result, double Pm, int th, vector<unsigned int> &random_vector);
void cai_mutation(single a, single &result, double Pm, int th, vector<unsigned int> &random_vector);
void undue_cai_mutation(single a, single &result, double Pm, int th, vector<unsigned int> &random_vector);
void lrcs_mutation(single a, single &result, double Pm, int th, vector<unsigned int> &random_vector);
void greedy_mutation(single a, single &result, double Pm, int th, vector<unsigned int> &random_vector);

/* módulos auxiliares */
string change_CDS(string cds, string new_codon, int index);
void update_CDSs(vector<string> CDSs, string new_cds, string curr_cds, vector<string> &new_vector);

#endif