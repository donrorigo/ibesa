#include <vector>
#include <set>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <cmath>
#include <assert.h> 
#include <stdio.h>
#include <functional>
#include "mhd.h"
#include "cai.h"
#include "lrcs.h"
#include "fitness.h"
#include "mutations.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <omp.h>
#include "datastructure.h"


#define KAPPA 0.05 /* variable de escala */ 
#define RHO 1.1 /* variable auxiliar de hipervolumen */
#define OLD 3 /* numero de mutaciones donde se considera la solucion estancada */
#define PROBABILITY 5

using namespace std;

vector<range> bounds; /* array de límites: [0] = HD, [1] = CAI, [2] = LRCS */
vector<single> population;
vector<single> solutions;
vector<vector<double> > indicators;
vector<unsigned int> random_vector;
vector<vector<string> > auxiliar_cdss;
typedef function<void(single &,single &, int, int, vector<unsigned int> &, vector<vector<string> > &)> mutation_function;
vector<mutation_function> greedy_mutations;
vector<mutation_function> optimum_mutations;    
static set<string> fitness_vector;
static set<string> alive_vector;

static void init(int total, string amino_sequence, int CDSs, int machos);
static void show_cdss(vector<string> CDSs);
static void sorting_population(vector<single> &vector_pop);
static void write_results(string code);
static void show_population();
static void show_single(single a);
static void stabilize_population();
static void create_superCAI(int numC, string sequence);