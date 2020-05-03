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
#include <atomic>


#define KAPPA 0.05 /* variable de escala */ 
#define RHO 2 /* variable auxiliar de hipervolumen */
#define OLD 2 /* numero de mutaciones donde se considera la solucion estancada */
#define RANDOMMUTATION 2
#define GREEDYMUTATION 40


using namespace std;

string code;
vector<single> population; /* vector para la poblacion */
vector< vector <single>> optimum_mutated; /* vector auxiliar para las mutaciones optimas */
vector<single> solutions; /* vector grande de soluciones */
vector<single> paretofront; /* vector final con el frente de pareto */
vector<range> bounds; /* array de límites: [0] = HD, [1] = CAI, [2] = LRCS */
vector<vector<double> > indicators;


vector<unsigned int> random_vector;
vector<vector<string> > auxiliar_cdss;
typedef function<void(single &,single &, int, int, vector<unsigned int> &, vector<vector<string> > &)> mutation_function;
vector<mutation_function> greedy_mutations;
vector<mutation_function> optimum_mutations;    
static set<string> nonrepeat;
atomic_bool dominated; 

static void init(int total, string amino_sequence, int CDSs, int machos);
static void show_cdss(vector<string> CDSs);
static void sorting_population();
static void write_results();
static void show_population(int x);
static void create_superCAI(int numC, string sequence);
static void export2utility();
static int three_mutations(int i, int th);

/* atributos de debuggeo */
map <int, std::pair<int, int> > howmany; /* primera posición machos */
map <int, int > oldest; 
map <int, int > nonutility; 
map <int, int > optimum_utility;