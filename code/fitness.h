#include "datastructure.h"
#include <algorithm>
#include <assert.h> 
#include <cmath>
#include <stdio.h>

#define KAPPA 0.05 /* variable de escala */ 
#define RHO 1.1 /* variable auxiliar de hioervolumen */

int dominates(single &a, single &b);
void compute_fitness(vector<single>  &population,vector<vector<double> > &indicators, vector<range> &bounds); 
static void compute_fitnesscomponents(vector<single>  &population, vector<vector<double> > &indicators, vector<range> &bounds); 
static void compute_bounds(vector<single>  &population, vector<range> &bounds); 
static double compute_indicator(vector<range> &bounds, single a, single b); 
static double compute_hypervolume(vector<range> &bounds, single a, single b, int d, bool flag); 