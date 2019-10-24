#include "fitness.h"

/*
    Comentario: Tercer objetivo a cumplir en la resolución de la recombinación homóloga. 
    Autor: Rodrigo Puerto Pedrera
    Fecha: mie jun 12, 2019
    Fuente: Art. de Eckart Zitzler and Simon Künzli - Indicator-Based Selection in Multiobjective Search
*/

int dominates(single &a, single &b)
/* Routine for usual non-domination checking. It will return the following values:
   1 if a dominates b
   -1 if b dominates a
   0 if both a and b are non-dominated */
{
    
    int flaga=0,flagb=0;

    if (a.objetives[0] > b.objetives[0]) flaga=1;
    else if (a.objetives[0] < b.objetives[0]) flagb=1;

    if (a.objetives[1] > b.objetives[1]) flaga=1;
    else if (a.objetives[1] < b.objetives[1]) flagb=1;
        
    if (a.objetives[2] < b.objetives[2]) flaga=1;
    else if (a.objetives[2] > b.objetives[2]) flagb=1;

    if (flaga == 1 && flagb == 0) return 1;
    if (flaga == 0 && flagb == 1) return -1;
    return 0;
}

void compute_bounds(vector<single> &population, vector<range> &bounds)
/*  computes the bounds of the current objectives */
{
     
    int i, j, size = population.size();

    #pragma omp for
    for (i = 0; i < 3; ++i)
    {
	    bounds[0].min = population[0].objetives[i];
	    bounds[0].max = population[0].objetives[i];
    }

    for (i = 1; i < size; ++i){
	    for (j = 0; j < 3; ++j){
	        if (bounds[j].min > population[i].objetives[j]) bounds[j].min = population[i].objetives[j];
	        if (bounds[j].max < population[i].objetives[j]) bounds[j].max = population[i].objetives[j];
	    }
    }
    
    return;
}

double compute_hypervolume(vector<range> &bounds, single a, single b, int d, bool flag)
/* calculates the hypervolume of that portion of the objective space that
   is dominated by individual a but not by individual b */
{

    double A, B, r, max;
    double volume = 0;

    r = RHO * (bounds[d - 1].max - bounds[d - 1].min);
    max = bounds[d - 1].min + r;
    A = a.objetives[d - 1];
    if (flag) B = max;
    else B = b.objetives[d - 1];

    if (d == 1){
	    if (A < B) volume = (B - A) / r;
	    else volume = 0;
    }
    else
    {
	    if (A < B)
        {
	        volume = compute_hypervolume(bounds, a, b, d - 1, true) * (B - A) / r;
	        volume += compute_hypervolume(bounds, a, b, d - 1, false) * (max - B) / r;
	    }else volume = compute_hypervolume(bounds, a, b, d - 1, false) * (max - A) / r;   
    }

    return volume;
}

double compute_indicator(vector<range> &bounds, single a, single b)
/* calcula el valor del indicador, únicamente; no confundir con el sumatorio*/
{
    return (dominates(a, b) == 1) ? -compute_hypervolume(bounds, a, b, 3, false) : compute_hypervolume(bounds, b, a, 3, false);    
}

void compute_fitnesscomponents(vector<single> &population, vector<vector<double> > &indicators, vector<range> &bounds)
/*  calculates indicators of the entire population */
{
    double maximum_absindicator = 0;
    int i, j;
    int size = population.size();

    /* determine indicator values and their maximum */
    for (i = 0; i < size; ++i){
        #pragma omp for schedule(guided)
	    for (j = 0; j < size; ++j){
	        indicators[i][j] = compute_indicator(bounds, population[i], population[j]);
            #pragma omp critical
            {
                if (maximum_absindicator < fabs(indicators[i][j])) 
                maximum_absindicator = fabs(indicators[i][j]);   
            }
	    }
    }

    /* calculate for each pair of invidiuals the corresponding fitness component */
    #pragma omp for schedule(guided) collapse(2)
    for (i = 0; i < size; ++i)
        for (j = 0; j < size; ++j) indicators[i][j] = exp((-indicators[i][j]/maximum_absindicator)/KAPPA);
        
    return;
}

void compute_fitness(vector<single> &population, vector<vector<double> > &indicators, vector<range> &bounds)
/* gives fitness values to the entire population */
{
    int i, j, size = population.size();
    double sum;

    compute_bounds(population, bounds);

    compute_fitnesscomponents(population, indicators, bounds);

    for (i = 0; i < size; ++i){
        sum = 0;
        for (j = 0; j < size; j++) if(i != j)sum += indicators[j][i];
        population[i].fitness = sum;
    }
    
    return;
}