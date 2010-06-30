#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H 1

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <string>

#include <arboretum/stMetricTree.h>
#include <arboretum/stPlainDiskPageManager.h>
#include <arboretum/stDiskPageManager.h>
#include <arboretum/stMemoryPageManager.h>
#include <arboretum/stSlimTree.h>
#include <arboretum/stMetricTree.h>

#include <arboretum/stUserLayerUtil.h>
#include <arboretum/stTypes.h>
#include <arboretum/stUtil.h>

#include "city.h"

typedef stMetricTree < TCity, TCityDistanceEvaluator > MetricTree;
typedef stSlimTree < TCity, TCityDistanceEvaluator > mySlimTree;

static const unsigned int C = 23;
static const unsigned int H = 4;
static const int FEATURE_SIZE = 5;
static const unsigned int NUMBER_OF_FEATURES = (H-1);
static const unsigned int GENOTYPE_SIZE = FEATURE_SIZE*NUMBER_OF_FEATURES;
static const unsigned int PHENOTYPE_SIZE = NUMBER_OF_FEATURES;

static const bool MAXIMIZATION = false;
static const int UPPER_LIMIT = C;
static const int LOWER_LIMIT = 0;
static const int DELTA = UPPER_LIMIT - LOWER_LIMIT;
static const unsigned int POPULATION_SIZE = 10;
static const unsigned int MAX_NUMBER_OF_GERATIONS = 10000;
static const double MUTATION_RATE = 0.5;
static const bool ELECTIVE = true;

#define DEBUG false
#define DEBUG_SELETION false
#define DEBUG_CROSSOVER false
#define DEBUG_MUTATION false

using namespace std;

typedef struct
{
	long genotype[GENOTYPE_SIZE];
	long phenotype[PHENOTYPE_SIZE];
	double fitness;
}Individual;

typedef struct
{
	Individual individuals[POPULATION_SIZE];
}Population;


static unsigned int number_of_gerations = 0;
static double media_of_fitness = 1;
static double old_media_of_fitness = media_of_fitness;
static double delta_fitness = 1;

void decoder(Individual *individual);
void coder(Individual *individual);
Individual * new_individual(long * features);
Population * new_population(void);
void print_individual(Individual *individual, const char *string);
void print_population(Population *population, const char *string, unsigned int number_of_geration);

void copy_individual(Individual* from, Individual* to);
void copy_population(Population* from, Population* to);
Population * new_empty_population(void);

void selection_by_tournament(Population* from, Population* to);
void crossover_uniform(Population* population);
void mutation(Population* population);

void set_phenotype(Population* population);
void set_genotype(Population *population);

Individual* better_individual(Population* population);
Individual* worse_individual(Population* population);
double get_media_of_population(Population* population);

void print_statistic(Population* pop, time_t begin, time_t end);

stResult<TCity> * main_genetic(mySlimTree* SlimTree, TCity * queryObject, stDistance range);

void evaluate_fitness(Population* population, mySlimTree* SlimTree, TCity * queryObject, stDistance range);
double objective_function(long* featues, mySlimTree* SlimTree, TCity * queryObject, stDistance range);

void infeasible(mySlimTree* SlimTree, Population* pop);

#endif
