#include "genetic_algorithm.h"

void decoder(Individual * individual)
{
	for(int i = 0; i < PHENOTYPE_SIZE ; i++)
	{
		long aux_number = 0;
        
		for(int j = i * FEATURE_SIZE; j < ((i * FEATURE_SIZE)+FEATURE_SIZE); ++j)
		{
			aux_number += pow(2,(abs((j % FEATURE_SIZE) - FEATURE_SIZE-1))) * individual->genotype[j];
		}
		/* TODO Enteder o Porquer */
		individual->phenotype[i] = aux_number/4;
	}
}

void coder(Individual * individual)
{
	long j = 0;
	for(int i = 0; i < PHENOTYPE_SIZE; ++i)
	{
   		long x, bit, nota;

   		x = individual->phenotype[i];

   		nota = pow(2, FEATURE_SIZE-1);
   		while(nota >= 1)
   		{
			
   			bit = x / nota;
   			individual->genotype[j] = bit;
			j++;
   			x = x % nota;
   			nota = nota / 2;
   		}
	}
}

Individual * new_individual(long * features)
{
	Individual * individual = (Individual *)malloc(sizeof(Individual));
	if(individual == NULL)
	{
		exit(1);
	}
	
	for(int i = 0; i < PHENOTYPE_SIZE; ++i)
	{
		individual->phenotype[i] = features[i];
	}
	return individual;
}

Population * new_population()
{
	Population * population;
	long features[PHENOTYPE_SIZE];
	if((population = (Population*)malloc(sizeof(Population))) != NULL)
	{
		for(int i = 0; i < POPULATION_SIZE; ++i)
		{
			for(int j = 0; j < PHENOTYPE_SIZE; ++j)
			{
				features[j] = rand() % (UPPER_LIMIT-LOWER_LIMIT) + LOWER_LIMIT;
			}
			population->individuals[i] = (*new_individual(features));
		}
	}
	else
	{
		exit(1);
	}
	return population;
}

Population * new_empty_population()
{
	Population * population;
	long features[PHENOTYPE_SIZE];
	if((population = (Population*)malloc(sizeof(Population))) != NULL)
	{
		for(int i = 0; i < POPULATION_SIZE; ++i)
		{
			for(int j = 0; j < PHENOTYPE_SIZE; ++j)
			{
				features[j] = 0;
			}
			population->individuals[i] = (*new_individual(features));
		}
	}
	else
	{
		exit(1);
	}
	return population;
}

void print_individual(Individual * individual, const char *string)
{
	cout << fixed;
	cout << "+---------------------------------------+" << endl;
	cout << "| " << string << ":" << endl;
	cout << "+---------------------------------------+" << endl;
	cout << "| Geneotype: ";
	for(int i = 0; i < GENOTYPE_SIZE; ++i)
	{
		if((i % FEATURE_SIZE) == 0)
		{
			cout << " ";
		}
		cout << individual->genotype[i];
	}
	cout << endl;
	cout << "| Pheneotype: ";
	for(int i = 0; i < PHENOTYPE_SIZE; ++i)
	{
		cout << individual->phenotype[i] << " ";
	}
	cout << endl;
	cout << "| Fitness: " << individual->fitness << endl;
	cout << "+---------------------------------------+" << endl;
}

void print_population(Population * population, const char *string, unsigned int number_of_geration = 0)
{
	cout << "+--------------------------------------------------------------------------+" << endl;
	cout << "| " << string << ", Number of Gerations: " << number_of_geration << endl;
	cout << "+--------------------------------------------------------------------------+" << endl;
	cout << "| NUM \tGeneotype: \t\tPheneotype: \tFitness: " << endl;
	cout << fixed;
	for(int i = 0; i < POPULATION_SIZE; ++i)
	{
		cout << "| [" << i << "]\t";
		for(int j = 0; j < GENOTYPE_SIZE; ++j)
		{
			if((j % FEATURE_SIZE) == 0)
			{
				cout << " ";
			}
			cout << population->individuals[i].genotype[j];
		}
		cout << "\t";
		for(int j = 0; j < PHENOTYPE_SIZE; ++j)
		{
			cout << population->individuals[i].phenotype[j] << " ";
		}
		cout << "\t\t";
		cout << population->individuals[i].fitness << endl;
	}
	cout << "+--------------------------------------------------------------------------+" << endl;
}
void copy_individual(Individual * from, Individual * to)
{
	for(int i = 0; i < GENOTYPE_SIZE; ++i)
	{
		to->genotype[i] = from->genotype[i];
	}
	for(int i = 0; i < PHENOTYPE_SIZE; ++i)
	{
		to->phenotype[i] = from->phenotype[i];
	}
	to->fitness = from->fitness;
}

void copy_population(Population * from, Population * to)
{
	for(int i = 0; i < POPULATION_SIZE; ++i)
	{
		copy_individual(&from->individuals[i], &to->individuals[i]);
	}
}

void selection_by_tournament(Population * from, Population * to)
{
	int index1, index2;
	for(int i = 0; i < POPULATION_SIZE; ++i)
	{
		index1 = rand() % POPULATION_SIZE;
		index2 = rand() % POPULATION_SIZE;

		if(from->individuals[index1].fitness >= from->individuals[index2].fitness)
		{
			copy_individual(&from->individuals[index1], &to->individuals[i]);
		} else
		{
			copy_individual(&from->individuals[index2], &to->individuals[i]);
		}
	}
}

void crossover_uniform(Population * population)
{
	long mask, index1, index2;
	long i, j, aux;
	for( i = 0; i < POPULATION_SIZE; i += 2)
	{
		index1 = i;
		index2 = i+1;

		for( j = 0; j < GENOTYPE_SIZE; j += 1)
		{
			mask = rand() % 2;
			if(mask == 0)
			{
				aux = population->individuals[index1].genotype[j];
				population->individuals[index1].genotype[j] = population->individuals[index2].genotype[j];
				population->individuals[index2].genotype[j] = aux;
			}
		}

		population->individuals[index1].fitness = 0;
		population->individuals[index2].fitness = 0;
		decoder(&population->individuals[index1]);
		decoder(&population->individuals[index2]);
	}
}

void mutation(Population *population)
{
	long mutation;
	long i, j;
	for(i = 0; i < POPULATION_SIZE; ++i)
	{
		for( j = 0; j < GENOTYPE_SIZE; ++j)
		{
			mutation = rand() % 100;
			if(mutation >= 0 && mutation < MUTATION_RATE)
			{
				if(population->individuals[i].genotype[j] == 0)
				{
					population->individuals[i].genotype[j] = 1;
				} else
				{
					population->individuals[i].genotype[j] = 0;
				}
			}
		}
	}
}

void set_phenotype(Population *population)
{
	for(int i = 0; i < POPULATION_SIZE; ++i)
	{
		decoder(&population->individuals[i]);
	}
}

void set_genotype(Population *population)
{
	for(int i = 0; i < POPULATION_SIZE; ++i)
	{
		coder(&population->individuals[i]);
	}
}

Individual* better_individual(Population* population)
{
	Individual * better_individual = &population->individuals[0];
	for (unsigned int i = 1; i < POPULATION_SIZE; i += 1)
	{
		if (better_individual->fitness < population->individuals[i].fitness)
		{
			better_individual = &population->individuals[i];
		}
	}
	return better_individual;
}

Individual* worse_individual(Population* population)
{
	Individual * worse_individual = &population->individuals[0];
	for (unsigned int i = 1; i < POPULATION_SIZE; i += 1)
	{
		if (worse_individual->fitness > population->individuals[i].fitness)
		{
			worse_individual = &population->individuals[i];
		}
	}
	return worse_individual;
}

double get_media_of_population(Population* population)
{
	double media = 0.0;
	for (unsigned int i = 0; i < POPULATION_SIZE; i += 1)
	{
		media += population->individuals[i].fitness;
	}
	media /= POPULATION_SIZE;
	return media;
}

void print_statistic(Population* pop, time_t begin, time_t end)
{
	cout << "+--------------------------------------------------------------------------+" << endl;
	cout << "| Statistics: " << endl;
	cout << "+--------------------------------------------------------------------------+" << endl;
	cout << "| Delta Fitness: " << delta_fitness << endl;
	print_individual(better_individual(pop), "Better Individual");
	print_individual(worse_individual(pop), "Worse Individual");
	cout << "| Timer in Second: " << difftime(end, begin) << endl;
	cout << "+--------------------------------------------------------------------------+" << endl;
}

stResult<TCity> * main_genetic(mySlimTree* SlimTree, TCity * queryObject, stDistance range)
{
	time_t begin, end;

	begin = time(NULL);
	
	srand(time(NULL));
	
	Population * pop = new_population();
	Population * aux_pop = new_empty_population();
	
	infeasible(SlimTree, pop);
	
	set_genotype(pop);
	
	evaluate_fitness(pop, SlimTree, queryObject, range);

	#if DEBUG
	print_population(pop ,"Pop", number_of_gerations);
	#endif

	while( (number_of_gerations < MAX_NUMBER_OF_GERATIONS) && (delta_fitness > 0.0001) )
	{
		number_of_gerations++;
		#if DEBUG
		cout << "Number of Gerations: "<< number_of_gerations << endl;
		#endif
		selection_by_tournament(pop, aux_pop);
		#if DEBUG
		cout << "Torneio" << endl;
		#endif
		crossover_uniform(aux_pop);
		#if DEBUG
		cout << "Crossover Uniforme" << endl;
		#endif
		mutation(aux_pop);
		#if DEBUG
		cout << "Mutação" << endl;
		#endif
		set_phenotype(aux_pop);
		#if DEBUG
		cout << "Fenotipo" << endl;
		#endif
		infeasible(SlimTree, aux_pop);
		#if DEBUG
		cout << "Infactibilidade" << endl;
		#endif
		set_genotype(aux_pop);
		#if DEBUG
		cout << "Fenotipo" << endl;
		#endif
		copy_population(aux_pop, pop);
		#if DEBUG
		cout << "Copiar aux pop" << endl;
		#endif
		evaluate_fitness(pop, SlimTree, queryObject, range);
		#if DEBUG
		cout << "Fitness" << endl;
		#endif
		
		old_media_of_fitness = media_of_fitness;
		media_of_fitness = get_media_of_population(pop);
		delta_fitness = fabs(media_of_fitness - old_media_of_fitness);
	}

	end = time(NULL);
	
	#if DEBUG
	print_population(pop, "PopulatioN", number_of_gerations);
	
	print_statistic(pop, begin, end);
	#endif
}

double objective_function(long * features, mySlimTree* SlimTree, TCity * queryObject, stDistance range)
{
	stPage * currPage;
	stSlimNode * currNode;
	stSlimLeafNode * leafNode;
	
	currPage = SlimTree->myPageManager->GetPage(SlimTree->GetRoot());
	currNode = stSlimNode::CreateNode(currPage);
	
	unsigned long numberOfEntries;
	double fitness = 0.0;
	for(int i = 0; i < PHENOTYPE_SIZE; ++i)
	{
		#if DEBUG
		cout << i << ": Feature: " << features[i] << endl;
		#endif
		stSlimIndexNode* indexNode = (stSlimIndexNode *)currNode;
		currPage = SlimTree->myPageManager->GetPage(indexNode->GetIndexEntry(features[i]).PageID);
		currNode = stSlimNode::CreateNode(currPage);
	}
	leafNode = (stSlimLeafNode *)currNode;
	TCity* city = new TCity();
	city->Unserialize(leafNode->GetObject(leafNode->GetRepresentativeEntry()), leafNode->GetObjectSize(leafNode->GetRepresentativeEntry()));
	TCityDistanceEvaluator * distanceEvaluator = new TCityDistanceEvaluator();
	stDistance distance = distanceEvaluator->GetDistance(queryObject, city);
	fitness = (double)distance;
	return fitness;
}

void evaluate_fitness(Population * population, mySlimTree* SlimTree, TCity * queryObject, stDistance range)
{
	double x, y;
	long i;
	for( i = 0; i < POPULATION_SIZE; ++i)
	{
		decoder(&population->individuals[i]);
		
		if(MAXIMIZATION)
		{
			population->individuals[i].fitness = objective_function(population->individuals[i].phenotype, SlimTree, queryObject, range);
		} else
		{
			population->individuals[i].fitness = (-1 * objective_function(population->individuals[i].phenotype, SlimTree, queryObject, range));
		}
	}
}

void infeasible(mySlimTree* SlimTree, Population* population)
{
	stPage * currPage;
	stSlimNode * currNode;
	
	unsigned long numberOfEntries;
	
	for(int i = 0; i < POPULATION_SIZE; i++)
	{
		currPage = SlimTree->myPageManager->GetPage(SlimTree->GetRoot());
		currNode = stSlimNode::CreateNode(currPage);
		for (int j = 0; j < PHENOTYPE_SIZE; j++)
		{
			stSlimIndexNode* indexNode = (stSlimIndexNode *)currNode;
			numberOfEntries = indexNode->GetNumberOfEntries();
		   	if(population->individuals[i].phenotype[j] >= numberOfEntries)
		   	{
		   		population->individuals[i].phenotype[j] = numberOfEntries-1;
		   	}
		   	currPage = SlimTree->myPageManager->GetPage(indexNode->GetIndexEntry(population->individuals[i].phenotype[j]).PageID);
		   	currNode = stSlimNode::CreateNode(currPage);
		}
		#if DEBUG
		cout << i << " -> ";
		register unsigned int j;
		for( j = 0; j < PHENOTYPE_SIZE; ++j)
		{
			cout << population->individuals[i].phenotype[j] << " ";
		}
		cout << endl;
		#endif
	}	
}