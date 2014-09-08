#include <malloc/malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// definition of the boolean type (not pre-defined in C)
typedef int bool;
#define TRUE 1
#define FALSE 0

#define MAX_INT 2147483647

#define POPULATION 40
#define PROB_CROSSOVER 0.7
#define PROB_MUTATION 0.01
#define MAX_GENERATIONS 200

// I use this three variables to calculate random numbers
// For obtain always the same numbers don't change the Seed. That allows me to run the program several times and
// always obtain the same results.
#define PRIME 65539
#define SCALE 0.4656612875e-9
unsigned long Seed = 1234;

#define Rand()  (( Seed = ( (Seed * PRIME) & MAX_INT) ) * SCALE ) // return a random number
#define Randint(low,high) ( (int) (low + (high-(low)+1) * Rand())) // return an integer random number in the range [low-high]

// struct where I will save the coordinates of each city
struct point {
	double x;
	double y;
};

// variables for calculate the execution time
static clock_t start_time;
static double elapsed;

char problem_name[25];
int tam;

// read input file
// I only save the name of the problem, the size of the input and the coordinates of the cities
struct point * read_input(FILE *input) {
	struct point *data;
	char buffer[80];
	long int j;

	fscanf(input,"%s", buffer);
	while(strcmp("NODE_COORD_SECTION", buffer) != 0) {
		if (strcmp("NAME", buffer) == 0) {
			fscanf(input,"%s", buffer);
			fscanf(input,"%s", buffer);
			strcpy(problem_name, buffer);
		}
		else if (strcmp("DIMENSION", buffer) == 0) {
			fscanf(input,"%s", buffer);
			fscanf(input,"%d", &tam);	
		}
	
    	fscanf(input,"%s", buffer);
	}

	data = malloc(sizeof(struct point) * tam);
	for (int i = 0; i < tam; i++) {
		fscanf(input, "%ld %lf %lf", &j, &data[i].x, &data[i].y);
	}

    fprintf(stderr, "Problem name: %s \n", problem_name);
	fprintf(stderr, "Problem size: %d \n", tam);

	return data;
}

// create a distance matrix between all the cities; As input it has the array of structs containing
// the coordinates of all the cities.
long int ** distances(struct point * data) {
	long int **matrix;
	double x, y;

	matrix = malloc(sizeof(long int *) * tam);


	for (int i = 0; i < tam; i++) {
		matrix[i] = malloc(sizeof(long int) * tam);
		
		for (int j = 0; j < tam; j++) {
			x = data[i].x - data[j].x;
			y = data[i].y - data[j].y;
			matrix[i][j] = (long int) sqrt(x*x + y*y);
		}
	}

	return matrix;
}

// calculate the tour length of a given tour, fitness function
long int tour_length(int *tour, long int **dist) {
	long int length = 0;

	for (int i = 0; i < tam-1; i++) {
		length += dist[tour[i]][tour[i+1]];
	}
	length += dist[tour[0]][tour[tam-1]]; // sum the distance between the first and last cities of the tour

	return length;
}

// generates a random permutation (random tour)
// input tour must have a valid tour before start this function
void random_tour(int * tour) {
	int aux, x;
	for (int i = 0; i < tam; i++) {
		aux = tour[i];
		x = Randint(0, tam-1);
		tour[i] = tour[x];
		tour[x] = aux;
	}
}

// changes the cities x and y in the tour
void two_opt_swap(int *tour, int x, int y) {
	int aux = tour[x];
	tour[x] = tour[y];
	tour[y] = aux;
}

void nearest_neighbor(int * tour, long int **dist) {
	int *visited = malloc(sizeof(bool) * tam);
	int num_visited = 0, actual, min_dist = MAX_INT, nearest;

	for (int i = 0; i < tam; i++)
		visited[i] = FALSE;

	actual = Randint(0, tam-1);
	tour[0] = actual;
	visited[actual] = TRUE;
	num_visited++;

	while (num_visited < tam) {
		for (int j = 0; j < tam; j++){
			if (dist[actual][j] < min_dist) {
				if (!visited[j]) {
					min_dist = dist[actual][j];
					nearest = j;
				}
			}
		}
		tour[num_visited] = nearest;
		visited[nearest] = TRUE;
		actual = nearest;
		num_visited++;
		min_dist = MAX_INT;
	}

	free(visited);
}

// selection operator, Tournament selection
void selectionOp(int **population, int **selection, long int **dist) {

	for(int i = 0; i < POPULATION; i++) {
		int r1 = Randint(0, POPULATION-1);
		int r2 = Randint(0, POPULATION-1);

		int fitness1 = tour_length(population[r1], dist);
		int fitness2 = tour_length(population[r2], dist);

		if (fitness2 > fitness1)
			for(int j = 0; j < tam; j++)
				selection[i][j] = population[r1][j];
		else
			for(int j = 0; j < tam; j++)
				selection[i][j] = population[r2][j];

	}

}

// save in best[] the best tour of the population
void elite(int **population, int best[], long int **dist) {

	int actual_best = -1, best_fitness = tour_length(best, dist);

	for(int i = 0; i < POPULATION; i++) {
		int actual = tour_length(population[i], dist);
		if(actual < best_fitness) {
			actual_best = i;
			best_fitness = actual;
		}
	}

	if(actual_best != -1)
		for(int i = 0; i < tam; i++)
			best[i] = population[actual_best][i];

}

// mutation operator
void mutation(int **population) {

	for(int i = 0; i < POPULATION; i++) {
		for (int j = 0; j < tam; j++) {
			int r = Rand();

			if(r < PROB_MUTATION) {
				int r2 = Randint(0, tam-1);
				two_opt_swap(population[i], j, r2);
			}
		}
	}
}

// PMX crossover operator
void crossoverPMX(int **selection, int *child) {

	int found0 = 0, found1 = 0, gotostep = 1, posi, encontrado = 0, aux = 0, pos, di = 0;

	for(int i = 0; i < POPULATION-1; i++) {

		di = 0;
		for (int j = 0; j < tam; j++)
			if(selection[i][j] != selection[i+1][j])
				di = 1;

		if(Rand() < PROB_CROSSOVER && di) {

			for(int j = 0; j < tam; j++)
				child[j] = -1;

			int r1 = Randint(0, tam-1);
			int r2 = Randint(0, tam-1);

			//while(r2 <= r1)
			//	r2 = Randint(0, tam-1);
			if(r2 <= r1)
				r2 = r1+1;
			if(r2 >= tam)
				r2 = tam-1;

			for(int j = r1; j <= r2; j++) {
				child[j] = selection[i][j];
			}

			for(int j = r1; j <= r2; j++) {
				int buscando = j;

				while(!encontrado && !found0) {

					for(int z = 0; z < tam && !aux; z++) {
						if(selection[i+1][buscando] == child[z])
							found0 = 1;
					}

					if(found0 == 0){
						for(int z = 0; z < tam; z++){
							if(selection[i][buscando] == selection[i+1][z]){
								found1 = z;
							}
						}

						if(found1 < r1 || found1 > r2) {
								encontrado = 1;
								posi = found1;
						}
						else {
							buscando = found1;
							aux = 1;
						}
					}

				}

				encontrado = 0;
				if(found0 == 0)
					child[posi] = selection[i+1][j];
				found0 = 0;
				aux = 0;

			}
		
			for(int j = 0; j < tam; j++) {
				if(child[j] == -1)
					child[j] = selection[i+1][j];
			}

			for(int j = 0; j < tam; j++)
				selection[i][j] = child[j];

		}

	}

}

int main(int argc, char *argv[]) {
	//start_timer
	start_time = clock();

	//init_program
	FILE *input;
	struct point *data;
	long int **dist, length, best_length;
	int **population, **selection, *best_tour, *actual_best;
	int generation = 0;

	input = fopen(argv[1], "r");
	if ( input == NULL ) {
		fprintf(stderr, "No TSP data file specified, abort\n");
		exit(1);
    }
    else
    	fprintf(stderr, "Reading TSP data file: %s \n", argv[1]);

    data = read_input(input);
    fclose(input);

    population = malloc(sizeof(int *) * POPULATION);
	for (int i = 0; i < POPULATION; i++)
		population[i] = malloc(sizeof(int) * tam);

	selection = malloc(sizeof(int *) * POPULATION);
	for (int i = 0; i < POPULATION; i++)
		selection[i] = malloc(sizeof(int) * tam);
	
	best_tour = malloc(sizeof(int) * tam);

    actual_best = malloc(sizeof(int) * tam);
	int *child = malloc(sizeof(int) * tam);

    //compute_distances
    dist = distances(data);

    //Generate initial population
    for (int i = 0; i < POPULATION; i++)
    	for(int j = 0; j < tam; j++)
    		population[i][j] = j;

    for(int i = 0; i < POPULATION; i++)
    	random_tour(population[i]);

    for(int i = 0; i < POPULATION; i++)
    	nearest_neighbor(population[i], dist);

	for (int i = 0; i < tam; i++)
    	actual_best[i] = population[0][i];
    for (int i = 0; i < tam; i++)
    	best_tour[i] = population[0][i];
    elite(population, best_tour, dist);

	do
    {
        selectionOp(population, selection, dist);

        crossoverPMX(selection, child);
        mutation(selection);

        for(int i = 0; i < POPULATION; i++)
        	for(int j = 0; j < tam; j++)
        		population[i][j] = selection[i][j];

        elite(population, actual_best, dist);

        if(tour_length(actual_best, dist) < tour_length(best_tour, dist)) {
        	for(int j = 0; j < tam; j++) {
        		best_tour[j] = actual_best[j];
        	}
        }

        generation++;

    } while (generation < MAX_GENERATIONS);

	//stop_timer
	elapsed = clock() - start_time;
	elapsed = elapsed / CLOCKS_PER_SEC;
	fprintf(stderr, "Execution time: %f \n", elapsed);
	fprintf(stderr, "Best Tour: %ld \n", tour_length(best_tour, dist));

	fprintf(stderr, "Tour: ");
	for(int i = 0; i < tam; i++)
		fprintf(stderr, "%d ", best_tour[i]);
	fprintf(stderr, "\n");


	//free
	free(data);
	for (int i = 0; i < tam; i++)
		free(dist[i]);
	free(dist);
	for (int i = 0; i < POPULATION; i++)
		free(population[i]);
	free(population);
	free(best_tour);
	free(child);


	return(0);
}