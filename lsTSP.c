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
long int tam;

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
			fscanf(input,"%ld", &tam);	
		}
	
    	fscanf(input,"%s", buffer);
	}

	data = malloc(sizeof(struct point) * tam);
	for (int i = 0; i < tam; i++) {
		fscanf(input, "%ld %lf %lf", &j, &data[i].x, &data[i].y);
	}

    fprintf(stderr, "Problem name: %s \n", problem_name);
	fprintf(stderr, "Problem size: %ld \n", tam);

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

// calculate the tour length of a given tour
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

// nearest neighbor heuristic
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

// changes the cities x and y in the tour
void two_opt_swap(int *tour, int x, int y) {
	int aux = tour[x];
	tour[x] = tour[y];
	tour[y] = aux;
}

// factorization of the lenght of the tour for be able to perform it faster
long int two_opt_length(int * tour, int x, int y, long int best_length, long int **dist) {
	long int actual_length = best_length;
	int x0 = x, y0 = y, x1 = x, y1 = y;

	if (x == 0)
		x0 = tam;
	if (x == tam - 1)
		x1 = -1;
	if (y == 0)
		y0 = tam;
	if (y == tam - 1)
		y1 = -1;

	int a = dist[tour[x0-1]][tour[x]] + dist[tour[x]][tour[x1+1]] + dist[tour[y0-1]][tour[y]] + dist[tour[y]][tour[y1+1]];
	int b = dist[tour[x0-1]][tour[y]] + dist[tour[y]][tour[x1+1]] + dist[tour[y0-1]][tour[x]] + dist[tour[x]][tour[y1+1]];

	actual_length = actual_length + a - b;

	return actual_length;
}

// 2-opt local search
int * two_opt_ls(int * tour, long int **dist) {
	int * best_tour = malloc(sizeof(int) * tam);
	long int best_length = tour_length(tour, dist);
	long int actual_length;
	bool change;

	int cont = 0;

	for (int i = 0; i < tam; i++)
		best_tour[i] = tour[i];

	do {
		change = FALSE;

		cont++;

		for (int i = 0; i < tam && !change; i++) {
			for (int j = i+1; j < tam && !change; j++) {
				two_opt_swap(tour, i, j);
				actual_length = two_opt_length(tour, i, j, best_length, dist);

				if (actual_length < best_length) {
					best_length = actual_length;

					//printf("New best: %ld \n", best_length);

					for(int a = 0; a < tam; a++)
						best_tour[a] = tour[a];
					
					change = TRUE;
				}
				else 
					two_opt_swap(tour, i, j);
			}
		}
	} while(change && cont < 40);

	return best_tour;
}

// function that calculates the k nearest neighbours of each city
void compute_nearest_k_neighbors(int ** NN, long int ** dist, int k) {
	long int prev = 0, actual = MAX_INT;
	int nearest;

	for(int i = 0; i < tam; i++) {
		for(int j = 0; j < k; j++) {
			for (int z = 0; z < tam; z++) {
				
				if(dist[i][z] < actual && dist[i][z] > prev) {
					actual = dist[i][z];
					nearest = z;
				}
			}
			prev = actual;
			NN[i][j] = nearest;
			actual = MAX_INT;
		}
		prev = 0;
	}

	/*for(int i = 0; i < tam; i++) {
    	for (int j = 0; j < k; j++)
    		printf("%d  ", NN[i][j]);
    	printf("\n");
    }*/
}

// 2-opt local search but only considers for change the k nearest neighbours of each city
int * two_opt_ls_k_neighbor(int * tour, long int **dist, int k) {
	int * best_tour = malloc(sizeof(int) * tam);
	long int best_length = tour_length(tour, dist);
	long int actual_length;
	bool change, found = FALSE;
	int elem;

	int **NN = malloc(sizeof(int *) * tam);

	for(int i = 0; i < tam; i++)
		NN[i] = malloc(sizeof(int) * k);

	compute_nearest_k_neighbors(NN, dist, k);

	int cont = 0;

	for (int i = 0; i < tam; i++)
		best_tour[i] = tour[i];

	do {
		change = FALSE;

		cont++;

		for (int i = 0; i < tam-1 && !change; i++) {
			for (int j = 0; j < k && !change; j++) {

				for(int b = 0; b < tam && !found; b++) {
					if (tour[b] == NN[tour[i]][j]) {
						found = TRUE;
						elem = b;
					}
				}

				found = FALSE;

				two_opt_swap(tour, i+1, elem);
				actual_length = two_opt_length(tour, i+1, elem, best_length, dist);

				if (actual_length < best_length) {
					best_length = actual_length;

					//printf("New best: %ld \n", best_length);

					for(int a = 0; a < tam; a++)
						best_tour[a] = tour[a];
					
					change = TRUE;
				}
				else 
					two_opt_swap(tour, i+1, elem);
			}
		}
	} while(change && cont < 40);

	for (int i = 0; i < tam; i++)
		free(NN[i]);
	free(NN);

	return best_tour;
}



int main(int argc, char *argv[]) {
	//start_timer

	//init_program
	FILE *input, *output;
	struct point *data;
	long int **dist, length, best_length;
	int *tour, *best_tour;
	int k = 10;

	output = fopen("output.txt", "w");
	input = fopen(argv[1], "r");
	if ( input == NULL ) {
		fprintf(stderr, "No TSP data file specified, abort\n");
		exit(1);
    }
    else
    	fprintf(stderr, "Reading TSP data file: %s \n", argv[1]);

    data = read_input(input);
    fclose(input);

    tour = malloc(sizeof(int) * tam);

    //compute_distances
    dist = distances(data);

    /*for(int i = 0; i < tam; i++) {
    	for (int j = 0; j < tam; j++)
    		printf("%ld  ", dist[i][j]);
    	printf("\n");
    }*/

 	start_time = clock();

	//compute_tour_ini
	for (int i = 0; i < tam; i++)
		tour[i] = i;
	length = tour_length(tour, dist);
    fprintf(stderr, "Initial tour ([0,1,2,...,n]) length: %ld \n", length);

    elapsed = clock() - start_time;
	elapsed = elapsed / CLOCKS_PER_SEC;
	fprintf(stderr, "Execution time: %f \n", elapsed);
	start_time = clock();

	//compute_tour_random
    random_tour(tour);
    length = tour_length(tour, dist);
    fprintf(stderr, "Random tour length: %ld \n", length);

    elapsed = clock() - start_time;
	elapsed = elapsed / CLOCKS_PER_SEC;
	fprintf(stderr, "Execution time: %f \n", elapsed);
 	start_time = clock();

	//compute_tour_nearest_neighbor
    nearest_neighbor(tour, dist);
    length = tour_length(tour, dist);
    fprintf(stderr, "Nearest neighbor tour length: %ld \n", length);

    elapsed = clock() - start_time;
	elapsed = elapsed / CLOCKS_PER_SEC;
	fprintf(stderr, "Execution time: %f \n", elapsed);
 	start_time = clock();

    /*random_tour(tour);
	best_tour = two_opt_ls(tour, dist);
    best_length = tour_length(best_tour, dist);
    fprintf(stderr, "2-opt local search tour length: %ld \n", best_length);*/

    //2-opt local search, consider only the k nearest neighbor space in each iteration
    //random_tour(tour);
    best_tour = two_opt_ls_k_neighbor(tour, dist, k);
	best_length = tour_length(best_tour, dist);
    fprintf(stderr, "2-opt local search tour length: %ld \n", best_length);
    

    fprintf(output, "Best tour length: %ld \n\nPermutation:\n", best_length);
    for (int i = 0; i < tam; i++) {
    	fprintf(output, "%d\n", best_tour[i]);
    }
    fclose(output);

	//stop_timer
	elapsed = clock() - start_time;
	elapsed = elapsed / CLOCKS_PER_SEC;
	fprintf(stderr, "Execution time: %f \n", elapsed);

	//free
	free(data);
	for (int i = 0; i < tam; i++)
		free(dist[i]);
	free(dist);
	free(tour);
	free(best_tour);

	return(0);
}