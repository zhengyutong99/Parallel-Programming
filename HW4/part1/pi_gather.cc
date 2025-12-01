#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---

    // TODO: MPI init
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    long long int number_of_tosses = tosses / world_size;
    long long int number_in_circle = 0;
    double x, y, distance_squared;
    unsigned int seed = (unsigned)time(NULL) + (unsigned)world_rank;

    for (long long int toss = 0; toss < number_of_tosses; toss++) {
        x = (double)rand_r(&seed) / RAND_MAX * 2.0 - 1.0;
        y = (double)rand_r(&seed) / RAND_MAX * 2.0 - 1.0;
        distance_squared = x * x + y * y;
        if (distance_squared <= 1) number_in_circle++;
    }

    long long int *gathered_counts = NULL;
    if (world_rank == 0) {
        gathered_counts = (long long int *)malloc(sizeof(long long int) * world_size);
    }

    MPI_Gather(&number_in_circle, 1, MPI_LONG_LONG_INT, gathered_counts, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

    // TODO: use MPI_Gather
    if (world_rank == 0)
    {
        // TODO: PI result
        long long int total_in_circle = 0;
        for (int i = 0; i < world_size; i++) {
            total_in_circle += gathered_counts[i];
        }

        pi_result = 4 * total_in_circle / ((double)tosses);

        free(gathered_counts);

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }
    
    MPI_Finalize();
    return 0;
}
