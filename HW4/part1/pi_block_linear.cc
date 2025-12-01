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

    // TODO: init MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    long long int number_of_tosses = tosses / world_size;
    long long int number_in_circle = 0;
    double x, y, ds;
    unsigned int seed = (unsigned)time(NULL) + (unsigned)world_rank;
    MPI_Status  status;

    for (long long int toss = 0; toss < number_of_tosses; toss++) {
        x = (double)rand_r(&seed) / RAND_MAX * 2.0 - 1.0;
        y = (double)rand_r(&seed) / RAND_MAX * 2.0 - 1.0;
        ds = x * x + y * y;
        if (ds <= 1) number_in_circle++;
    }

    if (world_rank > 0)
    {
        // TODO: handle workers
        MPI_Send(&number_in_circle, 1, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
    }
    else if (world_rank == 0)
    {
        // TODO: master
        long long int temp_count;
        for (int source = 1; source < world_size; source++) {
            MPI_Recv(&temp_count, 1, MPI_LONG_LONG_INT, source, 0, MPI_COMM_WORLD, &status);
            number_in_circle += temp_count;
        }
    }

    if (world_rank == 0)
    {
        // TODO: process PI result
        pi_result = 4 * number_in_circle / ((double)tosses);
        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
