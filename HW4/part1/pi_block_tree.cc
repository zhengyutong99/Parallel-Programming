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
    MPI_Status  status;

    unsigned int seed = (unsigned)time(NULL) + (unsigned)getpid();
    long long int local_tosses = tosses / world_size;
    long long int local_in_circle = 0;

    for (long long int toss = 0; toss < local_tosses; toss++) {
        double x = (double)rand_r(&seed) / RAND_MAX * 2 - 1;
        double y = (double)rand_r(&seed) / RAND_MAX * 2 - 1;
        if (x * x + y * y <= 1) {
            local_in_circle++;
        }
    }

    // TODO: binary tree reduction
    long long int global_in_circle = 0;

    for (int step = 1; step < world_size; step *= 2) {
        if (world_rank % (2 * step) == 0) {
            if (world_rank + step < world_size) {
                long long int other_in_circle;
                MPI_Recv(&other_in_circle, 1, MPI_LONG_LONG_INT, world_rank + step, 0, MPI_COMM_WORLD, &status);
                local_in_circle += other_in_circle;
            }
        } else {
            MPI_Send(&local_in_circle, 1, MPI_LONG_LONG_INT, world_rank - step, 0, MPI_COMM_WORLD);
            break;
        }
    }

    if (world_rank == 0)
    {
        // TODO: PI result
        global_in_circle = local_in_circle;
        pi_result = 4 * ((double)global_in_circle / (double)tosses);
        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
