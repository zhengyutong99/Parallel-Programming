#include <mpi.h>
#include <cstdio>
#include <cstdlib>

void construct_matrices(int *n_ptr, int *m_ptr, int *l_ptr, int **a_mat_ptr, int **b_mat_ptr) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        // Rank 0 reads the dimensions and data
        scanf("%d %d %d", n_ptr, m_ptr, l_ptr);
    }

    // Broadcast dimensions to all processes
    MPI_Bcast(n_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(m_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(l_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate memory in all processes
    *a_mat_ptr = (int *)malloc(sizeof(int) * (*n_ptr) * (*m_ptr));
    *b_mat_ptr = (int *)malloc(sizeof(int) * (*m_ptr) * (*l_ptr));

    if (rank == 0) {
        // Rank 0 reads matrix data
        for (int i = 0; i < (*n_ptr) * (*m_ptr); i++) {
            scanf("%d", *a_mat_ptr + i);
        }
        for (int i = 0; i < (*m_ptr) * (*l_ptr); i++) {
            scanf("%d", *b_mat_ptr + i);
        }
    }

    // Broadcast matrix data to all processes
    MPI_Bcast(*a_mat_ptr, (*n_ptr) * (*m_ptr), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(*b_mat_ptr, (*m_ptr) * (*l_ptr), MPI_INT, 0, MPI_COMM_WORLD);
    // if (rank == 0) {
    //     scanf("%d %d %d", n_ptr, m_ptr, l_ptr);
    //     *a_mat_ptr = (int *)malloc(sizeof(int) * (*n_ptr) * (*m_ptr));
    //     *b_mat_ptr = (int *)malloc(sizeof(int) * (*m_ptr) * (*l_ptr));
    //     for (int i = 0; i < (*n_ptr) * (*m_ptr); i++) {
    //         scanf("%d", *a_mat_ptr + i);
    //     }
    //     for (int i = 0; i < (*m_ptr) * (*l_ptr); i++) {
    //         scanf("%d", *b_mat_ptr + i);
    //     }
    // }
}

void matrix_multiply(const int n, const int m, const int l, const int *a_mat, const int *b_mat) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int *c_mat = (int *)malloc(sizeof(int) * n * l);
    if (c_mat == NULL) {
        fprintf(stderr, "Unable to allocate memory for matrix C.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Assuming a simple parallelization where each process computes the entire matrix multiplication
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            c_mat[i * l + j] = 0;
            for (int k = 0; k < m; k++) {
                c_mat[i * l + j] += a_mat[i * m + k] * b_mat[k * l + j];
            }
        }
    }

    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                printf("%d ", c_mat[i * l + j]);
            }
            printf("\n");
        }
    }

    free(c_mat);
}

void destruct_matrices(int *a_mat, int *b_mat) {
    free(a_mat);
    free(b_mat);
}
