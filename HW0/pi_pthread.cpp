#include <iostream>
#include <cstdlib>
#include <ctime>
#include <pthread.h>
#include <stdlib.h>
#define thread_num 100
using namespace std;

long long int number_in_circle[thread_num] = {0};
long long int number_of_tosses = 1e12;

void* norm(void* arg) {
    int thread_index = *((int*)arg);
    double x, y, ds;
    unsigned int local_seed = time(NULL);

    for (int toss = 0; toss < number_of_tosses / thread_num; toss++) {
        x = -1.0 + (rand_r(&local_seed) / (double)RAND_MAX) * 2.0;
        y = -1.0 + (rand_r(&local_seed) / (double)RAND_MAX) * 2.0;

        ds = x * x + y * y;
        if (ds <= 1)
            number_in_circle[thread_index] += 1;
    }
    pthread_exit((void*)0);
}

int main() {
    pthread_t threads[thread_num];
    int thread_indices[thread_num];
    
    for (int i = 0; i < thread_num; i++) {
        thread_indices[i] = i;
        pthread_create(&threads[i], NULL, norm, &thread_indices[i]);
    }

    for (int join_num = 0; join_num < thread_num; join_num++) {
        pthread_join(threads[join_num], NULL);
    }

    long long int total_in_circle = 0;
    for (int i = 0; i < thread_num; i++) {
        total_in_circle += number_in_circle[i];
    }

    double pi_estimate = 4 * static_cast<double>(total_in_circle) / number_of_tosses;
   
    cout.precision(10);
    cout << pi_estimate << endl;

    pthread_exit(NULL);
    return 0;
}
