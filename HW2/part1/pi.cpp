#include <iostream>
#include <cstdlib>
#include <ctime>
#include <pthread.h>
#include <stdlib.h>
using namespace std;

int thread_num;
long long int number_in_circle;
long long int number_of_tosses;
pthread_mutex_t mutex;

void* norm(void* arg) {
    // int thread_index = *((int*)arg);
    double x, y, ds;
    unsigned int local_seed = time(NULL);
    long long int tosses_per_thread = number_of_tosses / thread_num;
    long long int result = 0;

    for (int toss = 0; toss < tosses_per_thread; toss++) {
        x = -1.0 + (rand_r(&local_seed) / (double)RAND_MAX) * 2.0;
        y = -1.0 + (rand_r(&local_seed) / (double)RAND_MAX) * 2.0;

        ds = x * x + y * y;

        if (ds <= 1)
            result ++;
    }
    pthread_mutex_lock(&mutex);
    number_in_circle += result;
    pthread_mutex_unlock(&mutex);

    pthread_exit((void*)0);
}

int main(int argc, char* argv[]) {
    thread_num = atoi(argv[1]);
    number_of_tosses = atoll(argv[2]);

    pthread_mutex_init(&mutex, NULL);
    pthread_t threads[thread_num];
    int thread_indices[thread_num];
    // number_in_circle = new long long int[thread_num];
    
    for (int i = 0; i < thread_num; i++) {
        thread_indices[i] = i;
        pthread_create(&threads[i], NULL, norm, &thread_indices[i]);
    }

    // pthread_mutex_lock(&mutex);
    // long long int total_in_circle = 0;
    // for (int i = 0; i < thread_num; i++) {
    //     total_in_circle += number_in_circle[i];
    // }
    // pthread_mutex_unlock(&mutex);

    for (int join_num = 0; join_num < thread_num; join_num++) {
        pthread_join(threads[join_num], NULL);
    }

    if (number_of_tosses != 0) {
        // double pi_estimate = 4 * static_cast<double>(total_in_circle) / number_of_tosses;
        double pi_estimate = 4 * static_cast<double>(number_in_circle) / number_of_tosses;
        cout.precision(10);
        cout << pi_estimate << endl;
    } else {
        cout << "Error: number_of_tosses is zero." << endl;
    }

    pthread_mutex_destroy(&mutex);
    pthread_exit(NULL);
    return 0;
}
