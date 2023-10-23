#include <iostream>
#include <cstdlib>
#include <ctime>
#include <pthread.h>
#include <stdlib.h>
#define NUMOFTHREADS 2
using namespace std;

long long int number_in_circle1 = 0, number_in_circle2 = 0;
long long int number_of_tosses = 1e10;

// pthread_t callThd[NUMOFTHREADS];
// pthread_mutex_t mutexsum;

void* norm1(void*){
    double x1, y1, ds1;
    unsigned int local_seed = time(NULL);

    for (int toss1 = 0; toss1 < number_of_tosses/NUMOFTHREADS; toss1 ++)
    {
        x1 = -1.0 + (rand_r(&local_seed) / (double)RAND_MAX) * 2.0;
        y1 = -1.0 + (rand_r(&local_seed) / (double)RAND_MAX) * 2.0;

        ds1 = x1 * x1 + y1 * y1;
        if ( ds1 <= 1)
            number_in_circle1 += 1;
    }
    pthread_exit((void *)0);
}

void* norm2(void*){
    double x1, y1, ds1;
    unsigned int local_seed = time(NULL);

    for (int toss1 = 0; toss1 < number_of_tosses/NUMOFTHREADS; toss1 ++)
    {
        x1 = -1.0 + (rand_r(&local_seed) / (double)RAND_MAX) * 2.0;
        y1 = -1.0 + (rand_r(&local_seed) / (double)RAND_MAX) * 2.0;

        ds1 = x1 * x1 + y1 * y1;
        if ( ds1 <= 1)
            number_in_circle2 += 1;
    }
    pthread_exit((void *)0);
}

int main(){
    // pthread_mutex_init(&mutexsum, NULL);
    pthread_t t1, t2;
    // pthread_attr_t attr;
    // pthread_attr_init(&attr);
    // pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // for (int i = 0; i < NUMOFTHREADS; i++)
    // {
    //     pthread_create(&callThd[i], &attr, norm$i, NULL);
    // }

    pthread_create(&t1, NULL, norm1, NULL);
    pthread_create(&t2, NULL, norm2, NULL);
    pthread_join(t1, NULL);
    pthread_join(t2, NULL);
    // pthread_attr_destroy(&attr);

    // void *status;
    // for (int i = 0; i < NUMOFTHREADS; i++)
    // {
    //   // 等待每一個 thread 執行完畢
    //     pthread_join(t1, &status);
    // }

    double pi_estimate = 4 * static_cast<double>(number_in_circle1 + number_in_circle2) / number_of_tosses;
    printf("%.10lf \n", pi_estimate);
    pthread_exit(NULL);
    return 0;
}