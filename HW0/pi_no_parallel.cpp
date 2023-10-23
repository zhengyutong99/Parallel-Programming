#include <iostream>
#include <cstdlib>
#include <ctime>
#include <stdlib.h>
using namespace std;

long long int number_in_circle = 0;
long long int number_of_tosses = 1e8;

int main(){
    double x, y, distance_squared;
    srand((unsigned) time(NULL));
    
    for (int toss = 0; toss < number_of_tosses; toss ++) {
        x = -1.0 + (rand() / (double)RAND_MAX) * 2.0;
        y = -1.0 + (rand() / (double)RAND_MAX) * 2.0;

        distance_squared = x * x + y * y;
        if ( distance_squared <= 1)
            number_in_circle++;
    }
    double pi_estimate = 4 * static_cast<double>(number_in_circle) / number_of_tosses;

    printf("%.10lf \n", pi_estimate);
    return 0;
}