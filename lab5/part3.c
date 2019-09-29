#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "my_rand.h"

void Usage(char* prog_name);
int toss(int tosses);

int main(int argc, char* argv[]) {
  double global_result = 0.0;
  int thread_count, number_of_tosses;
  if (argc != 2)
    Usage(argv[0]);
  //thread_count = strtol(argv[1], NULL, 10);
  number_of_tosses = strtol(argv[1], NULL, 10);
  thread_count = omp_get_max_threads();
  printf("num threads:%d\n", thread_count);
  if (number_of_tosses % thread_count != 0) Usage(argv[0]);
  //printf("num threads:%d\n", thread_count);
# pragma omp parallel num_threads(thread_count) reduction(+: global_result)
  global_result += toss(number_of_tosses);
  double pi_estimate = 4 * global_result/((double)number_of_tosses);
  printf("%f\n", pi_estimate);
}

void Usage(char* prog_name) {
  fprintf(stderr, "usage: %s <number of threads>\n", prog_name);
  exit(0);
}
  
int toss(int tosses) {
  int i; 
  unsigned seed = (unsigned) time(NULL), p, q;
  double e1=1, e2=1, y;
  int toss;
  long number_in_circle = 0;
  int my_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();

  e1 = my_rand(&seed);
  e2 = my_rand(&seed);
  for(i = 0; i < 10 * my_rank; i++) {
    if (i % my_rank == 0) {
      e1 = my_drand(&e1);
      e2 = my_drand(&e2);
    } else {
      e1 = my_rand(&e1);
      e2 = my_rand(&e2);
    }
  }
  int n = tosses / thread_count;
  //printf("%d %d\n", my_rank, n);
  for(toss = 0; toss < n; toss++)
  {
    e1 = my_drand(&e1);
    e2 = my_drand(&e2);
    if (e1 < 0) {
      e1 *= -1;
    }
    if (e2 < 0) {
      e2 *= -1;
    }
    e1 *= 2;
    e2 *= 2;
    e1 -= 1;
    e2 -= 1;
    // printf("(%f, %f)\n", e1, e2);
    double distance_squared = e1 * e1 + e2 * e2;
    if(distance_squared <= 1) number_in_circle++;
  }
  printf("thread-%d: number in circle: %d\n", my_rank, number_in_circle);
  return number_in_circle;
}
