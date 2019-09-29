#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main(int argc, char* argv[]) {
  int     n, i;
  int     thread_count;

  printf("Enter iterations, num threads\n");
  scanf("%d %d", &n, &thread_count);
  //thread_count = omp_get_max_threads();
  printf("num threads: %d\n", thread_count);
  for (i = 0; i < thread_count; i++) {
    int local_a = (n / thread_count)* i;
    int local_b = (n / thread_count)* i + ((n / thread_count) - 1);
    int r = n % thread_count;
    if (r != 0) {
      if (i < r) {
	local_a += i;
	local_b += i + 1;
      } else {
	local_a += r;
	local_b += r;
      }
    }
    printf("Thread %2d: Iterations %3d - %3d\n", i, local_a, local_b);
  }
  return 0;
}
