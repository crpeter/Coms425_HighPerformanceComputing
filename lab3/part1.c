# include <stdio.h>
# include <stdlib.h>

int main() {
  srand(time(0));
  int rollHelp = 7;
  int i;
  int results[11];
  for (i = 0; i < 11; i++) {
    results[i] = 0;
  }
  for (i = 0; i < 36000; i++) {
    int die1 = rand() % rollHelp;
    int die2 = rand() % rollHelp;
    int sum = die1 + die2;
    results[sum-1]++;
  }
  for (i = 0; i < 11; i++) {
    printf("sum %2d: %5d percent: %f\n", i+2, results[i],\
	   ((double)results[i] / 36000.0));
  }
  printf("END\n");
  //printf("first: %d, second: %d\n", die1, die2);
  return 0;
}
