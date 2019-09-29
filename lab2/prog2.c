# include <stdio.h>
# include <math.h>


float calc_square(float n) {
  float xNext;
  float x = 1.0;
  int count = 0;
  while (count < 50) {
    float help = n/x;
    xNext = 0.5 * (x + help);
    // printf("i=%d:help:%f,xNext:%f", count, help, xNext);
    x = xNext;
    count++;
  }
  return xNext;
}


int main() {
  float n;
  printf("Enter a number:");
  scanf("%f", &n);
  float newt = calc_square(n);
  printf("newton square: %.4f, math.h square: %f\n", newt, sqrt(n));
}
