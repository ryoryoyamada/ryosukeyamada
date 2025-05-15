#include <stdio.h>

int main(void) {
  char molecule[6] = {'A', 'B', 'C', 'D', 'E', 'F'}; /* Name of molecule */
  int matrix[6][6] = {{-1, 0, 0, 0, 0, 0},           /* row of A */
                      {1, *, *, *, *, *},            /* row of B */
                      {0, *, *, *, *, *},            /* row of C */
                      {0, *, *, *, *, *},            /* row of D */
                      {0, *, *, *, *, *},            /* row of E */
                      {0, *, *, *, *, *}};           /* row of F */

  /* 以降、自分のプログラムを記述 */

  return 0;
}
