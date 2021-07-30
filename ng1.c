/*
c99  ng1.c `pkg-config --libs --cflags gsl`
*/

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

enum {LEFT, STAY, RIGHT};
enum
{ NS = 4, NA = 3 };

double P[NA][NS][NS] = {
  {
   {1, 0, 0, 0},
   {1, 0, 0, 0},
   {0, 1, 0, 0},
   {0, 0, 1, 0},
   },
  {
   {1, 0, 0, 0},
   {0, 1, 0, 0},
   {0, 0, 1, 0},
   {0, 0, 0, 1},
   },
  {
   {0, 1, 0, 0},
   {0, 0, 1, 0},
   {0, 0, 0, 1},
   {0, 0, 0, 1},
   },
};

double R[NS] = {
  1, 0, 0, 1,
};

double g = 4/5.0;

int
inv(double *a, double *b)
{
  int s;
  gsl_matrix_view A;
  gsl_matrix_view B;

  A = gsl_matrix_view_array(a, NS, NS);
  B = gsl_matrix_view_array(b, NS, NS);
  
  gsl_permutation *p = gsl_permutation_alloc(NS);
  gsl_linalg_LU_decomp(&A.matrix, p, &s);
  gsl_linalg_LU_invert(&A.matrix, p, &B.matrix);
  
  gsl_permutation_free(p);
}

int
main ()
{
  int i, j;
  double G[NS*NS];
  double Ginv[NS*NS];
  double V[NS];
  int p[NS] = {RIGHT, STAY, STAY, LEFT};

  for (i = 0; i < NS; i++)
    for (j = 0; j < NS; j++)
      G[NS * i + j] = (i == j) - g * P[p[i]][i][j];
  inv(G, Ginv);

  for (i = 0; i < NS; i++) {
    V[i] = 0;
    for (j = 0; j < NS; j++)
      V[i] += Ginv[NS * i + j] * R[j];
  }

  for (i = 0; i < NS; i++)
    printf("%6.2g ", V[i]);
  printf("\n");
}
