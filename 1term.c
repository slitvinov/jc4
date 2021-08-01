/*
  c99  1term.c `pkg-config --libs --cflags gsl`
*/

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>

enum { LEFT, RIGHT };
enum { NS = 5, NA = 2 };
const double P[NA][NS][NS] = {
    {
        {0, 0, 0, 0, 1},
        {1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0},
        {0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1},
    },
    {
        {0, 0, 0, 0, 1},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1},
        {0, 0, 0, 0, 1},
    },
};
const int A[NA] = {LEFT, RIGHT};
const char *Astr[NA] = {"<", ">"};
const double g = 4.0 / 5.0;

int inv(double *a, double *b) {
  int s;
  gsl_matrix_view A;
  gsl_matrix_view B;
  gsl_permutation *p;

  A = gsl_matrix_view_array(a, NS, NS);
  B = gsl_matrix_view_array(b, NS, NS);
  p = gsl_permutation_alloc(NS);
  gsl_linalg_LU_decomp(&A.matrix, p, &s);
  gsl_linalg_LU_invert(&A.matrix, p, &B.matrix);
  gsl_permutation_free(p);
}

void value(const double *R, const int *p, double *V) {
  double G[NS * NS];
  double Ginv[NS * NS];
  int i, j;
  for (i = 0; i < NS; i++)
    for (j = 0; j < NS; j++)
      G[NS * i + j] = (i == j) - g * P[p[i]][i][j];
  inv(G, Ginv);
  for (i = 0; i < NS; i++) {
    V[i] = 0;
    for (j = 0; j < NS; j++)
      V[i] += Ginv[NS * i + j] * R[j];
  }
}

void policy(const double *R, const double *V, int *p) {
  double Q, Qmax, amax;
  int a, j, i;
  for (i = 0; i < NS; i++) {
    Qmax = -DBL_MAX;
    for (a = 0; a < NA; a++) {
      Q = R[i];
      for (j = 0; j < NS; j++)
        Q += g * P[A[a]][i][j] * V[j];
      if (Q > Qmax) {
        Qmax = Q;
        amax = A[a];
      }
    }
    p[i] = amax;
  }
}

int main() {
  int i, j, t;
  double V[NS];
  int p[NS] = {LEFT, LEFT, LEFT, LEFT, LEFT};
  double R[NS] = {0, 0, 0, 1, 0};
  for (t = 0; t < 10; t++) {
    value(R, p, /**/ V);
    policy(R, V, p);
    for (i = 0; i < NS - 1; i++)
      printf("%6.2f%s ", V[i], Astr[p[i]]);
    printf("\n");
  }
}
