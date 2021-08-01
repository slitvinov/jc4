/*
c99  1.c `pkg-config --libs --cflags gsl`
*/

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include "util.inc"

enum { LEFT, STAY, RIGHT };
enum { NS = 4, NA = 3 };

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

int A[NA] = {LEFT, STAY, RIGHT};
const char *Astr[NA] = {"<", "o", ">"};

double R[NS] = {
    1,
    0,
    0,
    1,
};

double g = 4 / 5.0;

int main() {
  int i, j, t;
  double G[NS * NS];
  double Ginv[NS * NS];
  double V[NS];
  int p[NS] = {RIGHT, STAY, STAY, LEFT};

  for (t = 0; t < 10; t++) {

    for (i = 0; i < NS; i++)
      for (j = 0; j < NS; j++)
        G[NS * i + j] = (i == j) - g * P[p[i]][i][j];
    inv(NS, G, Ginv);

    for (i = 0; i < NS; i++) {
      V[i] = 0;
      for (j = 0; j < NS; j++)
        V[i] += Ginv[NS * i + j] * R[j];
    }

    double Q, Qmax, amax;
    int a;

    for (i = 0; i < NS; i++) {
      Qmax = -DBL_MAX;
      for (a = 0; a < NA; a++) {
        Q = R[i];
        for (j = 0; j < NS; j++) {
          Q += g * P[A[a]][i][j] * V[j];
        }
        if (Q > Qmax) {
          Qmax = Q;
          amax = A[a];
        }
      }
      p[i] = amax;
    }

    for (i = 0; i < NS; i++)
      printf("%6.2g:%s ", V[i], Astr[p[i]]);
    printf("\n");
  }
}
