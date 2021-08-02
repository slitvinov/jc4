/*
  c99  2term.c `pkg-config --libs --cflags gsl`
*/

#include "util.inc"
enum { LEFT, RIGHT };
enum { NS = 5, NA = 2 };
const double P[NA][NS][NS] = {
    {
        {0.0, 0.0, 0.0, 0.0, 1.0},
        {0.8, 0.1, 0.1, 0.0, 0.0},
        {0.0, 0.8, 0.1, 0.1, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 1.0},
    },
    {
        {0.0, 0.0, 0.0, 0.0, 1.0},
        {0.1, 0.1, 0.8, 0.0, 0.0},
        {0.0, 0.1, 0.1, 0.8, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 1.0},
    },
};
const int A[NA] = {LEFT, RIGHT};
const char *Astr[NA] = {"<", ">"};
const double g = 4.0 / 5.0;
const double R[NS] = {0, 0, 0, 1, 0};

void value(const double *R, const int *p, double *V) {
  double G[NS * NS];
  double Ginv[NS * NS];
  int i, j;
  for (i = 0; i < NS; i++)
    for (j = 0; j < NS; j++)
      G[NS * i + j] = (i == j) - g * P[p[i]][i][j];
  inv(NS, G, Ginv);
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

void get_reward(double alpha, double *R) {
  int i;
  for (i = 0; i < NS - 1; i++)
    R[i] = alpha * i;
  R[NS - 1] = 0;
}

void forward(const double *R, int *p) {
  int t;
  double V[NS];
  for (t = 0; t < 20; t++) {
    value(R, p, /**/ V);
    policy(R, V, /**/ p);
  }
}

int main() {
  int i, t;
  double Rphi[NS], Vp[NS], Vm[NS];
  int p[NS] = {LEFT, LEFT, LEFT, LEFT, LEFT};
  double phi = 0.0;
  double eps = 0.01;
  double grad;

  for (t = 0; t < 10; t++) {
    get_reward(phi + eps, Rphi);
    forward(Rphi, p);
    value(R, p, Vp);
    get_reward(phi - eps, Rphi);
    forward(Rphi, p);
    value(R, p, Vm);
    grad = 0;
    for (i = 0; i < NS - 1; i++)
      grad += Vp[i] - Vm[i];
    phi += grad * eps;
    printf("%g\n", grad);
  }
  for (i = 0; i < NS - 1; i++)
    printf("%6.2f%s ", Vp[i], Astr[p[i]]);
  printf("\n");
}
