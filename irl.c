/*
c99  irl.c `pkg-config --libs --cflags gsl`
*/

#include "util.inc"
#include <gsl/gsl_rng.h>

enum { LEFT, STAY, RIGHT };
enum { NS = 5, NA = 3, NP = 2, NT = 100, MAX_LEN = 100 };
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
	{0.1, 0.8, 0.1, 0.0, 0.0},
	{0.0, 0.1, 0.8, 0.1, 0.0},
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
const int A[NA] = {LEFT, STAY, RIGHT};
const char *Astr[NA] = {"<", "=", ">"};
const double g = 9.0 / 10.0;

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

void g_inv(const int *p, double *Ginv) {
  double G[NS * NS];
  int i, j;
  for (i = 0; i < NS; i++)
    for (j = 0; j < NS; j++)
      G[NS * i + j] = (i == j) - g * P[p[i]][i][j];
  inv(NS, G, Ginv);
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

void forward(const double *R, int *p) {
  int t;
  double V[NS];
  for (t = 0; t < 20; t++) {
    value(R, p, /**/ V);
    policy(R, V, /**/ p);
  }
}

double score(const int *p, const double *Ginv, const double *R, int *inside) {
  int a, i, j, k;
  double delta, delta_min, ans, eps;
  eps = 1e-8;
  ans = 0;
  for (i = 0; i < NS; i++) {
    delta_min = DBL_MAX;
    for (a = 0; a < NA; a++)
      if (p[a] != A[a]) {
	delta = 0;
	for (j = 0; j < NS; j++)
	  for (k = 0; k < NS; k++)
	    delta += (P[p[i]][i][j] - P[A[a]][i][j]) * Ginv[j * NS  + k] * R[k];
	if (delta < delta_min)
	  delta_min = delta;
      }
    if (delta_min < -eps) {
      *inside = 0;
      return DBL_MAX;
    }
    ans += delta_min;
  }
  *inside = 1;
  return ans;
}

double vabs(int n, double *a) {
  double ans;
  int i;
  ans = 0;
  for (i = 0; i < n; i++)
    ans += fabs(a[i]);
  return ans/n;
}

void copy(int n, const double *a, double *b)
{
  int i;
  for (i = 0; i < n; i++)
    b[i] = a[i];
}

int main() {
  double Ginv[NS * NS], sc, sc_max, step, lambda;
  int inside;
  int p[NS] = {LEFT,    LEFT, RIGHT,    STAY,     STAY};
  double R[NS], R_max[NS];
  g_inv(p, Ginv);

  R[NS - 1] = 0;
  sc_max = -DBL_MAX;
  step = 0.05;
  lambda = 0.0;
  for (R[0] = -1; R[0] <= 1; R[0] += step)  
    for (R[1] = -1; R[1] <= 1; R[1] += step)
      for (R[2] = -1; R[2] <= 1; R[2] += step)
	for (R[3] = -1; R[3] <= 1; R[3] += step) {
	  sc = score(p, Ginv, R, &inside) - lambda * vabs(NS, R);
	  if (inside && sc > sc_max) {
	    sc_max = sc;
	    copy(NS, R, R_max);
	  }
	}
  fprintf(stderr, "%.2f %.2f %.2f %.2f [%.2f]\n", R_max[0], R_max[1], R_max[2], R_max[3], sc_max);
}
