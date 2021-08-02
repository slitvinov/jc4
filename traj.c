/*
  c99  traj.c `pkg-config --libs --cflags gsl`
*/

#include "util.inc"
#include <gsl/gsl_rng.h>

enum { LEFT, STAY, RIGHT };
enum { NS = 5, NA = 3, NT = 10, MAX_LEN = 100 };
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
const double g = 4.0 / 5.0;

int traj[NT][MAX_LEN];
double Ret[NT];
double Prob[NT];

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

int rnd(gsl_rng *r, int lo, int hi)
{
  return lo + gsl_rng_uniform_int(r, hi - lo);
}

void gen_traj(void) {
  int t, s, s0, j, a;
  const gsl_rng_type * T;
  gsl_rng * r;

 T = gsl_rng_default;
 r = gsl_rng_alloc (T);

 for (t = 0; t < NT; t++) {
   s = rnd(r, 0, NS - 1);
   for (j = 0; ; j++) {
     traj[t][j] = s;
     fprintf(stderr, "s: %d:%d\n", t, s);
     if (s == NS - 1 || j == MAX_LEN - 1)
       break;
     do {
       a = rnd(r, 0, NA);
       s0 = rnd(r, 0, NS);
     } while (P[A[a]][s][s0] == 0.0);
     s = s0;
   }
 }
 gsl_rng_free(r);
}

void traj_ret(const double *R) {
  double ret;
  int t, j;
  for (t = 0; t < NT; t++) {
    ret = 0;
    for (j = 0; j < MAX_LEN && traj[t][j] != NS - 1; j++) {
      ret += R[traj[t][j]];
    }
    Ret[t] = ret;
    fprintf(stderr, "r: %d:%g\n", t, ret);
  }
}

void traj_prob(const int *p) {
  
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

int main() {
  int i;
  double V[NS];
  int p[NS] = {STAY, STAY, STAY, STAY, STAY};
  double R[NS] = {0, 0, 0, 1, 0};
  
  gen_traj();
  traj_ret(R);
  
  forward(R, p);
  value(R, p, V);
  for (i = 0; i < NS - 1; i++)
    printf("%6.2f%s ", V[i], Astr[p[i]]);
  printf("\n");
}
