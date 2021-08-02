/*
  c99  traj.c `pkg-config --libs --cflags gsl`
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

void traj_gen(void) {
  int t, s, s0, j, a;
  const gsl_rng_type * T;
  gsl_rng * r;

 T = gsl_rng_default;
 r = gsl_rng_alloc (T);

 for (t = 0; t < NT; t++) {
   s = rnd(r, 0, NS - 1);
   for (j = 0; ; j++) {
     traj[t][j] = s;
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
  }
}

void traj_prob(const int *p) {
  double prob;
  int t, j, a, x, y;
  for (t = 0; t < NT; t++) {
    prob = 1;
    for (j = 1; j < MAX_LEN && traj[t][j - 1] != NS - 1; j++) {
      x = traj[t][j - 1];
      y = traj[t][j];
      a = p[x];
      prob *= P[a][x][y];
    }
    Prob[t] = prob;
  }
}

double traj_score(void)
{
  double ans;
  int t;
  ans = 0;
  for (t = 0; t < NT; t++)
    ans += Ret[t] * Prob[t];
  return ans/NT;
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

void get_reward(const double *phi, double *R) {
  int i;
  for (i = 0; i < NS; i++)
    R[i] = 0;
  R[NS - 3] = phi[0];
  R[NS - 2] = phi[1];
}

int main() {
  int i;
  double S, Sm, phi[NP], V[NS];
  int p[NS] = {LEFT, LEFT, LEFT, LEFT, LEFT};
  const double R[NS] = {0, 0, 0, 1, 0};
  double Rth[NS];
  traj_gen();
  Sm = 0;
  for (phi[0] = -2.0; phi[0] < 2.0; phi[0] += 0.1) {
    for (phi[1] = -2.0; phi[1] < 2.0; phi[1] += 0.1) {
      get_reward(phi, Rth);
      forward(Rth, p);
      value(Rth, p, V);
      traj_prob(p);
      traj_ret(R);
      if ((S = traj_score()) > Sm) {
	Sm = S;
	fprintf(stderr, "score: %g %g       %g\n", phi[0], phi[1], Sm);
	for (i = 0; i < NS - 1; i++)
	  printf("%6.2f%s ", V[i], Astr[p[i]]);
	printf("\n");
      }
    }
  }
}
