#include "util.inc"

enum { LEFT, STAY, RIGHT };
enum { NONE = 0};
enum { NS = 5, NA = 3 };
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

int main() {
  int i, ip;
  const int *p;
  double V[NS];
  const int pp[][NS] = {
			 {NONE, STAY, STAY, NONE, NONE},
			 {NONE, LEFT, LEFT, NONE, NONE},
			 {NONE, RIGHT, RIGHT, NONE, NONE},
  };
  const char *s;
  const double R[NS] = {0, 0, 0, 1, 0};
  for (ip = 0; ip < (int)(sizeof pp / sizeof *pp); ip++)  {
    p = pp[ip];
    value(R, p, V);
    for (i = 0; i < NS - 1; i++) {
      s = i > 0 && i < NS - 2 ? Astr[p[i]] : "";
      printf("%6.2f%s ", V[i], s);
    }
    printf("\n");
  }
}
