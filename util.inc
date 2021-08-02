#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>

int inv(int n, double *a, double *b) {
  int s;
  gsl_matrix_view A;
  gsl_matrix_view B;
  gsl_permutation *p;

  A = gsl_matrix_view_array(a, n, n);
  B = gsl_matrix_view_array(b, n, n);
  p = gsl_permutation_alloc(n);
  gsl_linalg_LU_decomp(&A.matrix, p, &s);
  gsl_linalg_LU_invert(&A.matrix, p, &B.matrix);
  gsl_permutation_free(p);
}
