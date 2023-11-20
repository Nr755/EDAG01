#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int glob;
double epsilon = 0.000001;

struct simplex_t {
  int m;      /* Constraints. */
  int n;      /* Decision variables. */
  int *var;   /* 0..n 1 are nonbasic. */
  double **a; /* A. */
  double *b;  /* b. */
  double *x;  /* x. */
  double *c;  /* c. */
  double y;   /* y. */
};

int initial(struct simplex_t *s, int m, int n, double **a, double *b, double *c,
            double *x, double y, int *var);

int init(struct simplex_t *s, int m, int n, double **a, double *b, double *c,
         double *x, double y, int *var) {
  int i, k;
  s->m = m;
  s->n = n;
  s->a = a;
  s->b = b;
  s->c = c;
  s->x = x;
  s->y = y;
  s->var = var;

  if (s->var == NULL) {
    s->var = (int *)calloc((m + n + 1), sizeof(int));
    for (i = 0; i < m + n; i++) {
      s->var[i] = i;
    }
    for (k = 0, i = 1; i < m; i++) {
      if (b[i] < b[k]) {
        k = i;
      }
    }
    return k;
  } else {
    return 0;
  }
}

int select_nonbasic(struct simplex_t *s) {
  int i;
  for (i = 0; i < s->n; i++) {
    if (s->c[i] > 0) {
      return i;
    }
  }
  return -1;
}

void pivot(struct simplex_t *s, int row, int col) {
  double **a = s->a;
  double *b = s->b;
  double *c = s->c;
  int m = s->m;
  int n = s->n;
  int i, j, t;
  t = s->var[col];
  s->var[col] = s->var[n + row];
  s->var[n + row] = t;
  s->y = s->y + c[col] * b[row] / a[row][col];
  for (i = 0; i < n; i = i + 1) {
    if (i != col) {
      c[i] = c[i] - c[col] * a[row][i] / a[row][col];
      glob += 1;
    }
  }
  c[col] = -c[col] / a[row][col];
  for (i = 0; i < m; i = i + 1) {
    if (i != row) {
      b[i] = b[i] - a[i][col] * b[row] / a[row][col];
    }
  }
  for (i = 0; i < m; i = i + 1) {
    if (i != row) {
      for (j = 0; j < n; j = j + 1) {
        if (j != col) {
          a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
        }
      }
    }
  }
  for (i = 0; i < m; i = i + 1) {
    if (i != row) {
      a[i][col] = -a[i][col] / a[row][col];
    }
  }
  for (i = 0; i < n; i = i + 1) {
    if (i != col) {
      a[row][i] = a[row][i] / a[row][col];
    }
  }
  b[row] = b[row] / a[row][col];
  a[row][col] = 1 / a[row][col];
}

void prepare(struct simplex_t *s, int k) {
  int m = s->m;
  int n = s->n;
  int i;

  // Make room for xm+n at s->var[n] by moving s->var[n..n+m-1] one step to the
  // right.
  for (i = m + n; i > n; i--) {
    s->var[i] = s->var[i - 1];
  }
  s->var[n] = m + n;

  // Add xm+n to each constraint.
  n = n + 1;
  for (i = 0; i < m; i++) {
    s->a[i][n - 1] = 1;
  }

  s->x = (double *)malloc((m + n) * sizeof(double));
  s->c = (double *)malloc(n * sizeof(double));
  s->c[n - 1] = 1;
  s->n = n;

  pivot(s, k, n - 1);
}

double xsimplex(int m, int n, double **a, double *b, double *c, double *x,
                double y, int *var, int h) {
  struct simplex_t s;
  int i, row, col;
  if (!initial(&s, m, n, a, b, c, x, y, var)) {
    free(s.var);
    return NAN;
  }
  while ((col = select_nonbasic(&s)) >= 0) {
    row = -1;
    for (i = 0; i < m; i = i + 1) {
      if (a[i][col] > 0 &&
          (row < 0 || b[i] / a[i][col] < b[row] / a[row][col])) {
        row = i;
      }
    }
    if (row < 0) {
      free(s.var);
      return 1; // unbounded
    }
    pivot(&s, row, col);
  }
  if (h == 0) {
    for (i = 0; i < n; i = i + 1) {
      if (s.var[i] < n) {
        x[s.var[i]] = 0;
      }
    }
    for (i = 0; i < m; i = i + 1) {
      if (s.var[n + i] < n) {
        x[s.var[n + i]] = s.b[i];
      }
    }
    free(s.var);
  } else {
    for (i = 0; i < n; i = i + 1) {
      x[i] = 0;
    }
    for (i = n; i < n + m; i = i + 1) {
      x[i] = s.b[i - n];
    }
  }
  return s.y;
}

int initial(struct simplex_t *s, int m, int n, double **a, double *b, double *c,
            double *x, double y, int *var) {
  int i, j, k;
  double w;

  k = init(s, m, n, a, b, c, x, y, var);
  if (b[k] != 0) {
    return 1; // feasible
  }

  prepare(s, k);
  n = s->n;
  s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

  for (i = 0; i < m + n; i = i + 1) {
    if (s->var[i] == m + n - 1) {
      if (fabs(s->x[i]) > epsilon) {
        free(s->x);
        return 0; // infeasible
        free(s->c);
        return 0; // infeasible
      } else {
        break; // This i will be used on the next page.
      }
    }
  }

  if (i != n) {
    // xn+m is basic. find good nonbasic.
    for (j = k = 0; k < n; k = k + 1) {
      if (fabs(s->a[i - n][k]) > fabs(s->a[i - n][j])) {
        j = k;
      }
    }
    pivot(s, i - n, j);
    i = j;
  }

  if (i < n - 1) {
    // xn+m is nonbasic and not last. swap columns i and n-1
    k = s->var[i];
    s->var[i] = s->var[n - 1];
    s->var[n - 1] = k;
    for (k = 0; k < m; k = k + 1) {
      w = s->a[k][n - 1];
      s->a[k][n - 1] = s->a[k][i];
      s->a[k][i] = w;
    }
  } else {
    // xn+m is nonbasic and last. forget it.
    free(s);
    s->c = c;
    s->c = c;
    s->y = y;
    for (k = n - 1; k < n + m - 1; k = k + 1) {
      s->var[k] = s->var[k + 1];
    }
    n = s->n = s->n - 1;
    double *t = calloc(n, sizeof(double));
    for (k = 0; k < n; k = k + 1) {
      for (j = 0; j < n; j = j + 1) {
        if (k == s->var[j]) {
          // xk is nonbasic. add ck
          t[j] = t[j] + s->c[k];
          goto next_k;
        }
      }
      // xk is basic.
      for (j = 0; j < m; j = j + 1) {
        if (s->var[n + j] == k) {
          // xk is at row j
          break;
        }
      }
      s->y = s->y + s->c[k] * s->b[j];
      for (i = 0; i < n; i = i + 1) {
        t[i] = t[i] - s->c[k] * s->a[j][i];
      }
    next_k:;
    }
    for (i = 0; i < n; i = i + 1) {
      s->c[i] = t[i];
    }
    free(t);
    free(s->x);
  }

  return 1;
}

double simplex(int m, int n, double **a, double *b, double *c, double *x,
               double y) {
  return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

int main() {
  int i, j;
  int m;
  int n;
  scanf("%d %d", &m, &n);
  double *c;
  double **a;
  double *b;

  c = calloc(n, sizeof(double));
  a = calloc(m, sizeof(double *));
  for (i = 0; i < m; i += 1) {
    a[i] = calloc(n, sizeof(double));
  }
  b = calloc(m, sizeof(double));

  for (i = 0; i < n; i += 1) {
    scanf("%lf", &c[i]);
  }

  for (i = 0; i < m; i += 1) {
    for (j = 0; j < n; j += 1) {
      scanf("%lf", &a[i][j]);
    }
  }

  for (i = 0; i < m; i += 1) {
    scanf("%lf", &b[i]);
  }

  double *x = calloc(m, sizeof(double));

  double res = simplex(m, n, a, b, c, x, 0);
  printf("%lf", res);
  free(x);
  for (i = 0; i < m; i += 1) {
    free(a[i]);
  }
  free(a);
  free(b);
  free(c);
  return 0;
}