#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

int initial(struct simplex_t *s, int m, int n, double **a, double *b, double *c,
            double *x, double y, int *var) {
  int k;
  k = init(s, m, n, a, b, c, x, y, var);
  return 1;
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

double xsimplex(int m, int n, double **a, double *b, double *c, double *x,
                double y, int *var, int h) {
  struct simplex_t s;
  int i, row, col;
  if (!initial(&s, m, n, a, b, c, x, y, var)) {
    s.var = 0;
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
      s.var = 0;
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
    s.var = 0;
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