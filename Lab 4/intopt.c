#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int epsilon = 0.000001;
int glob = 1;

typedef struct node_t node_t;
struct node_t {
  int m;       // Constrains
  int n;       // Decision Variables
  int k;       // Parent branches on xk
  int h;       // Branch on xh
  double xh;   // xh
  double ak;   // Parent ak
  double bk;   // Parent bk
  double *min; // Lower Bounds
  double *max; // Upper Bounds
  double **a;  // A
  double *b;   // b
  double *c;   // c
  double *x;   // x
  double z;    // z
};

typedef struct set_t set_t;
struct set_t {
  set_t *succ;
  node_t *node;
};

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

struct node_t *succ(struct node_t *p, set_t **h, int m, int n, double *a,
                    double *b, double *c, int k, double ak, double bk,
                    double *zp, double *x);

double simplex(int m, int n, double **a, double *b, double *c, double *x,
               double y);

node_t *initial_node(int m, int n, double **a, double *b, double *c) {
  node_t *p = (node_t *)malloc(sizeof(node_t));
  p->a = (double **)malloc((m + 1) * sizeof(double *));
  for (int i = 0; i <= m; i++) {
    p->a[i] = (double *)malloc((n + 1) * sizeof(double));
  }
  p->b = (double *)malloc((m + 1) * sizeof(double));
  p->c = (double *)malloc((n + 1) * sizeof(double));
  p->x = (double *)malloc((n + 1) * sizeof(double));
  p->min = (double *)malloc(n * sizeof(double));
  p->max = (double *)malloc(n * sizeof(double));
  p->m = m;
  p->n = n;
  // Copy a, b, and c parameters to p
  for (int i = 0; i <= m; i++) {
    for (int j = 0; j <= n; j++) {
      p->a[i][j] = a[i][j];
    }
  }
  for (int i = 0; i <= m; i++) {
    p->b[i] = b[i];
  }
  for (int i = 0; i <= n; i++) {
    p->c[i] = c[i];
  }
  // Initialize min and max arrays
  for (int i = 0; i < n; i++) {
    p->min[i] = -1;
    p->max[i] = 1;
  }
  return p;
}

node_t *extend(node_t *p, int m, int n, double **a, double *b, double *c, int k,
               double ak, double bk) {
  node_t *q = (node_t *)malloc(sizeof(node_t));
  int i, j;
  q->k = k;
  q->ak = ak;
  q->bk = bk;
  if (ak > 0 && p->max[k] < 1) {
    q->m = p->m;
  } else if (ak < 0 && p->min[k] > 0) {
    q->m = p->m;
  } else {
    q->m = p->m + 1;
  }
  q->n = p->n;
  q->h = -1;
  q->a = (double **)malloc((q->m + 1) * sizeof(double *));
  for (i = 0; i <= q->m; i++) {
    q->a[i] = (double *)malloc((q->n + 1) * sizeof(double));
  }
  q->b = (double *)malloc((q->m + 1) * sizeof(double));
  q->c = (double *)malloc((q->n + 1) * sizeof(double));
  q->x = (double *)malloc((q->n + 1) * sizeof(double));
  q->min = (double *)malloc(n * sizeof(double));
  q->max = (double *)malloc(n * sizeof(double));
  memcpy(q->min, p->min, n * sizeof(double));
  memcpy(q->max, p->max, n * sizeof(double));
  for (i = 0; i <= m; i++) {
    for (j = 0; j <= n; j++) {
      q->a[i][j] = a[i][j];
    }
  }
  for (i = 0; i <= m; i++) {
    q->b[i] = b[i];
  }
  for (i = 0; i <= n; i++) {
    q->c[i] = c[i];
  }
  if (ak > 0) {
    if (q->max[k] == 1 || bk < q->max[k]) {
      q->max[k] = bk;
    } else if (q->min[k] == -1 || -bk > q->min[k]) {
      q->min[k] = -bk;
    }
  }
  for (i = m, j = 0; j < n; j++) {
    if (q->min[j] > -1) {
      q->a[i][j] = -1;
      q->b[i] = -q->min[j];
      i++;
    }
    if (q->max[j] < 1) {
      q->a[i][j] = 1;
      q->b[i] = q->max[j];
      i++;
    }
  }
  return q;
}

struct set_t *create_set(node_t *node) {
  struct set_t *set = (struct set_t *)malloc(sizeof(struct set_t));
  set->succ = NULL;
  set->node = node;
  return set;
}

void free_node(node_t *p) {
  int i;
  for (i = 0; i <= p->m; i++) {
    free(p->a[i]);
  }
  free(p->a);
  free(p->b);
  free(p->c);
  free(p->x);
  free(p->min);
  free(p->max);
  free(p);
}

void free_set(struct set_t *h) {
  free(h->node);
  free(h);
}

struct node_t *take_node(struct set_t *h) {
  if (h == NULL) {
    return NULL;
  }
  struct node_t *node = h->node;
  h->node = NULL;
  return node;
}

bool add_node(set_t *set, node_t *node) {
  while (set->succ != NULL) {
    set = set->succ;
  }
  set->succ = create_set(node);
  return true;
}

int is_integer(double *xp) {
  double x = *xp;
  double r = round(x); // ISO C lround
  if (fabs(r - x) < epsilon) {
    *xp = r;
    return 1;
  } else {
    return 0;
  }
}

int integer(struct node_t *p) {
  int i;
  for (i = 0; i < p->n; i = i + 1) {
    if (!is_integer(&p->x[i])) {
      return 0;
    }
  }
  return 1;
}

void bound(struct node_t *p, struct set_t **h, double *zp, double *x) {
  if (p->z > *zp) {
    *zp = p->z;
    memcpy(x, p->x, p->n * sizeof(double));
    struct set_t *curr = *h;
    struct set_t *prev = NULL;
    while (curr != NULL) {
      if (curr->node->z < p->z) {
        if (prev == NULL) {
          *h = curr->succ;
          free_node(curr->node);
          free_set(curr);
          curr = *h;
        } else {
          prev->succ = curr->succ;
          free_node(curr->node);
          free_set(curr);
          curr = prev->succ;
        }
      } else {
        prev = curr;
        curr = curr->succ;
      }
    }
  }
}

int is_finite(double x) {
  // ISO C function
  if (isnan(x) || fabs(x) == 1) {
    return 0;
  } else {
    return 1;
  }
}

int branch(struct node_t *q, double z) {
  double min, max;
  if (q->z < z) {
    return 0;
  }
  for (int h = 0; h < q->n; h = h + 1) {
    if (!is_integer(&q->x[h])) {
      if (q->min[h] == -1) {
        min = 0;
      } else {
        min = q->min[h];
      }
      max = q->max[h];
      if (floor(q->x[h]) < min || ceil(q->x[h]) > max) {
        continue;
      }
      q->h = h;
      q->xh = q->x[h];
      // delete each of a, b, c, x of q or recycle in other way
      return 1;
    }
  }
  return 0;
}

struct node_t *succ(struct node_t *p, set_t **h, int m, int n, double *a,
                    double *b, double *c, int k, double ak, double bk,
                    double *zp, double *x) {
  struct node_t *q = extend(p, m, n, a, b, c, k, ak, bk);
  if (q == NULL) {
    return q;
  }
  q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);
  if (is_finite(q->z)) {
    if (integer(q)) {
      bound(q, h, zp, x);
    } else if (branch(q, *zp)) {
      add_node(h, q);
      return q;
    }
  }
  free_node(q);
  return NULL;
}

double intopt(int m, int n, double **a, double *b, double *c, double *x) {
  struct node_t *p = initial_node(m, n, a, b, c);
  set_t *h = create_set(p);
  add_node(h, p);
  double z = -INFINITY; // best integer solution found so far
  p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);
  if (integer(p) || !is_finite(p->z)) {
    z = p->z;
    if (integer(p)) {
      memcpy(x, p->x, (p->n + 1) * sizeof(double));
    }
    free_node(p);
    free_set(h);
    return z;
  }
  branch(p, z);
  while (h != NULL) {
    p = take_node(h);
    succ(p, h, m, n, a, b, c, p->h, 1, floor(p->xh), &z, x);
    succ(p, h, m, n, a, b, c, p->h, -1, -ceil(p->xh), &z, x);
    free_node(p);
  }
  if (z == -1) {
    return NAN; // not-a-number
  } else {
    return z;
  }
}

// ------------------------------------------------------------------------------------------------------------------------

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

int main2() {
  int i, j;
  int m;
  int n;
  scanf("%d %d", &m, &n);
  double *c;
  double **a;
  double *b;

  int local_array[10];

  for (i = 0; i < 11; i += 1)
    local_array[i] = i;

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
  printf("%lf \n", res);
  free(x);
  for (i = 0; i < m; i += 1) {
    free(a[i]);
  }
  free(a);
  free(b);
  free(c);
  return 0;
}
