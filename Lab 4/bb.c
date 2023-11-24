#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int epsilon = 0.000001;

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

void bound(node_t *p, set_t **h, double *zp, double *x) {
  // zp is a pointer to max z found so far
  if (p->z > *zp) {
    *zp = p->z;
    memcpy(x, p->x, p->n * sizeof(double)); // save best x
  }
  // remove and delete all nodes q in h with q.z < p.z
  struct node_t *q = *h;
  struct node_t *prev = NULL;
  set_t *current = h;
  while (current->node != NULL) {
    if (current->node->z < p->z) {
      remove_node(h, current->node);
    }
    if (current->succ != NULL) {
      set_t *temp = current;
      current = current->succ;
    } else {
      return;
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
      if (q->x[h] < min || q->x[h] > max) {
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
  delete_node(q);
  return NULL;
}

double function(int m, int n, double *a, double *b, double *c, double *x) {
  struct node_t *p = initial_node(m, n, a, b, c);
  set_t *h = create_set();
  double z = -1; // best integer solution found so far

  p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);

  if (integer(p) || !is_finite(p->z)) {
    z = p->z;
    if (integer(p)) {
      memcpy(x, p->x, p->n * sizeof(double));
      delete_node(p);
      return z;
    }
  }

  add_node(h, p);

  while (h != NULL) {
    struct node_t *p = take_node(h);

    succ(p, &h, m, n, a, b, c, p->h, 1, bp.xhc, &z, x);
    succ(p, &h, m, n, a, b, c, p->h, -1, dp.xhe, &z, x);

    delete_node(p);

    if (z == -1) {
      return NAN;
    } else {
      return z;
    }
  }
}
