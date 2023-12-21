#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "error.h"
#include "poly.h"

struct poly_t {
  int coefficients;
  int degrees;
  poly_t *next;
};

int read_coefficients(const char *s, int *i) {
  int res;
  int isnegative;

  res = 0;
  isnegative = 0;
  while (s[*i] != '\0' && s[*i] != 'x') {
    if (isdigit(s[*i])) {
      res = res * 10 + s[*i] - '0';
    } else if (s[*i] == '-') {
      isnegative = 1;
    }
    ++(*i);
  }

  if (res == 0) {
    return 1;
  }

  if (isnegative) {
    res *= -1;
  }
  return res;
}

int read_degrees(const char *s, int *i) {
  int res;

  if (s[*i] != 'x') {
    return 0;
  } else if (s[++(*i)] != '^') {
    return 1;
  }

  res = 0;
  while (s[*i] != '\0' && s[*i] != ' ') {
    if (isdigit(s[*i])) {
      res = res * 10 + s[*i] - '0';
    }
    ++(*i);
  }

  return res;
}

poly_t *add_term(poly_t *poly, int coefficients, int degrees) {
  poly_t *new_term = malloc(sizeof(poly_t));
  if (new_term == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(EXIT_FAILURE);
  }
  new_term->coefficients = coefficients;
  new_term->degrees = degrees;
  new_term->next = NULL;

  if (poly == NULL) {
    return new_term;
  } else if (degrees > poly->degrees) {
    new_term->next = poly;
    return new_term;
  } else if (degrees == poly->degrees) {
    poly->coefficients += coefficients;
    free(new_term);
    return poly;
  } else {
    poly_t *current = poly;
    while (current->next != NULL && current->next->degrees > degrees) {
      current = current->next;
    }
    if (current->next != NULL && current->next->degrees == degrees) {
      current->next->coefficients += coefficients;
      free(new_term);
    } else {
      new_term->next = current->next;
      current->next = new_term;
    }
    return poly;
  }
}

poly_t *new_poly_from_string(const char *s) {
  poly_t *res;
  int i;

  if (*s == '\0') {
    return NULL;
  }

  res = malloc(sizeof(poly_t));

  i = 0;
  res->coefficients = read_coefficients(s, &i);
  res->degrees = read_degrees(s, &i);
  res->next = new_poly_from_string(s + i);

  return res;
}

void free_poly(poly_t *poly) {
  if (poly == NULL) {
    return;
  }

  free_poly(poly->next);
  free(poly);
}

poly_t *mul(poly_t *poly1, poly_t *poly2) {
  poly_t *result = NULL;

  for (poly_t *i = poly1; i != NULL; i = i->next) {
    for (poly_t *j = poly2; j != NULL; j = j->next) {
      result = add_term(result, i->coefficients * j->coefficients,
                        i->degrees + j->degrees);
    }
  }

  return result;
}

poly_t *sort_terms(poly_t *head) {
  if (head == NULL || head->next == NULL) {
    return head; // No need to sort
  }

  poly_t *sorted = NULL;

  // Extract all items from the original linked list and insert them
  // into the new list in sorted order
  while (head != NULL) {
    poly_t *current = head;
    head = head->next;

    if (sorted == NULL || current->degrees < sorted->degrees) {
      // Insert at the beginning
      current->next = sorted;
      sorted = current;
    } else {
      // Insert in the middle or at the end
      poly_t *temp = sorted;
      while (temp != NULL) {
        if (temp->next == NULL || current->degrees < temp->next->degrees) {
          current->next = temp->next;
          temp->next = current;
          break;
        }
        temp = temp->next;
      }
    }
  }

  return sorted;
}

void print_poly(poly_t *poly) {
  if (poly == NULL) {
    printf("0");
    return;
  }

  // Print the first term
  if (abs(poly->coefficients) != 1 || poly->degrees == 0) {
    printf("%d", poly->coefficients);
  } else if (poly->coefficients == -1) {
    printf("-"); // Print '-' for -1
  }

  if (poly->degrees != 0) {
    printf("x");
  }

  if (poly->degrees > 1) {
    printf("^%d", poly->degrees);
  }

  // Print the rest of the terms
  poly_t *current = poly->next;
  while (current != NULL) {
    if (current->coefficients > 0) {
      printf(" + ");
    } else {
      printf(" - ");
    }

    if (abs(current->coefficients) != 1 || current->degrees == 0) {
      printf("%d", abs(current->coefficients));
    }

    if (current->degrees != 0) {
      printf("x");
    }

    if (current->degrees > 1) {
      printf("^%d", current->degrees);
    }

    current = current->next;
  }

  printf("\n");
}
