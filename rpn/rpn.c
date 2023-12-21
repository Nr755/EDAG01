#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define STACK_SIZE 10

int stack[STACK_SIZE];
int top = -1;

static void push(int value) { stack[++top] = value; }

static int pop() {
  if (top == -1) {
    return 0;
  }
  return stack[top--];
}

int main(void) {
  int ch;
  int line = 1;

  while ((ch = getchar()) != EOF) {
    if (isdigit(ch)) {
      int num = ch - '0';
      ch = getchar();
      while (isdigit(ch)) {
        num = num * 10 + (ch - '0');
        ch = getchar();
      }
      ungetc(ch, stdin);
      if (top < STACK_SIZE - 1) {
        push(num);
      } else {
        ch = num + '0';
        goto clear;
      }
    } else if (ch == ' ') {
      continue;
    } else if ((ch == '+' || ch == '-' || ch == '*' || ch == '/') && top > 0) {
      int num2 = pop();
      int num1 = pop();
      switch (ch) {
      case '+':
        push(num1 + num2);
        break;
      case '-':
        push(num1 - num2);
        break;
      case '*':
        push(num1 * num2);
        break;
      case '/':
        if (num2 != 0) {
          push(num1 / num2);
        } else {
          goto clear;
          continue;
        }
        break;
      }
    } else if ((ch == '+' || ch == '-' || ch == '*' || ch == '/') && top <= 0) {
      goto clear;
      continue;
    } else if (ch == '\n') {
      if (top == 0) {
        printf("line %d: %d\n", line, pop());
        line++;
        while (top != -1) {
          pop();
        }
      } else {
        printf("line %d: error at %s\n", line, "\\n");
        line++;
        while (top != -1) {
          pop();
        }
      }
    } else {
    clear:
      printf("line %d: error at %c\n", line, ch);
      (line)++;
      while (top != -1) {
        pop();
      }
      while ((ch = getchar()) != '\n' && ch != EOF) {
      }
    }
  }
  return 0;
}
