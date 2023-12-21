#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_WORDS 400
#define MAX_WORD_LENGTH 12

long is_prime(long number) {
  if (number <= 1) {
    return 0;
  }

  for (long i = 2; i * i <= number; i++) {
    if (number % i == 0) {
      return 0;
    }
  }

  return 1;
}

int main() {
  char words[MAX_WORDS][MAX_WORD_LENGTH];
  long counts[MAX_WORDS] = {0};
  long num_words = 0;

  char word[MAX_WORD_LENGTH];
  long line_number = 1;

  while (fgets(word, sizeof(word), stdin) != NULL) {
    word[strcspn(word, "\n")] = 0;

    if (is_prime(line_number)) {
      printf("trying to delete %s: ", word);
      long found = 0;
      for (long i = 0; i < num_words; i++) {
        if (strcmp(words[i], word) == 0) {
          printf("deleted\n");
          for (long j = i; j < num_words - 1; j++) {
            strcpy(words[j], words[j + 1]);
            counts[j] = counts[j + 1];
          }
          num_words--;
          found = 1;
          break;
        }
      }
      if (!found) {
        printf("not found\n");
      }
    } else {
      long found = 0;
      for (long i = 0; i < num_words; i++) {
        if (strcmp(words[i], word) == 0) {
          printf("counted %s\n", word);
          counts[i]++;
          found = 1;
          break;
        }
      }
      if (!found) {
        printf("added %s\n", word);
        strcpy(words[num_words], word);
        counts[num_words] = 1;
        num_words++;
      }
    }

    line_number++;
  }

  long max_count = 0;
  char most_frequent_word[MAX_WORD_LENGTH] = {'\0'};
  for (long i = 0; i < num_words; i++) {
    if (counts[i] > max_count ||
        (counts[i] == max_count && strcmp(words[i], most_frequent_word) < 0)) {
      max_count = counts[i];
      strcpy(most_frequent_word, words[i]);
    }
  }

  printf("result: %s %ld\n", most_frequent_word, max_count);

  return 0;
}
