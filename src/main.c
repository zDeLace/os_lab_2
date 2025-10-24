// src/main.c  (обновлённый)
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include "../include/gauss.h" // поправил путь

void print_usage(char *prog) {
    printf("Usage: %s [-n size] [-t max_threads] [-s seed] [-r] [-v]\n", prog);
}

int main(int argc, char **argv) {
    int n = 1000; int max_threads = 4; unsigned int seed = 12345; int read_stdin = 0; int verbose = 0;
    int opt;
    while ((opt = getopt(argc, argv, "n:t:s:rv")) != -1) {
        switch (opt) {
            case 'n': n = atoi(optarg); break;
            case 't': max_threads = atoi(optarg); break;
            case 's': seed = (unsigned int)atoi(optarg); break;
            case 'r': read_stdin = 1; break;
            case 'v': verbose = 1; break;
            default: print_usage(argv[0]); return 1;
        }
    }
    if (n <= 0) { fprintf(stderr,"n must be >0\n"); return 1; }
    int ncols = n + 1;
    double *A = NULL;
    if (read_stdin) {
        int in_n = 0; if (scanf("%d", &in_n) != 1) { fprintf(stderr,"Failed to read n from stdin\n"); return 1; }
        if (in_n != n) { fprintf(stderr,"stdin n (%d) differs from -n (%d). Using stdin value.\n", in_n, n); n = in_n; ncols = n+1; }
        A = malloc((size_t)n * ncols * sizeof(double)); if (!A) { perror("malloc"); return 1; }
        for (size_t i = 0; i < (size_t)n * ncols; ++i) { if (scanf("%lf", &A[i]) != 1) { fprintf(stderr,"Failed to read matrix element %zu\n", i); return 1; } }
    } else { A = generate_random_system(n, seed); }

    if (verbose && n <= 20) {
        printf("Input matrix (first %d rows):\n", n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < ncols; ++j) printf("%10.4f ", A[(size_t)i * ncols + j]);
            printf("\n");
        }
    }

    printf("PID: %d\n", getpid()); fflush(stdout);

    double *A_copy = malloc((size_t)n * ncols * sizeof(double));
    if (!A_copy) { perror("malloc"); free(A); return 1; }
    memcpy(A_copy, A, (size_t)n * ncols * sizeof(double));

    double *x = malloc((size_t)n * sizeof(double));
    if (!x) { perror("malloc"); free(A); free(A_copy); return 1; }

    double elapsed = timed_solve(n, max_threads, A_copy, x);
    if (elapsed < 0) { free(A); free(A_copy); free(x); return 1; }
    printf("n=%d threads=%d elapsed=%.6f s\n", n, max_threads, elapsed);

    // Вычислим невязку ||Ax - b||_2 по оригинальной матрице A
    double sumsq = 0.0;
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        double *row = &A[(size_t)i * ncols];
        for (int j = 0; j < n; ++j) s += row[j] * x[j];
        double r = s - row[n];
        sumsq += r * r;
    }
    double resid = sqrt(sumsq);
    printf("residual ||Ax-b||_2 = %.6e\n", resid);

    if (verbose) {
        printf("Solution x (n=%d):\n", n);
        if (n <= 40) {
            for (int i = 0; i < n; ++i) printf("x[%d] = %.10g\n", i, x[i]);
        } else {
            // печатаем первые и последние элементы для компактности
            for (int i = 0; i < 5; ++i) printf("x[%d] = %.10g\n", i, x[i]);
            printf("...\n");
            for (int i = n-5; i < n; ++i) printf("x[%d] = %.10g\n", i, x[i]);
        }
    }

    free(A); free(A_copy); free(x);
    printf("Press Enter to exit...\n");
    getchar();

    return 0;
}
