// src/gauss.c
#define _POSIX_C_SOURCE 200809L
#include "../include/gauss.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <string.h>


typedef struct {
    double *A;
    int n;
    int k;
    int row_start;
    int row_end;
} WorkerArg;


static inline double *ROW(double *A, int ncols, int i) { return A + (size_t)i * ncols; }


void *worker_eliminate(void *arg) {
    WorkerArg *wa = (WorkerArg *)arg;
    double *A = wa->A;
    int n = wa->n;
    int k = wa->k;
    int ncols = n + 1;
    double *pivot_row = ROW(A, ncols, k);
    double pivot = pivot_row[k];
    if (fabs(pivot) < 1e-18) return NULL;
    for (int i = wa->row_start; i <= wa->row_end; ++i) {
        double *rowi = ROW(A, ncols, i);
        double factor = rowi[k] / pivot;
        rowi[k] = 0.0;
        for (int j = k + 1; j < ncols; ++j) rowi[j] -= factor * pivot_row[j];
    }
    return NULL;
}


void swap_rows(double *A, int ncols, int r1, int r2) {
    if (r1 == r2) return;
    double *row1 = ROW(A, ncols, r1);
    double *row2 = ROW(A, ncols, r2);
    for (int j = 0; j < ncols; ++j) {
        double t = row1[j]; row1[j] = row2[j]; row2[j] = t;
    }
}


double *generate_random_system(int n, unsigned int seed) {
    int ncols = n + 1;
    double *A = malloc((size_t)n * ncols * sizeof(double));
    if (!A) { perror("malloc"); exit(1); }
    srand(seed);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            double v = ((double)rand() / RAND_MAX) * 10.0 - 5.0;
            ROW(A, ncols, i)[j] = v;
            sum += fabs(v);
        }
        ROW(A, ncols, i)[i] += sum + 1.0;
        ROW(A, ncols, i)[n] = ((double)rand() / RAND_MAX) * 10.0 - 5.0;
    }
    return A;
}


int gaussian_elim_mt(double *A, int n, int max_threads) {
    int ncols = n + 1;
    for (int k = 0; k < n; ++k) {
        int pivot_row = k; double maxv = fabs(ROW(A, ncols, k)[k]);
        for (int i = k + 1; i < n; ++i) {
            double v = fabs(ROW(A, ncols, i)[k]);
            if (v > maxv) { maxv = v; pivot_row = i; }
        }
        if (maxv < 1e-18) { fprintf(stderr, "Matrix is singular or nearly singular at column %d\n", k); return -1; }
        if (pivot_row != k) swap_rows(A, ncols, pivot_row, k);


        int rows_to_elim = n - (k + 1);
        if (rows_to_elim <= 0) continue;


        int p = max_threads <= 0 ? 1 : max_threads;
        if (p > rows_to_elim) p = rows_to_elim;


        pthread_t *threads = malloc((size_t)p * sizeof(pthread_t));
        WorkerArg *args = malloc((size_t)p * sizeof(WorkerArg));
        if (!threads || !args) { perror("malloc"); exit(1); }


        int base = rows_to_elim / p; int rem = rows_to_elim % p; int cur = k + 1;
        for (int t = 0; t < p; ++t) {
            int cnt = base + (t < rem ? 1 : 0);
            if (cnt <= 0) { args[t].A = A; args[t].n = n; args[t].k = k; args[t].row_start = 1; args[t].row_end = 0; continue; }
            args[t].A = A; args[t].n = n; args[t].k = k; args[t].row_start = cur; args[t].row_end = cur + cnt - 1; cur += cnt;
            if (pthread_create(&threads[t], NULL, worker_eliminate, &args[t]) != 0) {
                perror("pthread_create"); worker_eliminate(&args[t]); threads[t] = 0;
            }
        }
        for (int t = 0; t < p; ++t) if (threads[t]) pthread_join(threads[t], NULL);
        free(threads); free(args);
    }
    return 0;
}


void back_substitution(double *A, int n, double *x) {
    int ncols = n + 1;
    for (int i = n - 1; i >= 0; --i) {
        double *row = ROW(A, ncols, i);
        double s = row[n];
        for (int j = i + 1; j < n; ++j) s -= row[j] * x[j];
        x[i] = s / row[i];
    }
}


#include <time.h>


double timed_solve(int n, int max_threads, double *A, double *x_out) {
    struct timespec t0, t1; clock_gettime(CLOCK_MONOTONIC, &t0);
    if (gaussian_elim_mt(A, n, max_threads) != 0) { fprintf(stderr, "Elimination failed\n"); return -1.0; }
    double *x = malloc((size_t)n * sizeof(double)); back_substitution(A, n, x);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    if (x_out) {
        memcpy(x_out, x, (size_t)n * sizeof(double));
    }
    free(x);
    return elapsed;
}