// include/gauss.h
#ifndef GAUSS_H
#define GAUSS_H


#include <stddef.h>


int gaussian_elim_mt(double *A, int n, int max_threads);
void back_substitution(double *A, int n, double *x);
void swap_rows(double *A, int ncols, int r1, int r2);


double *generate_random_system(int n, unsigned int seed);


double timed_solve(int n, int max_threads, double *A, double *x_out);


#endif // GAUSS_H