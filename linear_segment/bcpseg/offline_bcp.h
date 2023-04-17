#ifndef __OFFLINE_BCP_H__
#define __OFFLINE_BCP_H__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "src/labeled_aiarray/labeled_augmented_array.h"

//-------------------------------------------------------------------------------------

double const_prior(double r, int l);

double calculate_betaT(const double x[], int start, int end, double mean, int n, int length);

double calculate_prob(const double x[], int start, int end, double muT, double nuT, double scale, int length);

double gaussian_obs_log_likelihood(const double *data, int t, int s, int length);

double calculate_tmp_cond(double *Pcp, double *P, double *Q, double *g, int n, int j, int t);

double *offline_changepoint_detection(const double *data, int n, double truncate);

double *array_exp(double *arr, int n1);

void offline_bcp_segment_labeled(const double values[], labeled_aiarray_t *segments, const char *label, int length, double truncate, double cutoff);

#endif