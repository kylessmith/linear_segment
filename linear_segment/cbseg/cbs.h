#ifndef __CBS_H__
#define __CBS_H__

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include "src/labeled_aiarray/labeled_augmented_array.h"
//#include "utilities.c"

//-------------------------------------------------------------------------------------

typedef enum { false, true } bool;

typedef struct {
    int32_t start;
    int32_t end;
    double_t t;
} segment_t;

typedef struct {
    int32_t index;
    double_t value;
} index_pair_t;

typedef struct {
    int32_t max_start;
    int32_t max_end;
    double_t max_t;
    bool threshold;
} cbs_stat_t;


//-------------------------------------------------------------------------------------
//                                    cbs_structs.c
//=====================================================================================

segment_t segment_init(int32_t start, int32_t end, double_t t);

index_pair_t index_pair_init(int32_t index, double_t value);

cbs_stat_t cbs_stat_init(int32_t max_start, int32_t max_end, double_t max_t, bool threshold);

cbs_stat_t calculate_cbs_stat(double x[], int length);

double calculate_tstat(double x[], int i, int length);

cbs_stat_t calculate_cbs(double x[], int length, int shuffles, double p);


//-------------------------------------------------------------------------------------
//                                    cbs_labeled.c
//=====================================================================================

void rsegment_labeled(double x[], int start, int end, labeled_aiarray_t *segments, const char *label, int shuffles, double p);

void calculate_segment_labeled(double x[], int length, labeled_aiarray_t *segments, const char *label, int shuffles, double p);

void calculate_validate_labeled(double x[], int length,  labeled_aiarray_t *segments, labeled_aiarray_t *vsegments, const char *label, int shuffles, double p);


//-------------------------------------------------------------------------------------
//                                    utilities.c
//=====================================================================================

index_pair_t max_pair_index(index_pair_t i1, index_pair_t i2);

index_pair_t min_pair_index(index_pair_t i1, index_pair_t i2);

index_pair_t max_pair_value(index_pair_t i1, index_pair_t i2);

index_pair_t min_pair_value(index_pair_t i1, index_pair_t i2);

double calculate_mean(double x[], int length);

double *copy_array(double new_array[], double array[], int length);

void swap(double *a, double *b);

void shuffle_array(double array[], int length);

void printArray(double array[], int n);

double *slice_array(double new_array[], double array[], int start, int end);

double square(double x);


#endif