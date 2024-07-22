//=============================================================================
// Circular Binary Segmentation
// by Kyle S. Smith
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "cbs.h"
//#include "utilities.c"


//-----------------------------------------------------------------------------

segment_t segment_init(int32_t start, int32_t end, double_t t)
{   /* Initialize segment */

    // Initialize interval
    segment_t i;

    // Set attributes
    i.start = start;
    i.end = end;
    i.t = t;

    return i;
}


index_pair_t index_pair_init(int32_t index, double_t value)
{   /* Initialize index pair */

    // Initialize index_pair
    index_pair_t ip;

    // Set attributes
    ip.index = index;
    ip.value = value;

    return ip;
}


cbs_stat_t cbs_stat_init(int32_t max_start, int32_t max_end, double_t max_t, bool threshold)
{   /* Initialize cbs statistic */

    // Initialize cbs_stat
    cbs_stat_t c;

    // Set attributes
    c.max_start = max_start;
    c.max_end = max_end;
    c.max_t = max_t;
    c.threshold = threshold;

    return c;
}

//-----------------------------------------------------------------------------

cbs_stat_t calculate_cbs_stat(double x[], int length)
{   /* Calculate CBS statistic from array */

    // Initialize statistic
    cbs_stat_t cbs_stat = cbs_stat_init(0,0,0, false);
    
    double mean = calculate_mean(x, length);
    double sum = 0;
    index_pair_t i_max = index_pair_init(-1, DBL_MIN);
    index_pair_t i_min = index_pair_init(-1, DBL_MAX);
    index_pair_t e0;
    double x0;

    int i;
    for (i = 0; i < length; i++)
    {
        x0 = (x[i] - mean);
        sum += x0;
        //y[i] = sum;

        e0 = index_pair_init(i, sum);
        i_max = max_pair_value(i_max, e0);
        i_min = min_pair_value(i_min, e0);
    }

    index_pair_t i0 = min_pair_index(i_max, i_min);
    index_pair_t i1 = max_pair_index(i_max, i_min);

    cbs_stat.max_start = i0.index + 1;
    cbs_stat.max_end = i1.index + 1;
    cbs_stat.max_t = pow((i1.value - i0.value), 2.0) * length;
    cbs_stat.max_t = cbs_stat.max_t / (i1.index - i0.index);
    cbs_stat.max_t = cbs_stat.max_t / (length + 1 - i1.index + i0.index);

    return cbs_stat;
}


double calculate_tstat(double x[], int i, int length)
{   /* Calculate T statistic from array */

    double s0 = 0;
    double s1 = 0;

    int j;
    for (j = 0; j < length; j++)
    {
        if (j < i)
        {
            s0 += x[j];
        }
        else
        {
            s1 += x[j];
        }
    }
    // Find means
    s0 = s0 / i;
    s1 = s1 / (length - i);

    return (length - i) * i / length * square((s0 - s1));
}

//-----------------------------------------------------------------------------

cbs_stat_t calculate_cbs(double x[], int length, int shuffles, double p)
{   /* Calculate CBS */

    cbs_stat_t cbs_stat = calculate_cbs_stat(x, length);

    if ((cbs_stat.max_end - cbs_stat.max_start) == length)
    {
        cbs_stat.threshold = false;
        return cbs_stat;
    }

    if (cbs_stat.max_start < 5)
    {
        cbs_stat.max_start = 0;
    }

    if ((length - cbs_stat.max_end) < 5)
    {
        cbs_stat.max_end = length;
    }

    int thresh_count = 0;
    double alpha = shuffles * p;
    //double tmp_xt[length];
    double *tmp_xt = malloc(length * sizeof(double));
    double *xt = copy_array(tmp_xt, x, length);
    cbs_stat_t cbs_stat_shuffled;

    int i;
    for (i = 0; i < shuffles; i++)
    {
        shuffle_array(xt, length);
        cbs_stat_shuffled = calculate_cbs_stat(xt, length);

        if (cbs_stat_shuffled.max_t >= cbs_stat.max_t)
        {
            thresh_count += 1;
        }

        if (thresh_count > alpha)
        {   
            free(tmp_xt);
            cbs_stat.threshold = false;
            return cbs_stat;
        }
    }

    cbs_stat.threshold = true;
    free(tmp_xt);

    return cbs_stat;
}
