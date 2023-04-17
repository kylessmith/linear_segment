//=============================================================================
// Circular Binary Segmentation
// by Kyle S. Smith
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "cbs.h"
#include "src/labeled_aiarray/labeled_augmented_array.h"
//#include "utilities.c"


//-----------------------------------------------------------------------------

void rsegment_labeled(double x[], int start, int end, labeled_aiarray_t *segments, const char *label, int shuffles, double p)
{   
    int slice_length = end - start;
    //double tmp_slice_x[slice_length];
    double *tmp_slice_x = malloc(slice_length * sizeof(double));
    double *slice_x = slice_array(tmp_slice_x, x, start, end);
    cbs_stat_t cbs_stat = calculate_cbs(slice_x, slice_length, shuffles, p);

    free(tmp_slice_x);

    int cbs_length = cbs_stat.max_end - cbs_stat.max_start;

    if (cbs_stat.threshold == false || (cbs_length < 5) || (cbs_length == slice_length))
    {
        labeled_aiarray_add(segments, start, end, label);
    }
    else
    {
        if (cbs_stat.max_start > 0)
        {
            rsegment_labeled(x, start, start + cbs_stat.max_start, segments, label, shuffles, p);
        }
        
        if (cbs_length > 0)
        {
            rsegment_labeled(x, start + cbs_stat.max_start, start + cbs_stat.max_end, segments, label, shuffles, p);
        }
        
        if (start + cbs_stat.max_end < end)
        {
            rsegment_labeled(x, start + cbs_stat.max_end, end, segments, label, shuffles, p);
        }
    }

    return;
}


void calculate_segment_labeled(double x[], int length, labeled_aiarray_t *segments, const char *label, int shuffles, double p)
{
    int start = 0;
    int end = length;

    srand(time(NULL));

    rsegment_labeled(x, start, end, segments, label, shuffles, p);

    return;
}


void calculate_validate_labeled(double x[], int length,  labeled_aiarray_t *segments, labeled_aiarray_t *vsegments, const char *label_name, int shuffles, double p)
{
    // Get label
    uint32_t l = get_label(segments, label_name);
    label_t *label = &segments->labels[l];
    
    int i;
    int S[label->ail->nr + 1];
    for (i = 0; i < label->ail->nr; i++)
    {
        S[i] = label->ail->interval_list[i].start;
    }
    S[i] = length;

    int segment_start = 0;
    int left = 0;
    double t;

    int test;
    for (test = 0; test < label->ail->nr-1; test++)
    { 
        int slice_length = S[test+2] - S[left];
        double *tmp_slice_x = malloc(slice_length * sizeof(double));
        double *slice_x = slice_array(tmp_slice_x, x, S[left], S[test+2]);

        t = calculate_tstat(slice_x, S[test+1]-S[left], slice_length);

        double threshold = 0;
        int thresh_count = 0;
        int site = S[test+1] - S[left];

        //is this needed?
        double *tmp_xt = malloc(slice_length * sizeof(double));
        double *xt = copy_array(tmp_xt, slice_x, slice_length);
        bool flag = true;

        free(tmp_slice_x);

        int k;
        for (k = 0; k < shuffles; k++)
        {
            shuffle_array(xt, slice_length);
            threshold = calculate_tstat(xt, site, slice_length);

            if (threshold > t)
            {
                thresh_count += 1;
            }

            if (thresh_count >= (p * shuffles))
            {
                flag = false;
                break;
            }
        }

        if (flag == true)
        {
            labeled_aiarray_add(vsegments, segment_start, S[test+1], label_name);
            segment_start = S[test+1];
            left += 1;
        }

        free(tmp_xt);
    }
    labeled_aiarray_add(vsegments, segment_start, S[i], label_name);

    return;
}