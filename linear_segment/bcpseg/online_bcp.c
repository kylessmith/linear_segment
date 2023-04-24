//=============================================================================
// Bayesian Change Point Segmentation
// by Kyle S. Smith
//-----------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "online_bcp.h"
#include "src/labeled_aiarray/labeled_augmented_array.h"

//-----------------------------------------------------------------------------


studentT_t *studentT_init(double alpha, double beta, double kappa, double mu)
{   /* Initialize studentT_t object */

    // Initialize variables
    studentT_t *t = (studentT_t *)malloc(sizeof(studentT_t));
    t->alpha0 = alpha;
    t->alpha = (double *)malloc(sizeof(double));
    t->alpha[0] = alpha;
    
    t->beta0 = beta;
    t->beta = (double *)malloc(sizeof(double));
    t->beta[0] = beta;

    t->kappa0 = kappa;
    t->kappa = (double *)malloc(sizeof(double));
    t->kappa[0] = kappa;

    t->mu0 = mu;
    t->mu = (double *)malloc(sizeof(double));
    t->mu[0] = mu;

    t->n = 1;

    // Check if memory was allocated
    if (t == NULL)
    {
        fprintf (stderr, "Out of memory!!! (init)\n");
        exit(1);
    }

    return t;
}


double pdf(double x, double df)
{
    double numerator = exp(lgamma((df+1)/2)-lgamma(df/2));
    double denominator = sqrt(df*M_PI) * pow((1+pow(x,2)/df), ((df+1)/2));

    return numerator / denominator;
}


double *studentT_pdf(studentT_t *t, double x)
{
    // Initialize values
    double *pdf_values = (double *)malloc(sizeof(double) * t->n);
    // Check if memory was allocated
    if (pdf_values == NULL)
    {
        fprintf (stderr, "Out of memory!!! (studentT_pdf)\n");
        exit(1);
    }

    // Iterate over values
    int i;
    for (i = 0; i < t->n; i++)
    {
        double df = 2 * t->alpha[i];
        double loc = t->mu[i];
        double scale = sqrt(t->beta[i] * (t->kappa[i]+1) / (t->alpha[i] * t->kappa[i]));

        double y = (x - loc) / scale;
        pdf_values[i] = pdf(y, df) / scale;
    }

    return pdf_values;
}


void studentT_update_theta(studentT_t *t, double x)
{
    // Create new array
	double *new_mu = (double *)malloc(sizeof(double) * t->n+1);
    double *new_kappa = (double *)malloc(sizeof(double) * t->n+1);
    double *new_alpha = (double *)malloc(sizeof(double) * t->n+1);
    double *new_beta = (double *)malloc(sizeof(double) * t->n+1);

    // Check if memory was allocated
    if (new_mu == NULL || new_kappa == NULL || new_alpha == NULL || new_beta == NULL)
    {
        fprintf (stderr, "Out of memory!!! (studentT_update_theta)\n");
        exit(1);
    }

    // Set first values
    new_mu[0] = t->mu0;
    new_kappa[0] = t->kappa0;
    new_alpha[0] = t->alpha0;
    new_beta[0] = t->beta0;
	
	// Iterate over values
	int i;
	for (i = 0; i < t->n; i++)
	{
		new_mu[i+1] = (t->kappa[i] * t->mu[i] + x) / (t->kappa[i] + 1);
        new_kappa[i+1] = t->kappa[i] + 1.0;
        new_alpha[i+1] = t->alpha[i] + 0.5;
        new_beta[i+1] = t->beta[i] + (t->kappa[i] * pow((x - t->mu[i]), 2)) / (2.0 * (t->kappa[i] + 1.0));
	}

    // Update mu
    free(t->mu);
    t->mu = new_mu;
    free(t->kappa);
    t->kappa = new_kappa;
    free(t->alpha);
    t->alpha = new_alpha;
    free(t->beta);
    t->beta = new_beta;

    t->n++;
}


float *online_changepoint_detection(const double *data, int length, studentT_t *t, double hazard)
{
    // Initialize R
    float *R = (float *)malloc(sizeof(float) * ((length+1) * (length+1)));
    //double *R = (double *)malloc(sizeof(double) * ((length+1) * (length+1)));
    // Check if memory was allocated
    if (R == NULL)
    {
        fprintf (stderr, "Out of memory!!! (online_changepoint_detection)\n");
        exit(1);
    }
    R[0] = 1;
    
    // Iterate over data
    double x;
    int i;
    for (i = 0; i < length; i++)
    {
        x = data[i];
        // Evaluate the predictive distribution for the new datum under each of
        // the parameters.  This is the standard thing from Bayesian inference.
        double *predprobs = studentT_pdf(t, x);
       
        // Evaluate the growth probabilities - shift the probabilities down and to
        // the right, scaled by the hazard function and the predictive
        // probabilities.
        int j;
        for (j = 1; j < i+2; j++)
        {
            int R_index1 = (j * (length+1)) + (i+1);
            int R_index2 = ((j-1) * (length+1)) + i;
            R[R_index1] = R[R_index2] * predprobs[j-1] * (1-(1.0 / hazard));
        }
        
        // Evaluate the probability that there *was* a changepoint and we're
        // accumulating the mass back down at r = 0.
        double sum = 0.0;
        for (j = 0; j < i+1; j++)
        {
            int R_index = (j * (length+1)) + (i);
            sum += R[R_index] * predprobs[j] * (1.0 / hazard);
        }
        R[i+1] = sum;
        
        // Renormalize the run length probabilities for improved numerical
        // stability.
        sum = 0.0;
        for (j = 0; j < i+2; j++)
        {
            int R_index = (j * (length+1)) + (i+1);
            sum += R[R_index];
        }
        
        for (j = 0; j < i+2; j++)
        {
            int R_index = (j * (length+1)) + (i+1);
            R[R_index] = R[R_index] / sum;
        }

        // Update the parameter sets for each possible run length.
        studentT_update_theta(t, x);
        free(predprobs);
    }

    return R;    
}


void online_bcp_segment_labeled(const double *data, labeled_aiarray_t *segments, const char *label, int length, double cutoff, double hazard)
{
    // Calculate R
    studentT_t *t = studentT_init(0.1, .01, 1, 0);
    float *R = online_changepoint_detection(data, length, t, hazard);
    free(t);

    // Determine probabilities after 10
    int is_seg = 1;
    int start = 0;
    int i;
    for (i = 10; i < length; i++)
    {
        int R_index = (10 * (length+1)) + i;
        if (R[R_index] > cutoff)
        {
            if (is_seg == 0)
            {
                labeled_aiarray_add(segments, start, i-10, label);
                start = i-10;
                is_seg = 1;
            }
        } else {
            is_seg = 0;
        }
    }
    labeled_aiarray_add(segments, start, i, label);
    free(R);

    return;
}


double *online_bcp_probability(const double *data, int length, double cutoff, double hazard, int offset)
{
    // Calculate R
    studentT_t *t = studentT_init(0.1, .01, 1, 0);
    float *R = online_changepoint_detection(data, length, t, hazard);
    free(t);

    // Initialize probabilities
    double *probs = (double *)malloc(sizeof(double) * length);
    // Check if memory was allocated
    if (probs == NULL)
    {
        fprintf (stderr, "Out of memory!!! (online_bcp_probability)\n");
        exit(1);
    }

    // Assign probabilities
    int i;
    for (i = offset; i < length; i++)
    {
        int R_index = (offset * (length+1)) + i;
        probs[i-offset] = R[R_index];
    }
    // Set end to 0
    for (i = length-offset; i < length; i++)
    {
        probs[i] = 0.0;
    }
    free(R);

    return probs;
}


void segmentation_labeled(double *probability, labeled_aiarray_t *segments, const char *label, int length, double cutoff)
{
    // Initialize segments
    //labeled_aiarray_t *segments = labeled_aiarray_init();

    // Determine probabilities after cutoff
    int is_seg = 1;
    int start = 0;
    int i;
    for (i = 0; i < length; i++)
    {
        if (probability[i] > cutoff)
        {
            if (is_seg == 0)
            {
                labeled_aiarray_add(segments, start, i, label);
                start = i;
                is_seg = 1;
            }
        } else {
            is_seg = 0;
        }
    }
    labeled_aiarray_add(segments, start, i, label);

    //return segments;
    return;
}


void online_bcp_both_labeled(const double *forward_data, const double *reverse_data, labeled_aiarray_t *segments, const char *label, int length, double cutoff, double hazard, int offset)
{    

    // Determine probabilities
    double *probs_forward = online_bcp_probability(forward_data, length, cutoff, hazard, offset);
    double *probs_reverse = online_bcp_probability(reverse_data, length, cutoff, hazard, offset);

    // Iterate over probs and combine
    int i;
    // Set end
    for (i = 1; i < offset; i++)
    {
        probs_forward[length-i] = probs_reverse[i];
    }
    
    // Set middle
    for (i = offset; i < length-offset; i++)
    {
        probs_forward[i] = (probs_forward[i] + probs_reverse[length-i]) / 2.0;
    }

    // Free probs_reverse
    free(probs_reverse);

    // Segment
    segmentation_labeled(probs_forward, segments, label, length, cutoff);

    // Free probs_forward
    free(probs_forward);

    //return c_segments;
    return;
}