//=============================================================================
// Bayesian Change Point Segmentation
// by Kyle S. Smith
//-----------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "offline_bcp.h"
#include "src/labeled_aiarray/labeled_augmented_array.h"
#include "utilities.c"

//-----------------------------------------------------------------------------

double const_prior(double r, int l)
{
    // Constant prior

    return 1 / (double)l;
}


double calculate_betaT(const double x[], int start, int end, double mean, int n, int length)
{
    double sum = 0;
    end = fmin(end, length);

    int i;
    for (i = start; i < end; i++)
    {
        sum += pow((x[i] - mean), 2);
    }

    return 1.0 + 0.5 * sum + (((double)n)/(1.0 + (double)n)) * (pow(mean, 2) / 2.0);
}


double calculate_prob(const double x[], int start, int end, double muT, double nuT, double scale, int length)
{
    double sum = 0;
    end = fmin(end, length);

    int i;
    for (i = start; i < end; i++)
    {
        sum += log(1 + pow((x[i] - muT), 2) / (nuT * scale));
    }

    return sum;
}


double gaussian_obs_log_likelihood(const double *data, int t, int s, int length)
{
    s += 1;
    int n = s - t;
    double mean = calculate_mean(data, t, s, length);

    double muT = (n * mean) / (1 + n);
    double nuT = 1 + n;
    double alphaT = 1 + ((double)n / 2);
    double betaT = calculate_betaT(data, t, s, mean, n, length);
    double scale = (betaT * (nuT + 1)) / (alphaT * nuT);

    // splitting the PDF of the student distribution up is /much/ faster.
    // (~ factor 20) using sum over for loop is even more worthwhile
    double prob = calculate_prob(data, t, s, muT, nuT, scale, length);
    double lgA = lgamma((nuT + 1) / 2) - log(sqrt(M_PI * nuT * scale)) - lgamma(nuT/2);

    return n * lgA - (nuT + 1) / 2 * prob;
}


double calculate_tmp_cond(double *Pcp, double *P, double *Q, double *g, int n, int j, int t)
{
    int Pcp_index = (j-1) * (n-1) + (j-1);
    int P_index = (j*n) + t;
    int Q_index = j;

    double sum = 0;
    //double *tmp_cond;
    //tmp_cond = (double *)malloc(sizeof(double) * (t-(j-1)));

    int i;
    for (i = 0; i < t-(j-1); i++)
    {
        //tmp_cond[i] = Pcp[Pcp_index+i] + P[P_index+(i*n)] + Q[t+1] + g[i+1] - Q[Q_index+i];
        sum += exp1(Pcp[Pcp_index+i] + P[P_index+(i*n)] + Q[t+1] + g[i+1] - Q[Q_index+i]);
    }
    
    //return tmp_cond;
    return log(sum);
}


double *offline_changepoint_detection(const double *data, int n, double truncate)
{
    //
    //Compute the likelihood of changepoints on data.
    //Keyword arguments:
    //data                                -- the time series data
    //n                                   -- length of data
    //truncate                            -- the cutoff probability 10^truncate to stop computation for that changepoint log likelihood

    // Initialize variables
    double *Q;
    double *g;
    double *G;
    double *P;
    Q = (double *)malloc(sizeof(double) * n);
    //zero_array(Q, n);
    g = (double *)malloc(sizeof(double) * n);
    G = (double *)malloc(sizeof(double) * n);
    P = (double *)malloc(sizeof(double) * (n * n));

    // Save everything in log representation
    int t;
    for (t = 0; t < n; t++)
    {
		g[t] = log(const_prior(t, n+1));
        if (t == 0)
        {
            G[t] = g[t];
        } else {
            G[t] = logaddexp(G[t-1], g[t]);
        }
    }
    
    P[((n-1) * n) + (n-1)] = gaussian_obs_log_likelihood(data, n-1, n, n);
    Q[n-1] = P[((n-1) * n) + (n-1)];

    for (t = n-2; t >= 0; t--)
    {
        double P_next_cp = -INFINITY;
        int s;
        for (s = t; s < n-1; s++)
        {
            P[(t*n) + s] = gaussian_obs_log_likelihood(data, t, s+1, n);

            // compute recursion
            double summand = P[(t*n) + s] + Q[s + 1] + g[s + 1 - t];
            P_next_cp = logaddexp(P_next_cp, summand);

            // truncate sum to become approx. linear in time (see
            // Fearnhead, 2006, eq. (3))
            if ((summand - P_next_cp) < truncate)
            {
                break;
            }
        }
		
        P[(t*n) + (n-1)] = gaussian_obs_log_likelihood(data, t, n, n);

        // (1 - G) is numerical stable until G becomes numerically 1
        double antiG;
        if (G[n-1-t] < -1e-15)  // exp(-1e-15) = .99999...
        {
            antiG = log(1 - exp(G[n-1-t]));
        } else {
            // (1 - G) is approx. -log(G) for G close to 1
            antiG = log(-G[n-1-t]);
        }

        Q[t] = logaddexp(P_next_cp, P[(t*n) + n-1] + antiG);
    }
	
    double *Pcp = (double *)malloc(sizeof(double) * ((n-1) * (n-1)));
	
    for (t = 0; t < n-1; t++)
    {
        Pcp[t] = P[t] + Q[t + 1] + g[t] - Q[0];

        if (isnan(Pcp[t]))
        {
            Pcp[t] = -INFINITY;
        }
    }

    
    //double *tmp_cond;
    int j;
    for (j = 1; j < n-1; j++)
    {
        //printf("   j: %d\n", j);
        for (t = j; t < n-1; t++)
        {
            //tmp_cond = calculate_tmp_cond(Pcp, P, Q, g, n, j, t);
            //Pcp[(j*(n-1)) + t] = logsumexp(tmp_cond, t-(j-1));
            Pcp[(j*(n-1)) + t] = calculate_tmp_cond(Pcp, P, Q, g, n, j, t);

            if (isnan(Pcp[(j*(n-1)) + t]))
            {
                Pcp[(j*(n-1)) + t] = -INFINITY;
            }

            //free(tmp_cond);
        }
    }

    // Free arrays
    free(Q);
    free(g);
    free(G);
    free(P);

    return Pcp;
}


void offline_bcp_segment_labeled(const double values[], labeled_aiarray_t *segments, const char *label, int length, double truncate, double cutoff)
{

    // Calculate probabilities
    printf("changepoint detecting...\n");
    double *Pcp = offline_changepoint_detection(values, length, truncate);
    printf("done\n");
    double *probs = array_exp(Pcp, length-1);
    free(Pcp);

    // Determine segments
    printf("segmenting...\n");
    int is_seg = 1;
    int start = 0;
    int i;
    for (i = 0; i < length; i++)
    {
        if (probs[i] > cutoff)
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
    printf("done\n");
    labeled_aiarray_add(segments, start, i-1, label);
    free(probs);

    return;
}