#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "offline_bcp.h"

double lgamma(double x);

void display_1Darray(double *arr, int n)
{
	int i;
	for (i = 0; i < n; i++)
	{
		printf("%f  ", arr[i]);
	}
	printf("\n");
}

void display_2Darray(double *arr, int n1, int n2)
{
	int i;
	int j;
	for (i = 0; i < n2; i+=n1)
	{
		for (j = i; j < i+n1; j++)
		{
			printf("%f  ", arr[j]);
		}
		printf("\n");
	}
	printf("\n");
}


double logaddexp(double x1, double x2)
{
    return log(exp(x1) + exp(x2));
}


double calculate_mean(const double x[], int start, int end, int length)
{
    double sum = 0;
    double n = 0;
	end = fmin(end, length);

    int i;
    for (i = start; i < end; i++)
    {
        sum += x[i];
        n++;
    }

    return sum / n;
}


double logsumexp(double x[], int n)
{
    double sum = 0;

    int i;
    for (i = 0; i < n; i++)
    {
        sum += exp(x[i]);
    }

    return log(sum);
}


double *array_exp(double *arr, int n1)
{
	int i;
	int j;
    double *sum_arr = (double *)malloc(sizeof(double) * n1);
    for (i = 0; i < n1; i++)
	{
        sum_arr[i] = 0;
    }

	for (i = 0; i < n1; i++)
	{
		int start = i * n1;
        for (j = i; j < n1; j++)
		{
            sum_arr[j] += exp(arr[start+j]);
		}
	}
	
    return sum_arr;
}


inline double exp1(double x)
{
	x = 1.0 + x /256.0;
	x *= x; x *= x; x *= x; x *= x;
	x *= x; x *= x; x *= x; x *= x;

	return x;
}


double *array_concat(double *x, double y, int length)
{
	// Create new array
	double *new_arr = (double *)malloc(sizeof(double) * length+1);
	
	// Iterate over values
	int i;
	for (i = 0; i < length; i++)
	{
		new_arr[i] = x[i];
	}

	// Add new value
	new_arr[i] = y;

	return new_arr;
}