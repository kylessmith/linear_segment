#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cbs.h"


index_pair_t max_pair_index(index_pair_t i1, index_pair_t i2)
{   /* Return interval with max index */
    if (i1.index < i2.index)
    {
        return i2;
    }

    return i1;
}


index_pair_t min_pair_index(index_pair_t i1, index_pair_t i2)
{   /* Return interval with min index */
    if (i1.index < i2.index)
    {
        return i1;
    }

    return i2;
}


index_pair_t max_pair_value(index_pair_t i1, index_pair_t i2)
{   /* Return interval with max t */
    if (i1.value < i2.value)
    {
        return i2;
    }

    return i1;
}


index_pair_t min_pair_value(index_pair_t i1, index_pair_t i2)
{   /* Return interval with min t */
    if (i1.value < i2.value)
    {
        return i1;
    }

    return i2;
}


double calculate_mean(double x[], int length)
{
    double sum = 0;
    int i;
    for (i = 0; i < length; i++)
    {
        sum += x[i];
    }

    return sum / length;
}


double *copy_array(double new_array[], double array[], int length)
{    
    //int i;
    //for (i = 0; i < length; i++)
    //{
        //new_array[i] = array[i];
    //}
    memcpy(new_array, array, length * sizeof(double));

    return new_array;
}


// A utility function to swap to integers 
void swap(double *a, double *b) 
{ 
    double temp = *a; 
    *a = *b; 
    *b = temp; 
} 


void shuffle_array(double array[], int length)
{
    if (length > 1)
    {
        int i;
        for (i = length-1; i > 0; i--)
        {
            // Pick a random index from 0 to i 
            int j = rand() % (i+1);
            
            // Swap arr[i] with the element at random index 
            swap(&array[i], &array[j]);
        }
    }
}

// A utility function to print an array 
void printArray(double array[], int n) 
{ 
    int i;
    for (i = 0; i < n; i++) 
        printf("%f ", array[i]); 
    printf("\n"); 
} 


double *slice_array(double new_array[], double array[], int start, int end)
{
    int new_length = end - start;

    int i;
    for (i = 0; i < new_length; i++)
    {
        new_array[i] = array[i+start];
    }

    return new_array;
}

double square(double x)
{
    return x * x;
}