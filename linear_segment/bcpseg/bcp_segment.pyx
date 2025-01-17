import os
import pandas as pd
import numpy as np
cimport numpy as np
np.import_array()
cimport cython
from libc.stdint cimport uint32_t, int32_t, int64_t
from ailist import IntervalArray
#from AIList_core cimport AIList, ailist_t, ailist_init


def get_include():
    """
    Get file directory if C headers
    
    Arguments:
    ---------
        None
    Returns:
    ---------
        str (Directory to header files)
    """

    return os.path.split(os.path.realpath(__file__))[0]


cdef void _offline_bcpseg_labeled(const double[::1] values, labeled_aiarray_t *c_segments, const char *label, double truncate, double cutoff):
    # Find length of values
    cdef int length = values.size

    # Segment
    offline_bcp_segment_labeled(&values[0], c_segments, label, length, truncate, cutoff)


cdef void _online_bcpseg_labeled(const double[::1] values, labeled_aiarray_t *c_segments, const char *label, double cutoff, double hazard):
    # Find length of values
    cdef int length = values.size

    # Segment
    online_bcp_segment_labeled(&values[0], c_segments, label, length, cutoff, hazard)


cdef void _online_bcpseg_both_labeled(const double[::1] forward_values, const double[::1] reverse_values,
                                   labeled_aiarray_t *c_segments, const char *label, double cutoff, double hazard, int offset):
    # Find length of values
    cdef int length = forward_values.size

    # Segment
    online_bcp_both_labeled(&forward_values[0], &reverse_values[0], c_segments, label, length, cutoff, hazard, offset)


def bcpseg(np.ndarray values,
           np.ndarray labels = None,
           double truncate = -100,
           double cutoff=0.75,
           str method = "online_both",
           double hazard = 100,
           int offset = 10):
    """
    Implementation of a Bayesian Change Point Detection
    algorithm to segment values into Intervals

    Parameters
    ----------
        values : numpy.ndarray
			Floats to segment
        labels : numpy.ndarray
            Labels for breakup segmentation
        truncate : float
            Tolerance during offline segmentation [default=-100]
        cutoff : float
            Probability threshold for determining segment bounds [default=0.75]
        method : str
            Method to use: offline, online, online_both [default:'online']
        hazard : float
            Expected typical segment length [default:100]
        offset : int
            Number to skip before calculating probability [default:10]

    Returns
    -------
        segments : AIList
			Segment intervals
    """

    # Initilaize segments
    #cdef IntervalArray segments
    #cdef aiarray_t *c_segments
    cdef LabeledIntervalArray lsegments
    cdef labeled_aiarray_t *c_lsegments
    cdef const double[::1] label_values
    #cdef np.ndarray unique_labels
    cdef bytes label
    cdef const char *label_name

    # Determine if labels are provided
    is_array = False
    if labels is None:
        is_array = True
        labels = np.repeat("_IntervalArray", len(values))

    # Determine method
    if method == "online":
        # Find unique labels
        labels = labels.astype(bytes)
        unique_labels = pd.unique(labels)
        c_lsegments = labeled_aiarray_init()
        for label in unique_labels:
            label_values = values[labels==label]
            _online_bcpseg_labeled(label_values, c_lsegments, label, cutoff, hazard)
    
    elif method == "online_both":
        # Find unique labels
        labels = labels.astype(bytes)
        unique_labels = [bytes(l) for l in pd.unique(labels)]
        c_lsegments = labeled_aiarray_init()
        for label in unique_labels:
            label_values = values[labels==label]
            reverse_values = np.ascontiguousarray(np.flip(label_values))
            label_name = label
            if label_values.size > 3:
                _online_bcpseg_both_labeled(label_values, reverse_values, c_lsegments, label_name, cutoff, hazard, offset)
            else:
                labeled_aiarray_add(c_lsegments, 0, label_values.size, label_name)
    else:
        # Find unique labels
        labels = labels.astype(bytes)
        unique_labels = [bytes(l) for l in pd.unique(labels)]
        c_lsegments = labeled_aiarray_init()
        for label in unique_labels:
            label_values = values[labels==label]
            _offline_bcpseg_labeled(label_values, c_lsegments, label, truncate, cutoff)

    # Wrap segments
    lsegments = LabeledIntervalArray()
    lsegments.set_list(c_lsegments)
    if is_array:
        segments = IntervalArray()
        segments._wrap_laia(lsegments)

        return segments

    else:
        return lsegments
