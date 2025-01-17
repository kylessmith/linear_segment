import os
import pandas as pd
import numpy as np
cimport numpy as np
np.import_array()
from ailist import IntervalArray
from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_init


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


cdef cbs_stat_t _determine_cbs_stat(double[::1] values):
    cdef int length = values.size
    cdef cbs_stat_t cbs_stat = calculate_cbs_stat(&values[0], length)

    return cbs_stat

def determine_cbs_stat(double[::1] values):
    cdef cbs_stat_t cbs_stat = _determine_cbs_stat(values)

    return cbs_stat


cdef double _determine_t_stat(double[::1] values, int i):
    cdef int length = values.size
    cdef double tstat = calculate_tstat(&values[0], i, length)

    return tstat

def determine_t_stat(double[::1] values, int i):
    cdef double tstat = _determine_t_stat(values, i)

    return tstat


cdef cbs_stat_t _determine_cbs(double[::1] values, int shuffles, double p):
    cdef int length = values.size
    cdef cbs_stat_t cbs_stat = calculate_cbs(&values[0], length, shuffles, p)

    return cbs_stat

def determine_cbs(double[::1] values, int shuffles=1000, double p=0.05):
    cdef cbs_stat_t cbs_stat = _determine_cbs(values, shuffles, p)

    return cbs_stat


cdef void _segment_labeled(double[::1] values, labeled_aiarray_t *c_lsegments, char *label, int shuffles, double p):
    cdef int length = values.size
    calculate_segment_labeled(&values[0], length, c_lsegments, label, shuffles, p)

def init_segment(np.ndarray values, np.ndarray labels=None, int shuffles=1000, double p = 0.05):
    """
    Implementation of a Circular Biniary Segmentation
    algorithm to segment values into Intervals

    Parameters
    ----------
        values : numpy.ndarray
			Floats to segment
        labels : numpy.ndarray
            Labels for breakup segmentation
        shuffles : int
            Number of shuffles conduct
        p : double
            Pvalue cutoff

    Returns
    -------
        segments : IntervalArray or LabeledIntervalArray
			Segment intervals
    """

    # Determine if labels are provided
    is_array = False
    if labels is None:
        is_array = True
        labels = np.repeat("_IntervalArray", len(values))

    # Call C segmentation function
    cdef bytes label_name
    cdef labeled_aiarray_t *c_lsegments
    # Find unique labels
    unique_labels = pd.unique(labels)
    c_lsegments = labeled_aiarray_init()
    for label in unique_labels:
        label_values = values[labels==label]
        label_name = label.encode()
        _segment_labeled(label_values, c_lsegments, label_name, shuffles, p)

    # Wrap C intervals
    lsegments = LabeledIntervalArray()
    lsegments.set_list(c_lsegments)
    if is_array:
        segments = IntervalArray()
        segments._wrap_laia(lsegments)

        return segments

    else:
        return lsegments


cdef void _validate_labeled(double[::1] values, LabeledIntervalArray segments, labeled_aiarray_t *vsegments, char *label, int shuffles, double p):
    cdef int length = values.size
    calculate_validate_labeled(&values[0], length, segments.laia, vsegments, label, shuffles, p)
    
def validate(np.ndarray values,
             segments,
             labels= None,
             int shuffles=1000,
             double p = 0.05):
    """
    Implementation of a Circular Biniary Segmentation validation
    algorithm to validate segment Intervals

    Parameters
    ----------
        values : numpy.ndarray
			Floats to segment
        labels : numpy.ndarray
            Labels for breakup segmentation
        segments : IntervalArray or LabeledIntervalArray
            Predicted segments to validate
        shuffles : int
            Number of shuffles conduct
        p : double
            Pvalue cutoff

    Returns
    -------
        valid_segments : IntervalArray or LabeledIntervalArray
			Valid segment intervals
    """

    # Call Validate C function
    cdef int i
    cdef bytes label_name
    cdef labeled_aiarray_t *c_vlsegments
    cdef LabeledIntervalArray label_segments
    cdef LabeledIntervalArray vlsegments
    cdef double[::1] values_mem
    if isinstance(segments, IntervalArray):
        # Iterate over labels
        label_name = b"_IntervalArray"
        c_vlsegments = labeled_aiarray_init()
        label_segments = segments._laia.get("_IntervalArray")
        values_mem = values
        _validate_labeled(values_mem, label_segments, c_vlsegments, label_name, shuffles, p)
    else:
        # Iterate over labels
        c_vlsegments = labeled_aiarray_init()
        for label in segments.unique_labels:
            label_name = label.encode()
            label_segments = segments.get(label)
            values_mem = values[labels == label]
            _validate_labeled(values_mem, label_segments, c_vlsegments, label_name, shuffles, p)

    # Wrap C intervals
    vlsegments = LabeledIntervalArray()
    vlsegments.set_list(c_vlsegments)
    if isinstance(segments, IntervalArray):
        vsegments = IntervalArray()
        vsegments._wrap_laia(vlsegments)
        
        return vsegments

    else:
        return vlsegments


def cbs_segment(np.ndarray values,
                np.ndarray labels = None,
                int shuffles = 1000,
                double p = 0.05):
    """
    Implementation of a Circular Biniary Segmentation
    algorithm to segment values into Intervals

    Parameters
    ----------
        values : numpy.ndarray
			Floats to segment
        labels : numpy.ndarray
            Labels for breakup segmentation
        shuffles : int
            Number of shuffles conduct
        p : double
            Pvalue cutoff

    Returns
    -------
        segments : IntervalArray or LabeledIntervalArray
			Segment intervals
    """

    segs = init_segment(values, labels, shuffles, p)
    segs = validate(values, segs, labels, shuffles, p)

    return segs


