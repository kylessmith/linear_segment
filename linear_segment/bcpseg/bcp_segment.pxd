import numpy as np
cimport numpy as np
np.import_array()
cimport cython
from libc.stdint cimport uint32_t, int32_t, int64_t, uint16_t
from ailist import IntervalArray
from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_init, labeled_aiarray_add


cdef extern from "offline_bcp.c":
    # C is include here so that it doesn't need to be compiled externally
    pass

cdef extern from "offline_bcp.h":

    void offline_bcp_segment_labeled(const double values[], labeled_aiarray_t *segments, char *label, int length, double truncate, double cutoff) nogil


cdef extern from "online_bcp.c":
    # C is include here so that it doesn't need to be compiled externally
    pass

cdef extern from "online_bcp.h":

    void online_bcp_segment_labeled(const double values[], labeled_aiarray_t *segments, const char *label, int length, double cutoff, double hazard) nogil
    void online_bcp_both_labeled(const double *forward_data, const double *reverse_data, labeled_aiarray_t *segments, const char *label, int length, double cutoff, double hazard, int offset) nogil


cdef void _offline_bcpseg_labeled(const double[::1] values, labeled_aiarray_t *c_segments, const char *label, double truncate, double cutoff)
cdef void _online_bcpseg_labeled(const double[::1] values, labeled_aiarray_t *c_segments, const char *label, double cutoff, double hazard)
cdef void _online_bcpseg_both_labeled(const double[::1] foward_values, const double[::1] reverse_values,
                                      labeled_aiarray_t *c_segments, const char *label,
                                      double cutoff, double hazard, int offset)