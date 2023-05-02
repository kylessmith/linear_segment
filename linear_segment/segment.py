import numpy as np
from ailist import LabeledIntervalArray


# Local imports
from .cbseg import cbseg
from .bcpseg import bcp_segment


def segment(values: np.ndarray,
           labels: np.ndarray = None,
           truncate: float = -100,
           cutoff: float =0.75,
           method: str = "online_both",
           hazard: float = 100,
           offset: int = 10,
           shuffles: int = 1000,
           p: float = 0.05) -> LabeledIntervalArray:
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
        shuffles : int
            Number of shuffles conduct
        p : float
            Pvalue cutoff

    Returns
    -------
        segments : LabeledIntervalArray
			Segment intervals
    """

    if method == "cbs":
        segments = cbseg.cbs_segment(values,labels, shuffles, p)
        segments = cbseg.validate(values, segments,labels, shuffles, p)
    else:
        segments = bcp_segment.bcpseg(values, labels, truncate, cutoff, method, hazard, offset)

    return segments