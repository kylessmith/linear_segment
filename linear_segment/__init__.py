"""Init for BCP."""

from __future__ import absolute_import
from .segment import segment
from .bcpseg.bcp_segment import bcpseg
from .bcpseg.bcp_segment import get_include
from .cbseg.cbseg import determine_cbs_stat, determine_t_stat, determine_cbs, cbs_segment, validate, init_segment

# This is extracted automatically by the top-level setup.py.
__version__ = '1.2.1'
__author__ = "Kyle S. Smith"

__doc__ = """\

linear_segment
==============

.. currentmodule:: linear_segment

.. autofunction:: segment
    
"""