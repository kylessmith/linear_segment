# Linear Segmentation

[![Build Status](https://travis-ci.org/kylessmith/bcpseg.svg?branch=master)](https://travis-ci.org/kylessmith/linear_segmentation) [![PyPI version](https://badge.fury.io/py/bcpseg.svg)](https://badge.fury.io/py/linear_segmentation)
[![Coffee](https://img.shields.io/badge/-buy_me_a%C2%A0coffee-gray?logo=buy-me-a-coffee&color=ff69b4)](https://www.buymeacoffee.com/kylessmith)

linear_segmentation using Bayesian Change Point Segmentation or Circular Binary segmentation.


## Install

If you dont already have numpy and scipy installed, it is best to download
`Anaconda`, a python distribution that has them included.  
```
    https://continuum.io/downloads
```

Dependencies can be installed by:

```
    pip install -r requirements.txt
```

PyPI install, presuming you have all its requirements installed:
```
	pip install linear_segment
```

## Usage

```python
from linear_segment import segment
import numpy as np

# Create data
np.random.seed(10)
T = 50
x = np.zeros(T)
x[10:20] = 1.0
x[30:40] = 1.0

labels = np.repeat("a", T)   # "a" is a dummy label

# Calculate segments
segments = segment(x, labels, method="online_both", cutoff=0.3, offset=5)
print(segments)

LabeledIntervalArray
   (0-10, a)
   (10-20, a)
   (20-30, a)
   (30-40, a)
   (40-50, a)
segments = segment(x, labels, method="cbs", shuffles=200, p=0.05)
print(segments)

LabeledIntervalArray
   (0-10, a)
   (10-20, a)
   (20-30, a)
   (30-40, a)
   (40-50, a)

```

