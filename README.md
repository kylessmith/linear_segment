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
	pip install linear_segmentation
```

## Usage

```python
from linear_segmentation import segment
import numpy as np

# Create data
np.random.seed(10)
x = np.random.random(300000)
x[10000:20000] = x[10000:20000] + 0.1
x[25000:27000] = x[25000:27000] - 1
labels = np.repeat('a', len(x))

# Calculate segments
segments = segment(x, labels, method="online_both")
segments
segments = segment(x, labels, method="cbs")
segments

```

