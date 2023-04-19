Tutorial
=========

.. code-block:: python

	from linear_segment import segment
	import numpy as np
	
	# Segment
	values = np.random.random(1000)
	values[100:200] = values[100:200] + 2
	labels = np.repeat('a', len(values))
	segments = segmentvalues, labels)
	segments
	# Interval(1-100, 'a')
	# Interval(100-200, 'a')
	# Interval(200-1000, 'a')