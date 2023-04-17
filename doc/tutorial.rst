Tutorial
=========

.. code-block:: python

	from bcpseg import bcpseg
	import numpy as np
	
	# Segment
	values = np.random.random(1000)
	values[100:200] = values[100:200] + 2
	segments = bcpseg(values)
	segments
	# Interval(1-100)
	# Interval(100-200)
	# Interval(200-1000)