Installation
------------

To run the script, first compile dedup_umi.pyx, type:

.. code:: bash

   $ python setup.py build_ext --inplace

Which will produce a ``dedup_umi.xxx.so`` file in your local directory.

If you don't see it, it may be located in ``./build`` and move it to the
same location where ``count_umi.py`` locates.


Run:
------------

Show help message:

.. code:: bash

   $ python count_umi.py -h

Test run:
------------

.. code:: bash

   $ python count_umi.py \
		-s test_data/subset_trimmed_sgrna.fastq.gz \
		-b test_data/subset_trimmed_umi.fastq.gz \
    		-w test_data/SaturnV_One_Two.csv \
		-o count_umi_output.txt
