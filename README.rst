<<<<<<< HEAD
Installation
------------
=======
# INSTALL:
>>>>>>> 1c1a0451258715b42be2e0f3ed7f391d344b4256

To run the script, first compile dedup_umi.pyx, type:

.. code:: bash

   $ python setup.py build_ext --inplace

Which will produce a 'dedup_umi.xxx.so' file in your local directory.

If you don't see it, it may be located in ./build and move it to the
same location where 'count_umi.py' locates.


<<<<<<< HEAD
Run:
------------
=======
# RUN:
>>>>>>> 1c1a0451258715b42be2e0f3ed7f391d344b4256

Show help message:

.. code:: bash
<<<<<<< HEAD

=======
>>>>>>> 1c1a0451258715b42be2e0f3ed7f391d344b4256
   $ python count_umi.py -h

Test run:
------------

.. code:: bash
<<<<<<< HEAD

=======
   
>>>>>>> 1c1a0451258715b42be2e0f3ed7f391d344b4256
   $ python count_umi.py \
		-s test_data/subset_trimmed_sgrna.fastq.gz \
		-b test_data/subset_trimmed_umi.fastq.gz \
    		-w test_data/SaturnV_One_Two.csv \
		-o count_umi_output.txt
