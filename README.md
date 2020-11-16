# INSTALL:
========
To run the script, first compile dedup_umi.pyx, type:

   $python setup.py build_ext --inplace

Which will produce a 'dedup_umi.xxx.so' file in your local directory. 

If you don't see it, it may be located in ./build and move it to the 
same location where 'count_umi.py' locates. 


# RUN:
=======
Show help message:

   $python count_umi.py -h

Test run:

   $python count_umi.py \
		-s test_data/subset_trimmed_sgrna.fastq.gz \
		-b test_data/subset_trimmed_umi.fastq.gz \
    		-w test_data/SaturnV_One_Two.csv \
		-o count_umi_output.txt


