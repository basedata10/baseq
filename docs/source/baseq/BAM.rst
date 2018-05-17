.. _BAM:

BAM
====

Functions
______________

- Read bam file, stats the bamfile (reads, mapping ratio...);
- Get the depth for a genomic region (and visualization);
- Get the reads overlapped with a genomic region;

Design
_______
Most of the function develop based on "samtools". The version should be >=1.3.0

- samtools depth: to get the coverage depth;
- samtools view chrN:start-end : to get the overlapped reads;

Class
_______
.. autoclass:: baseq.bam.BAMTYPE
   :members: