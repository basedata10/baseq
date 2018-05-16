.. _Drops:

.. module::baseq.drops

DropSeq
==============

Commands
^^^^^^^^^
.. click:: baseq.drops.cmd:cli
   :prog: baseq-Drop
   :show-nested:

APIs
^^^^^

Extract, Count barcode
""""""""""""""""""""""
.. automodule:: baseq.drops.barcode.count
   :members:

Correct, stats barcode
"""""""""""""""""""""""
.. automodule:: baseq.drops.barcode.stats
   :members:

Barcode Split
^^^^^^^^^^^^^^
.. automodule:: baseq.drops.barcode.split
   :members:

Reads Tagging
^^^^^^^^^^^^^^

Alternative Poly Adenelation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: baseq.drops.apa.scaner.scan
.. autofunction:: baseq.drops.apa.samples.APA_usage
.. autofunction:: baseq.drops.apa.genes.scan_genes
.. autofunction:: baseq.drops.apa.UTR.scan_utr
.. autofunction:: baseq.drops.apa.samples.APA_usage

.. automodule:: baseq.drops
   :members:
   :undoc-members: