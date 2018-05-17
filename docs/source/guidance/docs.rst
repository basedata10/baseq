.. _cmd:

开发文档
========
A document will helps you and others better understands your code.
文档会帮助代码的开发，完善以及协作。
我们使用sphinx帮助打包生成文档。文档为 rst_ 格式。
我们可以参考 request_ 库的文档，以及 github_ 源码。

.. _rst: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#literal-blocks
.. _request: http://docs.python-requests.org/en/master/#
.. _github: https://raw.githubusercontent.com/requests/requests/master/docs/index.rst

rst语法
----------

章节
""""

列表包括有序列表和无序列表
::
     Section Title
    ===============

    The Section Title
    -----------------

    The Section Title
    ^^^^^^^^^^^^^^^^^

    The Section Title
    """"""""""""""""""

列表
"""""

* This is a bulleted list.
* It has two items, the second
  item uses two lines.

1. This is a numbered list.
2. It has two items too.

#. This is a numbered list.
#. It has two items too.

必须在上面和下面加一个空行。
::
    * This is a bulleted list.
    * It has two items, the second
      item uses two lines.

    1. This is a numbered list.
    2. It has two items too.

    #. This is a numbered list.
    #. It has two items too.

代码
""""


图像
""""

.. image:: ../_static/images.jpeg

注意必须有一个空行在上面。
::

    .. image:: ../_static/images.jpeg

表格
""""
生成表格的简单方式
::
    =====  =====  ======
       Inputs     Output
    ------------  ------
      A      B    A or B
    =====  =====  ======
    False  False  False
    True   False  True
    False  True   True
    True   True   True
    =====  =====  ======

=====  =====  ======
   Inputs     Output
------------  ------
  A      B    A or B
=====  =====  ======
False  False  False
True   False  True
False  True   True
True   True   True
=====  =====  ======

链接
""""
链接文字末尾加"_"，前后有空格；另起一行记录链接的内容；
::
    我们使用sphinx帮助打包生成文档。文档为 rst_ 格式。

    .. _rst: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#literal-blocks

.. _rst: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#literal-blocks

代码文档
----------

automodule
""""""""""""""
打包一个模块的说明文档。
::
    .. autoclass:: baseq.bam.BAMTYPE
        :members:

autofunction
""""""""""""""""
打包一个函数的说明文档。
::
    .. autofunction:: baseq.bam.BAMTYPE

函数文档的示例
::
    Barcode split into 16 files according to the valid barcode in the bcstats files.

    #. Determine whether the last base mutates;
    #. Filter by whitelist;

    :param protocol: 10X/Dropseq/inDrop.
    :param name: barcode_count.
    :param bcstats: Valid Barcode.
    :param output: (./bc_stats.txt)

    Return:
        The splitted reads will be write to XXXX/split.AA.fa

autoclass
""""""""""""
打包一个类的说明文档。
::
    .. autoclass:: baseq.bam.BAMTYPE

click command
"""""""""""""""
对于click command，使用如下方式生成文档
::
    .. click:: baseq.drops.cmd:cli
       :prog: baseq-Drop
       :show-nested:



External hyperlinks, like Python_.
.. _Python: http://www.python.org/