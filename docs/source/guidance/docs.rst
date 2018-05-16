.. _cmd:

Write a Doc
============
A document will helps you and others better understands your code.
文档会帮助代码的开发，完善以及协作。
我们使用sphinx帮助打包生成文档。文档为 rst_ 格式。
我们可以参考 request_ 库的文档，以及 github_ 源码。

.. _rst: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#literal-blocks
.. _request: http://docs.python-requests.org/en/master/#
.. _github: https://raw.githubusercontent.com/requests/requests/master/docs/index.rst

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

图像
""""

.. image:: ../_static/images.jpeg

注意必须有一个空行在上面。
::

    .. image:: ../_static/images.jpeg

链接
""""

我们使用sphinx帮助打包生成文档。文档为 rst_ 格式。
::
    我们使用sphinx帮助打包生成文档。文档为 rst_ 格式。

    .. _rst: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#literal-blocks

.. _rst: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#literal-blocks

automodule
""""""""""

打包一个模块的说明文档。
::
    .. autoclass:: baseq.bam.BAMTYPE
        :members:


External hyperlinks, like Python_.
.. _Python: http://www.python.org/