.. _cmd:

文档格式/rst
============
A document will helps you and others better understands your code.
文档会帮助代码的开发，完善以及协作。
我们使用sphinx帮助打包生成文档。文档为 rst_ 格式，参考 文档_ 。
我们可以参考 request_ 库的文档，以及 github_ 源码。

除了了解rst的语法格式，一个好的函数文档开需要了解其规范和要求。
可以参考Google的函数文档 范例_ 。

.. _范例 : http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
.. _rst: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#literal-blocks
.. _request: http://docs.python-requests.org/en/master/#
.. _github: https://raw.githubusercontent.com/requests/requests/master/docs/index.rst
.. _文档：http://www.sphinx-doc.org/en/stable/rest.html#tables

rst语法
----------

文字格式
"""""""""
**text** 加重
""
    **text**

章节
""""
章节的标定确定了页面的框架。在rst中，章节下面应该加上不短于标题的符号串。章节是分层次的，不同层次的章节由下面的符号串进行区分，一般推荐如下设计：
::
    Header
    =======

    Section L2
    -----------

    Section L3
    ^^^^^^^^^^^^

    Section L4/Paragraph
    """""""""""""""""""""

列表
"""""

* This is a bulleted list.
* It has two items, the second
  item uses two lines.

1. This is a numbered list.
2. It has two items too.

#. This is a numbered list.
#. It has two items too.

注意，需要在上面和下面加一个空行。
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
参考 sphinx_ 的说明。
::
    .. code-block:: python
    :emphasize-lines: 3,5

        def some_function():
           interesting = False
           print 'This line is highlighted.'
           print 'This one is not...'
           print '...but this one is.'


.. code-block:: python
   :emphasize-lines: 3,5

   def some_function():
       interesting = False
       print 'This line is highlighted.'
       print 'This one is not...'
       print '...but this one is.'

.. _sphinx: http://www.sphinx-doc.org/en/stable/markup/code.html

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

推荐使用csv table，可以指定标题，列名以及宽度，注意指令上需要空行；
::
    .. csv-table:: A CSV Table
       :header: "Treat", "Quantity", "Description"
       :widths: 15, 10, 30

       "Albatross", 2.99, "On a stick!"
       "Crunchy Frog", 1.49, "If we took the bones out, it wouldn't be
       crunchy, now would it?"
       "Gannet Ripple", 1.99, "On a stick!"

.. csv-table:: A CSV Table
   :header: "Treat", "Quantity", "Description"
   :widths: 15, 10, 30

   "Albatross", 2.99, "On a stick!"
   "Crunchy Frog", 1.49, "If we took the bones out, it wouldn't be
   crunchy, now would it?"
   "Gannet Ripple", 1.99, "On a stick!"

链接
""""
链接文字末尾加"_"，前后有空格；另起一行记录链接的内容；
::
    我们使用sphinx帮助打包生成文档。文档为 rst_ 格式。

    .. _rst: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#literal-blocks

.. _rst: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#literal-blocks