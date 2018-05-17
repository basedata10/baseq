.. _cmd:

文档规范
=============
A document will helps you and others better understands your code.
文档会帮助代码的开发，完善以及协作。
我们可以参考 request_ 库的文档，以及 github_ 源码。

除了了解rst的语法格式，一个好的函数文档开需要了解其规范和要求。
可以参考Google的函数文档 范例_ 。

.. _范例 : http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
.. _rst: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#literal-blocks
.. _request: http://docs.python-requests.org/en/master/#
.. _github: https://raw.githubusercontent.com/requests/requests/master/docs/index.rst

代码文档
---------------
清楚、准确、全面的描述一个函数的功能，实现方式以及调用方法，结果等需要满足一定的规范，下面的一些规则会帮助你写出更好的文档：

-
-
-

举个例子：
::
    def module_level_function(param1, param2=None, *args, **kwargs):
        """This is an example of a module level function.

        Function parameters should be documented in the ``Args`` section. The name
        of each parameter is required. The type and description of each parameter
        is optional, but should be included if not obvious.

        If \*args or \*\*kwargs are accepted,
        they should be listed as ``*args`` and ``**kwargs``.

        The format for a parameter is::

            name (type): description
                The description may span multiple lines. Following
                lines should be indented. The "(type)" is optional.

                Multiple paragraphs are supported in parameter
                descriptions.

        Args:
            param1 (int): The first parameter.
            param2 (:obj:`str`, optional): The second parameter. Defaults to None.
                Second line of description should be indented.
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.

        Returns:
            bool: True if successful, False otherwise.

            The return type is optional and may be specified at the beginning of
            the ``Returns`` section followed by a colon.

            The ``Returns`` section may span multiple lines and paragraphs.
            Following lines should be indented to match the first line.

            The ``Returns`` section supports any reStructuredText formatting,
            including literal blocks::

                {
                    'param1': param1,
                    'param2': param2
                }

        Raises:
            AttributeError: The ``Raises`` section is a list of all exceptions
                that are relevant to the interface.
            ValueError: If `param2` is equal to `param1`.

        """
    if param1 == param2:
        raise ValueError('param1 may not be equal to param2')
    return True

插入自动文档
------------------

automodule
""""""""""""""
打包一个模块的说明文档。
::
    .. automodule:: baseq.bam.BAMTYPE
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