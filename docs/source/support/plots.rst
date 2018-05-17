绘图 Matplotlib
===============
常见的Python绘图库是 Matplotlib_ ，点击查看教程。
seaborn_ 在matplotlib的基础上进行开发，简化了调用方式，更容易上手，可以无缝集成到使用matplotlib的代码中。

.. _Matplotlib: https://matplotlib.org/users/pyplot_tutorial.html
.. _seaborn: https://seaborn.pydata.org/

常见流程
--------
步骤包括：导入包，初始化，绘制以及保存。
::
    #导入
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    #初始化
    plt.figure(figsize=(10, 10))

    #绘制
    plt.plot(x, y)

    #保存
    plt.savefig("Normalize.png")

常见类型
--------
散点图
::
    #空心圆，蓝色边线
    plt.scatter(list_x, list_y, facecolors='none', edgecolors='b')
直方图
::
    plt.hist(df.norm_by_GC, bins=300, density=True)

结构控制
-----------
绘制多幅图
::
    #绘制2X2的图，大小为10X10；
    plt.figure(figsize=(10, 10))
    #左上
    plt.subplot(2, 2, 1)
    plt.plot()
    #右上
    plt.subplot(2, 2, 2)
    plt.plot()
    ...
    plt.savefig("Multiple.png")