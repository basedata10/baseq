生成 PPT
===============
关于在Python中使用 **pptx** 读取和生成PPT，教程_ 。快速入手_

.. _教程: https://python-pptx.readthedocs.io/en/latest/
.. _快速入手: https://python-pptx.readthedocs.io/en/latest/user/quickstart.html

创建，写入以及保存
--------------------

.. image:: https://python-pptx.readthedocs.io/en/latest/_images/hello-world.png

代码
::
    from pptx import Presentation

    prs = Presentation()
    title_slide_layout = prs.slide_layouts[0]
    slide = prs.slides.add_slide(title_slide_layout)
    title = slide.shapes.title
    subtitle = slide.placeholders[1]

    title.text = "Hello, World!"
    subtitle.text = "python-pptx was here!"

    prs.save('test.pptx')


页面
--------------------

.. image:: https://python-pptx.readthedocs.io/en/latest/_images/bullet-slide.png

代码
::
    from pptx import Presentation

    prs = Presentation()
    bullet_slide_layout = prs.slide_layouts[1]

    slide = prs.slides.add_slide(bullet_slide_layout)
    shapes = slide.shapes

    title_shape = shapes.title
    body_shape = shapes.placeholders[1]

    title_shape.text = 'Adding a Bullet Slide'

    tf = body_shape.text_frame
    tf.text = 'Find the bullet slide layout'

    p = tf.add_paragraph()
    p.text = 'Use _TextFrame.text for first bullet'
    p.level = 1

    p = tf.add_paragraph()
    p.text = 'Use _TextFrame.add_paragraph() for subsequent bullets'
    p.level = 2

    prs.save('test.pptx')

表格
---------

.. image:: https://python-pptx.readthedocs.io/en/latest/_images/add-picture.png

代码
::
    from pptx import Presentation
    from pptx.util import Inches

    img_path = 'monty-truth.png'

    prs = Presentation()
    blank_slide_layout = prs.slide_layouts[6]
    slide = prs.slides.add_slide(blank_slide_layout)

    left = top = Inches(1)
    pic = slide.shapes.add_picture(img_path, left, top)

    left = Inches(5)
    height = Inches(5.5)
    pic = slide.shapes.add_picture(img_path, left, top, height=height)

    prs.save('test.pptx')

