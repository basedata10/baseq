生成 Excel
===============
关于在Python中读取和生成Excel，教程_ 。我们使用 xlsxwriter_ 包操作Excel。

.. _教程: https://www.datacamp.com/community/tutorials/python-excel-tutorial
.. _xlsxwriter: http://xlsxwriter.readthedocs.io/

创建，写入cell以及保存
--------------------
格式以及其他
::
    import xlsxwriter
    workbook = xlsxwriter.Workbook('QC.xlsx')
    workbook.formats[0].set_font_size(12)
    workbook.formats[0].set_font_name('arial')
    format_main = workbook.add_format({'bold': False, 'font_size': 12, 'font_name': 'arial'})
    format_header = workbook.add_format({'bold': True, 'font_size': 15, 'font_name': 'arial'})
    #prepare Page...
    qcpage = workbook.add_worksheet("Report")
    qcpage.set_column('D:D', 40)
    qcpage.set_column('E:E', 40)
    qcpage.write('A1', 'Sample', format_header)
    qcpage.write('B1', 'MeanQuality', format_header)
    qcpage.write('C1', 'BiasIndex', format_header)
    qcpage.write('D1', 'BasePlot', format_header)
    qcpage.write('E1', 'QualityPlot', format_header)

    #build the Excel...
    for idx, sample in enumerate(samples):
        print(idx, sample)
        result = fastq_basecontent_quality(sample[0], sample[1])
        qcpage.set_row(idx+1, 120)
        qcpage.write(idx+1, 0, sample[0], format_main)
        qcpage.write(idx+1, 1, result[2], format_main)
        qcpage.write(idx+1, 2, result[3], format_main)
        qcpage.insert_image(idx+1, 3, result[0], {"x_scale":0.7, "y_scale":0.7, 'x_offset': 5, 'y_offset': 5})
        qcpage.insert_image(idx+1, 4, result[1], {"x_scale":0.7, "y_scale":0.7, 'x_offset': 5, 'y_offset': 5})

    workbook.close()

pandas数据框
---------------
创建多个sheet
::
    writer = pd.ExcelWriter('VCF_stats.xlsx', engine='xlsxwriter', options={'font_name':'arial'})
    pd.DataFrame(results, columns=["sample", "counts", "mean_depth", "GT_01", "GT_11"]).to_excel(writer, sheet_name='VCF')
    pd.DataFrame(MAF, columns=["sample"] + [str(round(x/100, 2)) for x in range(50)]).to_excel(writer, sheet_name='MAF')
    writer.book.formats[0].set_font_name('arial')
    writer.save()


