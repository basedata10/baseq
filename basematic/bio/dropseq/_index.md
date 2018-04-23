# 统计Cell barcode

## 检测Cellbarcode

- DropSeq: 提取出Pairend1中1-16位置的碱基
- inDrop：其barcode分为两部分，中间使用W1分隔，使用该序列定位两部分barcode的位置，并进行切割；
- DropSeq: 提取出Pairend1中1-12位置的碱基

## 返回结果

为一个dict, 键为barcode的序列，值为其在数据中的出现次数； 

- read_counts：总的reads数；
- cb_counts: 总的barcode数；
- reads_counts_filterd: 去除杂项barcode后的reads数；
- cb_filtered: 去除杂项cb后的barcode数；
- CB counts：一个dict包含，过滤后的barcode的counts信息
- barcode sequences : count