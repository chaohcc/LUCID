## 这个目录，在2025年3月30日开始，添加delete 功能。 不包含delete 功能的目录，添加了压缩包：
backup_without_delete_xxxx 中

delete 功能的实现，需要在每一个leaf piece （含m_piece\b_piece）中添加delete_flag
类似于location bitmap

添加delete 功能之后，需要考虑：
1）adaptive LRU 的更新
2）delete key 在leaf piece 和在page 中的回收问题。

Features

inner rebuild as the number of new models retrained
reorganize the leaf level
when the retrained new models > 2, use link, that insert the model from number 2 into the buffer piece, and so on.


merge the read page into memory as the holes in the page reach the threshold

use a flattened structure to reduce the tracking usage
    pageid, num_in_page, bitmap, offsets, all use the data type of pageid (the variate with the largest length - sizeof ), alough they should have different data types in the first glance.
     for example, page id (unsigned int), num in page (unsigned int),
     bitmap (bool, organized as unsigned int with every 32 bits)
     offset (unsigned short, situate the high and low addresses, respectively)


model_based insert in buffer
等比例缩放leaf piece 的buffer piece
when the num_keys in buffer == buffer->data_capacity_
retrain


