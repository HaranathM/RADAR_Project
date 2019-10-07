# **RADAR Project**

1. Training Cells and Guard Cells:
* Training Cells of `2*4` size and Guard Cells are of '1-by-2' size are considered after trail and error.
* As more number of training cells(like '6-by-8' or '8-by-8' ) are giving more detections, a lesser number is considered

2. Offset Value
* The offset is randomly taken as 7dB after trail and error. 

* A value less than 5dB is of no use in this case. It's giving more false alarms for the corresponding Trail and Guard Cell values.
* While a value more than 8dB is making all the CUT values less than the threshold values.

3. Non-Threshold cells at edges

* We made sure that the the iterations are going to happen only around the Cell Under Test. 
* In other words, the margins in the row and column directions are limited by (Gr+Tr) and (Gd+Td) values respectively.
* The values of threshold cells at the edges are treated as zeros and they are added to the final CFAR signal matrix. 