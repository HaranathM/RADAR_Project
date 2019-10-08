# **RADAR Project**

1. *Training Cells and Guard Cells*
* Training Cells with very small values (say `2 * 4` size) and Guard Cells are of `1 * 2` size yield false alarms.
* As more number of training cells(like `8 * 8` or more) and guard cells ( `5 * 5` or more)are giving accurate, some value in that range is considered

2. *Offset Value*
* The offset is randomly taken as 7dB after trail and error. 

* A value less than 5dB is of no use in this case. It's giving more false alarms for the corresponding Trail and Guard Cell values.
* While a value more than 8dB is making all the CUT values less than the threshold values.

3. *Non-Threshold cells at edges*

* We made sure that the the iterations are going to happen only around the Cell Under Test. 
* In other words, the margins in the row and column directions are limited by (Gr+Tr) and (Gd+Td) values respectively.
* The values of threshold cells at the edges are treated as zeros and they are added to the final CFAR signal matrix.
* The other way and the easier one is to take a zeros matrix equivalent to the size of RDM and then replace its indices(same as CUT indices of RDM) with ones, whenever CUT>threshold values. 
 