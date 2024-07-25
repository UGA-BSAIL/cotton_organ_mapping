This code mainly realizes the vertical and horizontal distribution of cotton plants.

1. For vertical distribution, the main code is: vertical_distribution.m. It obtains the boll count and boll distribution at plant height from the prediction results of PointNet++.

2. For horizontal distribution, its main codes are Horizontal.m and branch_type.m. Horizontal.m combines the segmentation results of TreeQSM and PointNet++ to obtain the distribution of cotton bolls on each branch relative to the main stem. Branch_type.m determines the branch type.

The rest of the related environment configurations are in their own folders.