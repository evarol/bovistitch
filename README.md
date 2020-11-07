# bovistitch


Bovistitch is a tool to stitch and vignette correct multi slide images using linearly approximated optimal transport. There are no parameters to tune other than the granularity of histogram estimation.

Requires the ImarisReader github repo from https://github.com/PeterBeemiller/ImarisReader

Use:

Run bovistitch.m, select ImarisReader path, and then select the ims files you want stitched. The output stitched/vignette corrected volumes will be in the same path as the source images.


Examples:

8x8 planar FOVs

![Demo](https://github.com/evarol/bovistitch/blob/master/fig_1.png)

![Demo](https://github.com/evarol/bovistitch/blob/master/fig_2.png)

![Demo](https://github.com/evarol/bovistitch/blob/master/fig_3.png)


2x2 FOV x 483 z-slices (z-maximum projection view):

![Demo](https://github.com/evarol/bovistitch/blob/master/fig_4.png)

![Demo](https://github.com/evarol/bovistitch/blob/master/fig_5.png)

Contributors: B. Rao, E. Varol
