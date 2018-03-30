# vcg_viz_tool
Jupyter notebook based tool to visualize clinical and simulated VCGs

## Notes
* This will not run on github or nbviewer. It needs to be downloaded or cloned
and run through jupyter notebook.
* Add the node files in the .zip folder to the proper folders in your database

e.g. LBBB_dyssynch/BiV1/BiV1_nodes.txt

* Everything is in the Continuity frame of reference so the projection labels
(Transverse, Frontal, Left Sagittal) are not accurate.  If I have the rotations,
I’d like to display everything in the clinical reference frame as it may make more sense.

* The most buggy feature is the four click buttons to move the pacing location spatially.
Basically, it errors when you get to the edge of the region sometimes.
I don’t think it is worth fixing this until after the current deadlines.
* Also, BiV2 may be particularly buggy because the pacing locations were different.