# NEON crown area
**Calculate crown area of woody plants as seen from above based on height and diameter measurements**

The dataset **NEON.DOM.SITE.DP1.10098.001 - Woody plant vegetation structure** provides structure measurements, including height, canopy diameter, and stem diameter, as well as mapped position of individual woody plants. However, crowns can overlap and they can also be fully contained within other crowns.

![Example tree crowns](./ABBY_001.pdf)

In this script, overlapping crown areas are assigned to the taller individual, or split among individuals with the same height. Crown areas of smaller individuals completely covered by taller individuals are omitted. 

