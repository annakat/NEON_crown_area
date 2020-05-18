# NEON_crown_area
Calculate crown area as seen from above from tree height and crown diameter measurements

The dataset **NEON.DOM.SITE.DP1.10098.001 - Woody plant vegetation structure** provides structure measurements, including height, canopy diameter, and stem diameter, as well as mapped position of individual woody plants. However, tree crowns can overlap and crowns can also be contained within crowns.

![Example tree crowns](./ABBY_001.pdf)

The scripts calculates crown area as seen from above. Overlapping crown areas are assigned to the taller tree, or split among trees with the same height. Crowns area of smaller trees completely covered by taller trees are omitted. 

