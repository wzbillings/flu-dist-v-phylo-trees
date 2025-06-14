# Comparing distance-based trees and phylogenetic trees for influenza

Developing universal influenza vaccines will require improved understanding of how influenza variants differ from each other. We find that temporal distances perform poorly overall, but even sequence distances which match phylogenetic distances well do not match cartographic distances based on actual immune response data.

This was my term project for Justin Bahl's molecular epidemiology course at the University of Georgia. It contains novel results and will be included as an appendix to my dissertation (though probably never published formally in an academic journal). Since it will be in my dissertation, this repo contains the code and data to reproduce my results.

In this repo I:
- did some rough antigenic distance calculations (they have now been improved upon and expanded, contact me, Andreas Handel, or read my dissertation for details). This includes cartographic distances based on serological data from a study conducted by Ted M. Ross.
- used the pairwise distance matrices between influenza A strains to build neighbor-joining trees.
- built maximum likelihood trees for the sequences we included.
- compared the phylogenetic ML trees to the neighbor joining trees.

In general, cartography and evolutionary trees seem to tell a different story -- which is not straightforward. We expect antigenic innovations to evolve due to selective pressure in the population, but the genetic evolutionary tree cannot be reconstructed by cartography. Hamming distance, temporal distance, and $p$-Epitope distance beuild much more similar trees to the phylogenetic tree, implying that cartography is again capturing something unique (or maybe there are unanswered questions about cartography).

Anyways, please contact me if you are interested in this work!
