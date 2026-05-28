# Comparing distance-based trees and phylogenetic trees for influenza

Contributors: Savannah L. Miller, Murphy H. John, Amanda L. Skarlupka, Ted M. Ross, Justin Bahl, Andreas Handel

Developing universal influenza vaccines will require improved understanding of how influenza variants differ from each other. We find that temporal distances perform poorly overall, but even sequence distances which match phylogenetic distances well do not match cartographic distances based on actual immune response data.

This was my term project for Justin Bahl's molecular epidemiology course at the University of Georgia. A short draft was published as an Appendix to my dissertation. We are hoping to expand this work into a published manuscript since it provides additional detail about the driving forces behind some of my previous work on antigenic distance.

In this repo I:
- did some rough antigenic distance calculations (they have now been improved upon and expanded, contact me, Andreas Handel, or read my dissertation for details). This includes cartographic distances based on serological data from a study conducted by Ted M. Ross.
- used the pairwise distance matrices between influenza A strains to build neighbor-joining trees.
- built maximum likelihood trees for the sequences we included.
- compared the phylogenetic ML trees to the neighbor joining trees.

In general, cartography and evolutionary trees seem to tell a different story -- which is not straightforward. We expect antigenic innovations to evolve due to selective pressure in the population, but the genetic evolutionary tree cannot be reconstructed by cartography. Hamming distance, temporal distance, and $p$-Epitope distance beuild much more similar trees to the phylogenetic tree, implying that cartography is again capturing something unique (or maybe there are unanswered questions about cartography).

Anyways, please contact me if you are interested in this work!
