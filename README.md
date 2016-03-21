# mlehot
Software for identifying recombination hotspots using LDhot

This software package allows for the estimation of recombination hotspots using different implementations of LDhot (Wall and Stevison 2016; Auton et al. 2012, 2014).  Our new implementation improves the power to detect recombination hotspots over previous implementations.  Briefly, LDhot uses a pairwise composite-likelihood approach based on the work of Hudson (2001).  For a region of interest, LDhot tests whether the central 2 Kb sub-region (of a larger region) has a higher underlying recombination rate than the rest of the region by forming a likelihood-ratio test statistic.  Critical values for this test statistic are estimated using coalescent simulations with a constant recombination rate.  Most of the improved power of our method comes from having a smaller window size, though some increased power is due to a slightly different methodology for calling hotspots.
A gzipped tar file contains all source code, executables, data files and a more extensive readme.