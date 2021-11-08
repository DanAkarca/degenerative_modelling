# degenerative_modelling
A simple toolbox to run degenerative network models

The algorithm works by initialising a network, before pruning the connections away according to some degenerative mechanism.
It is underpinned by some cost of pruning alongside some value to pruning the edge away.
The degenerative rules are akin to those described by Betzel et al. 2016 (Neuroimage) and Akarca et al 2021 (Nature Communications).
They are parameterised by scalars, eta and gamma which determines the extent to which the components of the degeneration influence pruning. 
