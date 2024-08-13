Cochlear nucleus is the first stage in the human auditory system that receives auditory input. This repository consists of the biophysically detailed neural network model of bushy cells of ventral cochlear nucleus with implementations of gap junctions and inhibition coming from D-stellate and tuberculoventral cells. More information about the model can be found in accompanying paper.

Once the repository is pulled, open MATLAB, cd to 'BEZ2018_ANF_model' folder and run 'mexANmodel.m' code.

The 'bushy_cell_gap_junctions_main' code is the main code that produces the membrane potential traces of globular and spherical bushy cells, D-stellate cells and tuberculoventral cells. The 'finding_threshold_exct' is used to find the threshold EPSC level depending on the different gap junction strength and the number of cells reside in a fully connected cluster of bushy cells.

Ian C. Bruce (ibruce@ieee.org), Melih Yayli (yaylim@mcmaster.ca) - August 2024