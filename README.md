Aplysia Attractor Analysis
============================

Suite of MATLAB code for the eLife paper
"A spiral attractor drives rhythmic locomotion"
Angela Bruno, Bill Frost & Mark Humphries

[DOI and link to come]

All experimental data is available at:

Please cite the data as:
[CRCNS link and DOI to come]

---
### Summary

We provide here all the code we used:

Figures/

All code to produce all panels of main text and supplemental figures
Note that this code produces more panels than were used in the final paper
It also includes a range of helper functions for plotting

Paths to intermediate results are coded for the locations on the main author's PC:
edit these to point at the paths on your system

Analysis/

All code to produce the analysis results used in all the figures. 
Note that this is our research code, with minimal reformatting, so it does considerably more analyses than are in the paper

Functions/ 

The set of helper functions called by one or more scripts

Data/

The experimental data in the format used in the research code

---
### Dependencies:

Needs:
Chronux toolbox (for spike-train spectra): http://chronux.org/

Uses (relevant functions are included in this toolbox):
Consensus clustering (for detection of ensembles): https://github.com/mdhumphries/SpikeTrainCommunitiesToolBox

Neural ensemble analysis (for classification of ensembles): https://github.com/mdhumphries/NeuralEnsembleAnalysis

--- 

### Getting started

Many scripts were written to separately analyse each of the three stimulation sets (da01, da02, da03). Each will thus analyse the set (e.g. da01) specified by the user across all 10 preparations. 
Some summary scripts load all 30 recordings together 

The top-level function DataSet_Properties.m does the basic analysis used by the majority of other scripts. Run this first.

Some scripts depend on the output of others in their directories. The order of dependencies is:

Classify_attractors/
1. Analyse_by_Recurrence_Points -> Check_All_Real_Eigenvalues.m
2. Analyse_by_Recurrence_Points -> recurrence_point_results_Viz.m -> Are_jumps_perturbations_or_transitions2.m -> Does_Oscillation_Period_Change.m
3. Analyse_by_Recurrence_Points -> recurrence_point_results_Viz.m -> Are_sequential_attractors_the_same.m
4. Analyse_by_Recurrence_Points -> recurrence_point_results_Viz.m -> Are_sequential_attractors_the_same_CommonAxes.m -> Quantifying_Participation_by_PCA_Contribution.m -> RobustnessOfAttractor_to_single_neuron_variability.m
5. Analyse_by_Recurrence_Points -> recurrence_point_results_Viz.m -> Are_sequential_attractors_the_same_CommonAxes.m -> Quantifying_Participation_by_PCA_Contribution.m -> What_is_Participation.m
6. Analyse_by_Recurrence_Points -> recurrence_point_results_Viz.m -> Are_sequential_attractors_the_same_CommonAxes.m -> Initial_Conditions.m
7. RecurrencePlot_LineCounting.m -> Quantify_attractors_by_recurrence_plot_lines.m

Ensembles/
1. StaticEnsembles.m -> all others
2. Analyse_Spike_Train_Properties.m -> Types_of_Ensemble_Across_Dataset.m -> Consistency_of_Types_Across_Programs.m

P10_analysis/

All: Needs output of Dimensionality_of_Projections.m from /Classify_attractors/
1. DecodeP10_statespace.m -> Best_P10_StateSpace_GLMmodels.m

### Bugs

Found an error?
Path wrong or missing? 
Function not found?

Flag all issues by clicking the "Issues" tab at the top of this page (next to "Code"). Raise a new Issue with your error and we will fix as we can. 




