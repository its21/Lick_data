An example of how I reduce the data and find the velocities for our second night of observations.
I only fit the targets with multiple exposures in this example (just because I saved the targets with 1 exposure with a different nomenclature so I run another program for those, but it’s basically the same procedure, just the name of the files change)

I. DATA REDUCTION: Reduce the data 
in pyraf
data_reduction.py (read the beginning of the program)


II. TEMPLATE FITTING: In data_reduced are the results of data_reduction.py.
(first you need to download the templates! see the README in folder T3500-7500)
in python
%run template_fit_n2.py


III. PLOTTING:
%run plot_best_templates_n2.py (it finishes with an error but that’s because we have not reduced all the target in night 2.. should find an elegant way to end the program)


Note. calibration_errors.py calculates the calibration errors, I didn’t clean the program because the outputs are already saved so you shouldn’t need to run it, but it shows how I calculated them
The rest of the programs are not very useful, I just use them for plotting.