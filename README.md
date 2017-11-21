# Analyze Simulation data
# Simulation data is generated for 3 different boundary conditions
# Right end free to slide to get the thermalized length of the ribbon
# Right end clamped at the thermalized length
# Right end clamped at T=0 configuration, Stretched case

# Makefile is present in all the three situations

# SliderData
# Reads traj.gsd file and generates analyzeSlider.log file 
# Use analyzeSlider.log to check for thermalization and thermal equilibrium length of ribbon

# ThermalClampedData
# Reads traj_thermal.gsd file for each run and outputs different files with analyzed data
# analyze.log - System observables at each time step
# width.bin - Ribbon height averaged over width (y) in the last half of simulation for each run
# backbone.bin - Ribbon backbone height averaged in the last half of simulation for each run
# cnode.bin - Time series data for the central node of the backbone for each run
# hgt_prof_real.dat - Avg Height Squared data for heat maps

# StretchedData
# Same as thermal case


