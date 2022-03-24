# Delaunay density diangostic (for MSD and grad-MSD rates)
   as described in the paper
   Data-driven geometric scale detection via Delaunay interpolation
   by Andrew Gillette and Eugene Kur
   Version 1.0, March 2022

## Usage

1) Clone the DelaunaySparse repository:

   $ git clone https://github.com/vtopt/DelaunaySparse.git

2) Install the DelaunaySparse software.  

   The software may be able to self-install with the included python bindings.  Try:

   $ cd [directory with DelaunaySparse repository]/DelaunaySparse/python
   $ python example.py 

   If not, try

   $ cd ../src
   $ make

   to troubleshoot the compliation of DELAUNAYSPARSE.  
   MacOS users who encounter the error "ld: library not found for -lSystem" may wish to try

   $ xcode-select --install

   to update the command line tools required for installation.

3) Copy the two python scripts and the .sh file into the DelaunaySparse python subfolder:

   $ cd [directory with DelaunaySparse repository]/DelaunaySparse/python
   $ cp [relevant path prefix]/delauany_density_diagnostic.py .
   $ cp [relevant path prefix]/generate_ddd_figures.py .
   $ cp [relevant path prefix]/run_ddd_trials.sh

4) Run the driver script:

   $ ./run_ddd_trials.sh

   This script will run a total of 100 trials of the delaunay_density_diagnostic.py script, 
      save the results as csv files, then call generate_ddd_figures.py
      to generate a png figure.  Details can be found in the header of the script.

5) If the figure generates correctly, run

   $ python delaunay_density_diagnostic.py --help

   to see the command line options that can be added to the driver script for 
   user-specified experiments.

   