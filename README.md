Delaunay Density Diagnostic
----------------
   Version 1.0, March 2022

   This code implements algorithms described in:\
   **Data-driven geometric scale detection via Delaunay interpolation**\
   Andrew Gillette and Eugene Kur, 2022 \
   https://arxiv.org/abs/2203.05685


Usage
----------------

1. Clone the DelaunaySparse repository:
   ~~~~
   git clone https://github.com/vtopt/DelaunaySparse.git
   ~~~~
2. Install the DelaunaySparse software.  \
   The software may be able to self-install with the included python bindings.  Try:
   ~~~~
   cd [directory with DelaunaySparse repository]/DelaunaySparse/python
   python example.py
   ~~~~
   If not, try:
   ~~~~
   cd ../src
   make
   ~~~~
   to troubleshoot the compilation of DelaunaySparse.  
   MacOS users who encounter the error `ld: library not found for -lSystem` may wish to try:
   ~~~~
   xcode-select --install
   ~~~~
   to update the command line tools required for installation.  In addition, it may be necessary to run the following before compiling:
   ~~~~
   export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"
   ~~~~

3. Once DelaunaySparse is successfully installed, copy the `.py` scripts and the `.sh` file from the `ddd` repository into the `DelaunaySparse/python` subfolder:
    ~~~~
   cd [directory with DelaunaySparse repository]/DelaunaySparse/python
   cp [relevant path prefix]/ddd/*.py .
   cp [relevant path prefix]/ddd/*.sh .
   ~~~~
4. Run the driver script:
   ~~~~
   ./run_ddd_trials.sh
   ~~~~
   The above script will run a total of 100 trials of the `delaunay_density_diagnostic.py` script,
      save the results as `.csv` files, then call `generate_ddd_figures.py`
      to generate a `.png` figure.  Details can be found in the header of  `run_ddd_trials.sh`.

5. If the figure generates correctly, run
   ~~~~
   python delaunay_density_diagnostic.py --help
   ~~~~
   to see the command line options that can be added to the driver script for
   user-specified experiments.


Authors
----------------
The Delaunay density diagnostic code was created by Andrew Gillette, gillette7@llnl.gov, with input from Eugene Kur, kur1@llnl.gov.

Citation information
----------------
If you are referring to this code in a publication, please cite the following paper:

Andrew Gillette and Eugene Kur.  *Data-driven geometric scale detection via Delaunay interpolation*.  Submitted.  2022.  LLNL-JRNL-832470.

~~~~
@article{GK2022,
  author = Gillette, Andrew and Kur, Eugene},
  title = {Data-driven geometric scale detection via Delaunay interpolation},
  journal = {Submitted. Preprint at arXiv:2203.05685},
  year = {2022},
}
~~~~




License
----------------

Delaunay density diagnostic is distributed under the terms of the MIT license.

All new contributions must be made under the MIT license.

See [LICENSE](https://github.com/ddd/blob/main/LICENSE) and
[NOTICE](https://github.com/ddd/blob/main/NOTICE) for details.

SPDX-License-Identifier: (MIT)

LLNL-CODE-833036
