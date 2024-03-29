# Elementary computations for studying Mathieu unit cells.

author: Bernard Riemann, [ORCiD: 0000-0002-5102-9546](https://orcid.org/0000-0002-5102-9546)

used in:

B. Riemann, "Mathieu unit cell as a template for low emittance lattices", [Phys. Rev. Accel. Beams 24, 031601 (2021), doi:10.1103/PhysRevAccelBeams.24.031601](https://doi.org/10.1103/PhysRevAccelBeams.24.031601)

If you publish material using this software, please cite the above reference.

## Requirements

The computations require Python, Cython, and Scipy (see [requirements.txt](requirements.txt)). It is recommended to use a modern Linux distribution, although a setup in Windows should be possible in principle.

## Usage

You can run the shell script `make_figures.sh`. Just enter

    sh make_figures.sh
    
into the terminal, which will generate all required figures as PDFs.

When called the first time, this might take a while due to the required computations. These computations are however stored, thus reruns will be much faster. 

If you change something and want to update a figure, you can try

    python main.py figurename.pdf
   
This works in some cases, but not all, as some calls to `main.py` produce multiple PDFs at once. When in doubt, just run `make_figures.sh` on everything.


## Known issues

* the `example.opa` file cannot be read directly into OPA sometimes: this is an encoding problem. Using dos encoding should fix this in principle -- however, the easiest way is duplicating an existing opa file and copy-pasting the example content into that file.

* in the `FourierCell` sub-module, the sign of the `b1` variable is defined opposite to the case in the paper. This is corrected in all PDFs that are generated.

* likewise in the `FourierCell` sub-module, the sign of the `k1` variable is defined opposite to the case in the paper. This is corrected in all PDFs that are generated.


## Internal structure

* `main.py` is a python script that generates PDFs (figures) on demand

* the `fourierCell` Cython sub-module (`fourierCell.pyx`) holds the core functionality that is not plotting-related.
 
* The `microMacroCell` Cython sub-module (`microMacroCell.pyx`) stems from my earlier iteration on these scripts. It is still used in a few places, but mainly to compute initial k values for `TuneMap`.

