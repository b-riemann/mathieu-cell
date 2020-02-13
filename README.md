# Elementary computations for studying Mathieu unit cells.

* author: Bernard Riemann, [ORCID iD: 0000-0002-5102-9546](https://orcid.org/0000-0002-5102-9546)

* used in the publication

> B. Riemann, "An elementary low emittance lattice study using Mathieu unit cells", manuscript submitted to Phys. Rev. Accel. Beams (2020)

If you publish material using this software, please cite the above reference.

## Requirements

* Python 3.6 or higher (tested on 3.6.3)

* Cython (tested on 0.26.1)

* Scipy (tested on 0.19.1)

* Matplotlib 2.1 (tested on 2.1.0)

## Usage

You can run

    ./make_figures.sh
    
which will generate all required figures. Afterwards, you can call 

    python main.py figurename.pdf
   
to update a given figure.

## Known issues

* the `example.opa` file cannot be read directly into OPA sometimes: this is an encoding problem. Using dos encoding should fix this in principle -- however, the easiest way is duplicating an existing opa file and copy-pasting the contents.

* in the `FourierCell` module, the sign of the `b1` variable is defined opposite to the case in the paper. This is corrected in all PDFs that are generated.

* likewise in the `FourierCell` module, the sign of the `k1` variable is defined opposite to the case in the paper. This is corrected in all PDFs that are generated.


### internal structure

* `main.py` is a python script that generates PDFs (figures) on demand

* the `fourierCell` Cython module (`fourierCell.pyx`) holds the core functionality that is not plotting-related.
 
* The `microMacroCell` Cython module (`microMacroCell.pyx`) stems from my earlier iteration on these scripts. It is still used in a few places, but mainly to compute initial k values for `TuneMap`.

