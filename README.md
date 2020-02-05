# mathieu-cell

* elementary computations for studying Mathieu unit cells.

* author: Bernard Riemann, [ORCID iD: 0000-0002-5102-9546](https://orcid.org/0000-0002-5102-9546)

* used in the publications

> B. Riemann, "An elementary low emittance lattice study using Mathieu unit cells", manuscript submitted to Phys. Rev. Accel. Beams

> B. Riemann, "An elementary low-emittance ring design study using Mathieu unit cells", abstract submitted to IPAC20.

## Usage

You can run

    ./make_figures.sh
    
which will generate all required figures. Afterwards, you can call 

    python main.py figurename.pdf
   
to update a given figure.

## Known issues

* the `example.opa` file cannot be read directly into OPA sometimes: this is an encoding problem. using dos encoding should fix this; the easiest way is duplicating an existing opa file and copy-pasting the contents.

* in the `FourierCell` module, the sign of the b1 variable is defined opposite to the case in the paper. This is corrected in all PDFs that are generated.

* likewise in the `FourierCell` module, the sign of the k1 variable is defined opposite to the case in the paper. This is corrected in all PDFs that are generated.


### internal structure

* `main.py` is a python script that generates PDFs (figures) on demand

* the `fourierCell` Cython module (`fourierCell.pyx`) holds the core functionality that is not plotting-related.
 
* The `microMacroCell` Cython module (`microMacroCell.pyx`) stems from my earlier iteration on these scripts. It is still used in a few places, but mostly only initial k values for TuneMap.

