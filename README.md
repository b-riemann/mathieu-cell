# mathieu-cell

* elementary computations for studying Mathieu unit cells.

* author: Bernard Riemann, [ORCID iD: 0000-0002-5102-9546](https://orcid.org/0000-0002-5102-9546)

* used in the publication

> B. Riemann, "An elementary low-emittance ring design study using Mathieu unit cells", abstract submitted to IPAC20.


## Usage

You can run

    ./make_figures.sh
    
which will generate all required figures. Afterwards, you can call 

    python main.py figurename.pdf
   
to update a given figure.


### internal structure
 
* it seems that that the `microMacroCell` module is the older module i made. At the moment, `main.py` does not use `microMacroCell` in its dependencies.

* the `fourierCell` module has less functionality but is easier to use.