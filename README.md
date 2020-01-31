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

## Known issues

* in the `FourierCell` module, the sign of the b1 variable is defined opposite to the case in the paper

* the `example.opa` file cannot be read directly into OPA sometimes: this is an encosing problem. using dos encoding should fix this; the easiest way is duplicating an existing opa file and copy-pasting the contents.

### internal structure
 
* The `microMacroCell` module stems from my earlier iteration on these scripts. It is stull used, but only to generate TuneMap k values

* the `fourierCell` module holds the remaining core functionality.


