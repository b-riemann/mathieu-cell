#!/bin/sh
python setupMicroMacroCell.py build_ext --inplace
python setupFourierCell.py build_ext --inplace
python main.py islands.pdf necktie.pdf chroma.pdf F.pdf example.pdf # jacow
python main.py Islands.pdf sextuVals.pdf G.pdf b1scanA b1scanB Example.pdf extendedExample.pdf H.pdf # + prab
