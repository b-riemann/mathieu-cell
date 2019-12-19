#!/bin/sh
python setupMicroMacroCell.py build_ext --inplace
python setupFourierCell.py build_ext --inplace
python main.py islands.pdf necktie.pdf chroma.pdf F.pdf example.pdf # jacow
python main.py sextuVals.pdf # prab
