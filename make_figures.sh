#!/bin/sh
python setupMicroMacroCell.py build_ext --inplace
python setupFourierCell.py build_ext --inplace
python main.py islands.pdf necktie.pdf tunemap-check.pdf
