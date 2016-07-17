#!/bin/bash

f2py -c geodesiclm.pyf -L/data/mark/lib -lgeodesiclm  -llapack -lblas -lgfortran
