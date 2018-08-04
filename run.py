#!/usr/bin/env python

import interfaceff2gro as iff
import os
import sys


Class = ' 1 '
frc_file = 'cvff_interface_v1_5'
word = '.car'

iff.create_lammps_data_file(Class, frc_file, word)


