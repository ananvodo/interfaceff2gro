# interfaceff2gro


@Author: Andres Vodopivec

@Date: 2020-12-08


SCRIPT CLEAN UP AND EASY TO UNDERSTAND EACH STEP. 

VERY IMPORTANT: CORRECTED PCFF_INTERFACE.FRC FILE AND SIGMA(CHARMM) TO SIGMA(GMX) CONVERSION.

interfaceff2gro is a Python script that creates input files for GROMACS using Hendrik Heinz research group Interface forcefield (doi:10.1021/la3038846, [Interface_ff](https://bionanostructures.com/interface-md/)) - so far it only works with CHARMM forcefield (for supported atomtypes please see prm file that is included [here](https://github.com/kolmank/interfaceff2gro/tree/master/forcefield)). 

It uses msi2lmp utility (without the -ignore, and a corrected version of pcff_interface.frc by Krishan Kanhaiya) from [Lammps](https://github.com/lammps/lammps) software package to transform mdf and car files (Materials Studio files) to Lammps data file and [MDAnalysis](https://www.mdanalysis.org/) Python library.

interfacesff2gro creates 3 files: a gro file having coordinates information of your structure, an itp file having topology information and system.top file being a master topology file. 

## Getting Started

The script works with Python 3.6. The easiest way to build Python environment is by using [Conda](https://conda.io/docs/). The source file of msi2lmp are included in [msi2lmp folder](https://github.com/kolmank/interfaceff2gro/tree/master/msi2lmp) so you do not need to download them.

### Prerequisites

Create Python 3.6 environment and install MDAnalysis. The easiest way is to use Conda.

```
conda config --add channels conda-forge
conda create -n envname python=3.6 mdanalysis
source activate envname
```

### Installing

Download the files to you computer:

```
git clone https://github.com/ananvodo/interfaceff2gro.git

```

Enter msi2lmp directory and compile it. According to the authors of the code, the files have been prepared to be compiled using gcc so check if gcc is your default compiler. If you plan to use different compiler, you have to modify the files to make it work. Frc file having atomtypes, bondtypes and angletypes of Interface forcfield has been added to frc_files directory. To compile please use the commands below:

```
cd interfaceff2gro/msi2lmp/src/
make
```

Enter the main directory (/interfaceff2gro) and make interfaceff2gro.py executable:

```
cd ..
cd ..
chmod +x interfaceff2gro.py
```

The script is ready to go!

## Usage 

To use the script, copy car and mdf files to the main directory of the script (you can copy it directly from MODEL_DATABASE directory of Interface forcefield). IMPORTANT! Both files have to have the same name, otherwise the script wont work. Use the command below:

```
./interfaceff2gro.py nameofyourcarfile.car
```

As a result, you should get three files: nameofyourcarfile.gro, nameofyourcarfile.itp and system.top. Examples can be found [here](https://github.com/kolmank/interfaceff2gro/tree/master/examples).

PLEASE be aware that the python script calls msi2lmp.exe without the -ignore flag. If your molecule has missing parameters the python code will fail. If being aware you have missing paramerters and still wish to continue please add -ignore to the very end of the line 167 of the python code.

## Known limitations

* It only works with rectangular simulation boxes.
* Supports only atoms included in prm file, which is [here](https://github.com/kolmank/interfaceff2gro/tree/master/forcefield). According to my knowledge, it should work with silica, metal, cement (only c3a and c3s), clay and hydroxyapatite models that can be found in the database of Interface forcefield.

## Acknowledgments

* Hendrik Heinz research group  that created Interface forcefield (doi:10.1021/la3038846, [Interface_ff](https://bionanostructures.com/interface-md/)).
* [Lammps](https://github.com/lammps/lammps) software package - msi2lmp used in this project is part of it.
* [MDAnalysis](https://www.mdanalysis.org/) Python library.

#### Please feel free to use this program. If you find anything useful in this program or decide to use it for your research, please fork it, click on the star, and add it in your references or citations.
