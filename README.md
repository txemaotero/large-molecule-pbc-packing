# large-molecule-pbc-packing

This repository provides a python script to generate random molecular
configurations taking into account periodic boundary conditions for large
molecules. [Packmol](http://m3g.iqm.unicamp.br/packmol/home.shtml) software is
required for running this tool so make sure you have it installed and cite it
if you use this script. [Gromacs](http://www.gromacs.org) editconf utility is
also required. If you do not have Gromacs installed a basic installation should
be enough (un linux-based systems run `apt install gromacs`).

## How it works?

The packing process starts putting one of the large molecules in an empty box.
Then the molecule is randomly displaced which may lead in a configuration with
part of the molecule outside the box so periodic boundary conditions are applied
to bring back the atoms to the simulation box. After that, a new large molecule
is put into that box. All these steps are repeated until a simulation box with
all the large molecules is obtained. Finally, additional solvent molecules are
packed inside that box if it is required.
