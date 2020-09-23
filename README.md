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

## Dependencies

As it was mentioned, [Packmol](http://m3g.iqm.unicamp.br/packmol/home.shtml)
and [Gromacs](http://www.gromacs.org) are required as well as
[python](https://www.python.org) 3.6 or grater and the
[numpy](https://numpy.org) python package. One may ask why Gromacs is needed.
The answer is that we use the `SystemGro` utility from `compyna` module to
move the molecules and apply pbcs to move the atoms to the simulation box.

## How to use it

A basic usage of this tool would be to run in a terminal:

```
python pbc-packing.py input_info.json out_dir_name
```

The first argument (`input_info.json`) is the path to the file with all the
necessary information to perform the packing and the second one is a
directory which will be created to store all the outputs. The last is
optional and by default it will be set to `packing`. To write the
`inpunt_info.json` file see the structure of the
[example](example_input.json) that can be found in this repository. In the
next lines the meaning of each line in this example is commented:

```js
{
    // Specify the paths to the packmol and gromacs executables, this is, the
    // command that you use to run the programs
    "packmol_executable": "packmol",
    "gromacs_executable": "gmx",
    // "box" is a list of floats corresponding to the dimension of the box in
    // the x, y and z directions. If the list has just one value a cubic box
    // will be created.
    "box": [60],
    // "large_molecules" is a dictionary specifying all the kind of large
    // molecules that you want to put into the simulation box using pbcs.
    "large_molecules": {
        // For each kind of molecule you have to specify a name (PEO in this
        // example).
        "PEO": {
            // If "replace" is false the large molecules that will be put into
            // the box will be selected from the "paths" list without repetition
            // otherwise, if it is true the same molecule can be randomly
            // selected multiple times.
            "replace": true,
            // "amount" is the number of large molecule of this kind that you
            // want to put into the box.
            "amount": 20,
            // "paths" is a list with all the available molecular configurations
            // of this kind of large molecule that can be randomly selected to
            // put into the box. NOTE: It is highly recomended to write absolute
            // paths to the files. Also note that if "replace" is false "paths"
            // must have at least "amount" different elements. If you only have
            // one molecular configuration of this kind of large molecule you
            // can set "replace" to true and just specify one path in "paths".
            "paths": [
                "../molecules/large/peo_0.pdb",
                "../molecules/large/peo_1.pdb",
                "../molecules/large/peo_2.pdb",
                "../molecules/large/peo_3.pdb",
                "../molecules/large/peo_4.pdb",
                "../molecules/large/peo_5.pdb",
                "../molecules/large/peo_6.pdb",
                "../molecules/large/peo_7.pdb",
                "../molecules/large/peo_8.pdb",
                "../molecules/large/peo_9.pdb",
                "../molecules/large/peo_10.pdb",
                "../molecules/large/peo_11.pdb",
                "../molecules/large/peo_12.pdb",
                "../molecules/large/peo_13.pdb",
                "../molecules/large/peo_14.pdb",
                "../molecules/large/peo_15.pdb",
                "../molecules/large/peo_16.pdb",
                "../molecules/large/peo_17.pdb",
                "../molecules/large/peo_18.pdb",
                "../molecules/large/peo_19.pdb"
            ]
        }
    },
    // At the end of the large molecule packing process you can optionally
    // solvate the obtained box with other molecules which should be specified
    // in the "solvent" section. Note all this solvent molecules will be put in
    // the box at the same time without any kind of additional displacement.
    "solvent": {
        // For each solvent molecule specify the path to the pdb file with the
        // molecular configuration and the amount of them you want in to be in
        // the box. Again, we recomend to use absolute paths.
        "../molecules/ace.pdb": 200,
        "../molecules/ea.pdb": 200
    }
}
```

The script simply initializes an instance of the `PBCPacking` class with the
input parameters and then call its `run_packing` method. This class is also
documented so if you have any doubt do not hesitate to check it out. 


## Limitations

In the current version, this tool only allows to generate orthoedric simulation
boxes. On the other hand, although a box without solvent molecules can be
generated at least one large molecule is required.
