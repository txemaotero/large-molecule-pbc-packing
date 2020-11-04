import os
import json
import subprocess
from sys import path
from pathlib import Path

import numpy as np

from MDAnalysis import Universe


class PBCPacking:
    """
    Class to manage the packing of large molecules for MD simulations. 

    Due to the large size of the molecules, this class makes some tricks in
    the packing process to simulate a packing with periodic boundary conditions
    (pbcs). See run_packing method documentation. To perform the packing just
    initialize an instance of the class and then call the `run_packing` method. 

    Parameters
    ----------
    finput : str
        The path to the file with all the information needed to perform the
        packing. See the README for a detailed explanation of this file
        content.
    out_dir : str, Optional
        The path to the directory that will be created to store the outputs. By
        default the class will create a "packing" directory in the current one.

    """

    def __init__(self, finput: str, out_dir: str = 'packing'):
        with open(finput) as _file:
            self.input_info = json.load(_file)
        self._unify_input_info() 

        self.out_dir = os.path.abspath(out_dir)
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self._initial_path = os.path.abspath('.')

        # A list with al the paths to pack with pbcs
        self._large_mol_list = sum([list(np.random.choice(value['paths'],
                                                          value['amount'],
                                                          replace=value['replace']))
                                    for value in self.input_info["large_molecules"].values()],
                                   [])
        self.n_large = len(self._large_mol_list)
        self.box_str = '0. 0. 0. ' + ' '.join([f'{val:g}' for val in self.input_info['box']])
        self.box_side = np.array(self.input_info['box'])

    def _unify_input_info(self):
        if 'large_molecules' not in self.input_info:
            raise IOError('At least one large molecule should be specify.')
        self.input_info.setdefault('packmol_executable', 'packmol')
        self.input_info.setdefault('gromacs_executable', 'gmx')
        self.input_info.setdefault('solvent', {})
        if 'box' not in self.input_info or not self.input_info['box']:
            raise IOError('You must specify simulation box shape.')

        self.input_info['box'] = box = self.input_info['box']
        if len(box) == 1:
            box *= 3
        if len(box) != 3:
            raise ValueError('"box" keyword must have 1 or 3 elements.')

    def run_packing(self, remove_tmp: bool = True):
        """
        Generates the random configuration in the out_dir.

        Packmol is used to put one of the provided large molecules
        configuration in an empty box. Then the molecule in this box is
        displaced a random vector and pbcs are applied to place every atom in
        the simulation box. This new box is used as input in a new packmol
        call and a new large molecule is boxed in the free space of the box.
        A new random displacement is applied now and the process is repeated
        until all the polymer molecules are in the box. Finally, the solvent
        molecules are packed in the box with all the large molecules. The
        final configuration can be found in the "boxed.gro" or "boxed.pdb"
        files.

        NOTE: Make sure you have gromacs and packmol installed in your
        system.

        Parameters
        ----------
        remove_tmp : bool, Optional
            If True (default) the files from the intermediate steps will be
            removed at the end.

        """
        os.chdir(self.out_dir)

        # Pack one large
        self.write_box_first_inp()
        self._pack_and_fix('box_first')

        # Special case for only one large
        if self._large_mol_list:
            self.write_box_one_more_large_inp()

            for _ in range(self.n_large - 1):
                os.remove('final.pdb')
                self._pack_and_fix('box_one_more')

        os.rename('initial.gro', 'final.gro')
        os.rename('initial.pdb', 'final.pdb')

        self.write_box_solvent_inp()
        self._pack_and_fix('box', box_out_basename='boxed', out_basename='boxed',
                           move=False)
        if remove_tmp:
            for _file in Path('.').glob('#*#'):
                _file.unlink()
            for _file in Path('.').glob('*.log'):
                _file.unlink()
            os.remove('final.gro')
            os.remove('final.pdb')
        os.chdir('..')

    def _pack_and_fix(self, box_inp_basename: str,
                      box_out_basename: str = 'final',
                      out_basename: str ='initial', move: bool = True):
        with open(f'{box_inp_basename}.inp') as box_inp, \
             open(f'{box_inp_basename}.log', 'w') as box_log:
            subprocess.call([self.input_info['packmol_executable']],
                            stdin=box_inp, stdout=box_log)
        if not os.path.exists(f'{box_out_basename}.pdb'):
            raise IOError(f'Packmol raises an error. Check the {box_inp_basename}.log file.')
        self._editconf(f'{box_out_basename}.pdb', f'{box_out_basename}.gro')
        self.move_and_add_box(f'{box_out_basename}.gro', f'{out_basename}.gro', move)
        self._editconf(f'{out_basename}.gro', f'{out_basename}.pdb')

    def move_and_add_box(self, initial: str, final: str, move: bool = True):
        """
        Moves molecules a random vector, applies pbcs and writes box dimensions.

        Parameters
        ----------
        initial : str
            Path with the gro to modify.
        final : str
            Path with to write the modified gro.
        move : bool, optional
            If True all the molecules in the system will be displaced a random
            vector and then pbcs will be applied to bring back the atoms to the
            box. Each of the random vector components is selected with a
            homogeneous distribution from 0 to the corresponding box side.
            Setting to False this parameter is useful to write the box
            dimensions in the initial file.

        """
        universe = Universe(initial)
        universe.dimensions = [*self.box_side, 90, 90, 90]
        if move:
            universe.atoms.positions += self.box_side * np.random.random(3)
            universe.atoms.pack_into_box()
        universe.atoms.write(final)

    def _editconf(self, initial: str, final: str):
        command = f'{self.input_info["gromacs_executable"]} editconf -f {initial} -o {final}'.split()
        out = open(os.path.splitext(initial)[0]+'_editconf.log', 'w')
        try:
            subprocess.check_call(command, stdout=out, stderr=out)
        except subprocess.CalledProcessError:
            raise IOError('Error in gromacs editconf')
        finally:
            out.close()

    def write_box_first_inp(self):
        """
        Writes the Packmol input to pack the first large molecule.
        """
        large_mol = self._large_mol_list.pop(0)
        text = f"""
tolerance 0.8
filetype pdb
output final.pdb

structure {large_mol}
    number 1
    inside box  {self.box_str}
end structure
        """
        with open('box_first.inp', 'w') as _file:
            _file.write(text)

    def write_box_one_more_large_inp(self):
        """
        Writes the Packmol input to pack a large molecule in an existing box.
        """
        large_conf = self._large_mol_list.pop(0)
        text = f"""
tolerance 0.8
filetype pdb
output final.pdb

structure ./initial.pdb
    number 1
    fixed   0.     0.     0.    0.    0.    0.
end structure

structure {large_conf}
    number 1
    inside box  {self.box_str}
end structure
        """
        with open('box_one_more.inp', 'w') as _file:
            _file.write(text)

    def write_box_solvent_inp(self):
        """
        Writes the Packmol input to pack the solvent in the existing box.
        """
        text = f"""
tolerance 1.2
filetype pdb
output boxed.pdb

structure ./final.pdb
    number 1
    fixed   0.     0.     0.    0.    0.    0.
end structure

    """
        for path, amount in self.input_info['solvent'].items():
            text += f"""
structure {path}
    number {amount}
    inside box  {self.box_str}
end structure

"""
        with open('box.inp', 'w') as _file:
            _file.write(text)


if __name__ == '__main__':
    import sys
    import shutil

    info_file = sys.argv[1]
    try:
        out_dir = sys.argv[2]
    except IndexError:
        out_dir ='packing'

    packing = PBCPacking(info_file, out_dir=out_dir)
    packing.run_packing()
