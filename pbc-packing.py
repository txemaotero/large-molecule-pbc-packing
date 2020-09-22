import os
import json
import subprocess
from sys import path
from pathlib import Path

import numpy as np

from compyna import SystemGro


class PBCPacking:
    """
    Class to manage simulations with large molecules like polymers.

    Due to the large size of the molecules, this class makes some tricks in
    the packing process to simulate a packing with pbcs (see run_box method
    documentation). 

    """

    def __init__(self, finput: str, out_dir: str = 'packing'):
        with open(finput) as _file:
            self.inpunt_info = json.load(_file)
        self._unify_input_info() 

        self.out_dir = os.path.abspath(out_dir)
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self._initial_path = os.path.abspath('.')

        # A list with al the paths to pack with pbcs
        self._large_mol_list = sum([list(np.random.choice(value['paths'],
                                                          value['amount'],
                                                          replace=value['replace']))
                                    for value in self.inpunt_info["large_molecules"].values()],
                                   [])
        self.n_large = len(self._large_mol_list)
        self.box_str = '0. 0. 0. ' + ' '.join([f'{val:g}' for val in self.inpunt_info['box']])
        self.box_side = np.array(self.inpunt_info['box'])

    def _unify_input_info(self):
        self.inpunt_info.setdefault('packmol_executable', 'packmol')
        self.inpunt_info.setdefault('gromacs_executable', 'gmx')
        self.inpunt_info.setdefault('large_molecules', {})
        self.inpunt_info.setdefault('solvent', {})
        if 'box' not in self.inpunt_info or not self.inpunt_info['box']:
            raise IOError('You must specify simulation box shape.')

        self.inpunt_info['box'] = box = self.inpunt_info['box']
        if len(box) == 1:
            box *= 3
        if len(box) != 3:
            raise ValueError('"box" keyword must have 1 or 3 elements.')

    def run_packing(self, remove_tmp: bool = True):
        """
        Generates the inital configuration using Packmol in a box subdir.

        Packmol is used to put one of the provided polymer configuration in
        an empyt box. Then the molecule in this box is displaced a random
        vector and pbcs are applied to place every atom in the simulation
        box. This new box is used as input in a new packmol and a new polymer
        is boxed in the free space of the box. The same displacemnt is
        applied now and the process is repeated untill all the polymer
        molecules are in the box. Finally, the solvent molecules are packed
        in the box with all the polimers.

        The final configuration is in the "boxed.gro" file.

        NOTE: packmol comand should be available in the system.
        """
        os.chdir(self.out_dir)

        # Pack one poly
        self.write_box_first_inp()
        self._pack_and_fix('box_first')

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
            subprocess.call([self.inpunt_info['packmol_executable']],
                            stdin=box_inp, stdout=box_log)
        if not os.path.exists(f'{box_out_basename}.pdb'):
            raise IOError(f'Packmol raises an error. Check the {box_inp_basename}.log file.')
        self._editconf(f'{box_out_basename}.pdb', f'{box_out_basename}.gro')
        self.move_and_add_box(f'{box_out_basename}.gro', f'{out_basename}.gro', move)
        self._editconf(f'{out_basename}.gro', f'{out_basename}.pdb')

    def move_and_add_box(self, initial: str, final: str, move: bool = True):
        """
        Displaces the simulation box a random vector, applies pbc and writes
        box dimensions.

        Parameters
        ----------
        initial: str
            Path with the gro to modify.
        final: str
            Path with to write the modified gro.

        """
        system = SystemGro(initial)
        l_box_nm = self.box_side * 0.1
        np.fill_diagonal(system.box_vectors, l_box_nm)
        if move:
            disp = l_box_nm * np.random.random(3)
            system.move_molecules(disp)
            system.inside_box()
        system.write_gro(final)

    def _editconf(self, initial: str, final: str):
        command = f'{self.inpunt_info["gromacs_executable"]} editconf -f {initial} -o {final}'.split()
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
        Writes the Packmol input to pack ths solvent in the existing box.
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
        for path, amount in self.inpunt_info['solvent'].items():
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

    shutil.rmtree(out_dir)
    packing = PBCPacking(info_file, out_dir=out_dir)
    packing.run_packing()
