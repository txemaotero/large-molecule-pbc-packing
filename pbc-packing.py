import os
import json
import subprocess
import warnings
from sys import path
from pathlib import Path
from typing import Tuple

import numpy as np

from MDAnalysis import Universe


def find_max_side_box_inside(path):
    """
    Finds the maximum side of the box that contains the molecule in path.
    """
    positions = Universe(path).atoms.positions
    return max(positions.max(0) - positions.min(0))


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
    INP_HEADER = """
tolerance 1.2
avoid_overlap yes
randominitialpoint
filetype pdb
nloop0 10000
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
        self.box_side = np.array(self.input_info['box'])
        self.solvent_box_offsets = np.zeros(6)
        self.long_mol_box_offsets = np.zeros(6)
        self._long_mol_max_side = np.max([find_max_side_box_inside(os.path.join(self.out_dir, path))
                                           for path in self._large_mol_list])
        self._inp_file = None    # This will be used to save the packmol input
        self._init_offsets()

    @property
    def solvent_box_str(self):
        box = self.solvent_box_offsets.copy()
        box[3:] += self.box_side
        return ' '.join([f'{val:g}' for val in box])

    @property
    def long_mol_box_str(self):
        box = self.long_mol_box_offsets.copy()
        box[3:] += self.box_side
        return ' '.join([f'{val:g}' for val in box])

    @property
    def random_long_mol_subbox(self):
        inside = self.long_mol_box_offsets.copy()
        inside[3:] += self.box_side
        for axis, pbc in enumerate(self.input_info['pbc']):
            if pbc:
                continue
            if os.path.exists('./initial.pdb'):
                hist, edges = np.histogram(Universe('./initial.pdb').atoms.positions[:, axis], 
                                           bins=50, range=(inside[axis], inside[axis+3]))
                edges_mid = (edges[1:] + edges[:-1]) / 2
                hist = hist.max() - hist
                hist = hist / hist.sum()
                random_left = max(inside[axis],
                                  np.random.choice(edges_mid, p=hist) - self._long_mol_max_side/2)
                random_left = min(random_left,
                                  inside[axis+3]-self._long_mol_max_side)
            else:
                random_left = np.random.random() * (inside[axis+3] - inside[axis])
                random_left += inside[axis]
            inside[axis] = random_left
            inside[axis+3] = random_left + self._long_mol_max_side
        return ' '.join([f'{val:g}' for val in inside])

    def _init_offsets(self):
        # 1 angstroms on the box sides that are replicated to avoid overlap when
        # pbc applied only for solvent.
        pbcs = np.array(self.input_info['pbc'])
        self.solvent_box_offsets[:3] += pbcs
        self.solvent_box_offsets[3:] -= pbcs
        # 2 angstroms on the box sides that are not replicated to avoid overlap
        # whith interfaces or walls.
        not_pbcs = 2*np.float_(np.logical_not(pbcs))
        self.solvent_box_offsets[:3] += not_pbcs
        self.solvent_box_offsets[3:] -= not_pbcs

        self.long_mol_box_offsets[:3] = not_pbcs
        self.long_mol_box_offsets[3:] = -not_pbcs

    def _unify_input_info(self):
        if 'large_molecules' not in self.input_info:
            raise IOError('At least one large molecule should be specified.')
        self.input_info.setdefault('packmol_executable', 'packmol')
        self.input_info.setdefault('solvent', {})
        if 'box' not in self.input_info or not self.input_info['box']:
            raise IOError('You must specify simulation box shape.')

        self.input_info['box'] = box = self.input_info['box']
        if len(box) == 1:
            box *= 3
        if len(box) != 3:
            raise ValueError('"box" keyword must have 1 or 3 elements.')

        self.input_info.setdefault('pbc', [1, 1, 1])
        if not isinstance(self.input_info['pbc'], list) or not all(i in [0, 1] for i in self.input_info['pbc']):
            raise IOError('Invalid pbc sepecified, use a 3 element list with 1 or 0.')
        self.input_info['pbc'] = tuple(self.input_info['pbc'])

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
        final configuration can be found in the "boxed.pdb" file.

        NOTE: Make sure you have packmol installed in your system.

        Parameters
        ----------
        remove_tmp : bool, Optional
            If True (default) the files from the intermediate steps will be
            removed at the end.

        """
        os.chdir(self.out_dir)
        # Ignore MDAnalysis warnings
        warnings.filterwarnings('ignore')

        # Pack one large
        print(f'Packing large molecules (1/{self.n_large})', end='\r')
        self._inp_file = open('box.inp', 'w+')
        self.write_box_first_inp()
        self._pack_and_fix('box_first', pbc=self.input_info['pbc'])

        # Special case for only one large
        if self._large_mol_list:
            for i in range(self.n_large - 1):
                print(f'\rPacking large molecules ({i+1}/{self.n_large})', end='\r')
                os.remove('final.pdb')
                self.write_box_one_more_large_inp()
                self._pack_and_fix('box_one_more', pbc=self.input_info['pbc'])

        print(f'\rPacking large molecules ({self.n_large}/{self.n_large})', end='\r')
        os.rename('initial.pdb', 'final.pdb')

        print('\nPacking the solvent molecules')
        self.write_box_solvent_inp()
        self._pack_and_fix('box', box_out_basename='boxed', out_basename='boxed',
                           move=False)
        if remove_tmp:
            print('Removing temporary files')
            for _file in Path('.').glob('#*#'):
                _file.unlink()
            for _file in Path('.').glob('*.log'):
                _file.unlink()
            os.remove('final.pdb')
        self._inp_file.close()
        os.remove('box.inp')
        warnings.filterwarnings('default')
        os.chdir('..')

    def _pack_and_fix(self, box_inp_basename: str,
                      box_out_basename: str = 'final',
                      out_basename: str = 'initial', move: bool = True,
                      pbc: Tuple =(1, 1, 1)):
        with open(f'{box_inp_basename}.log', 'w') as box_log:
            self._inp_file.seek(0)
            subprocess.call([self.input_info['packmol_executable']],
                            stdin=self._inp_file, stdout=box_log)
        if not os.path.exists(f'{box_out_basename}.pdb'):
            raise IOError(f'Packmol raises an error. Check the {box_inp_basename}.log file.')
        self._inp_file.seek(0)
        self._inp_file.truncate()
        self.move_and_add_box(f'{box_out_basename}.pdb', f'{out_basename}.pdb', move, pbc)

    def move_and_add_box(self, initial: str, final: str, move: bool = True,
                         pbc: Tuple = (1, 1, 1)):
        """
        Moves molecules a random vector, applies pbcs and writes box dimensions.

        Parameters
        ----------
        initial : str
            Path with the pdb to modify.
        final : str
            Path where the modified pdb will be written.
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
            maximum_displ = self.box_side * self.input_info['pbc']
            universe.atoms.positions += maximum_displ * np.random.random(3)
            universe.atoms.pack_into_box()
            
            # Rotation in caso of walls to avoid acumulation in one side.
            # universe.atoms.positions = self.rotate_box(universe.atoms.positions, np.array([1, 0, 0]), np.pi)

        universe.atoms.write(final)

    def rotate_box(self, positions: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
        """
        Applies a rotation to a given MDAnalysis Universe

        Parameters
        ----------
        positions :  np.ndarray
            Array containing the positions of all atoms in the box.
        axis : np.ndarray
            Axis around which the rotation will be performed.
        angle : np.ndarray
            Angle of rotation.

        Returns
        -------
        np.ndarray
            The rotated positions.
        """
        positions -= np.array([*self.box_side]) * np.array([1/2., 1/2., 1/2.])
        
        rot = get_rotation_matrix(axis, angle)
        positions = rot.dot(positions.T).T

        positions += np.array([*self.box_side]) * np.array([1/2., 1/2., 1/2.])
        return positions

    def write_box_first_inp(self):
        """
        Writes the Packmol input to pack the first large molecule.
        """
        large_mol = self._large_mol_list.pop(0)
        text = self.INP_HEADER + f"""
output final.pdb

structure {large_mol}
    number 1
    resnumbers 2
    inside box  {self.long_mol_box_str}
end structure
        """
        self._inp_file.write(text)

    def write_box_one_more_large_inp(self):
        """
        Writes the Packmol input to pack a large molecule in an existing box.
        """
        large_conf = self._large_mol_list.pop(0)
        text = self.INP_HEADER + f"""
output final.pdb

structure ./initial.pdb
    number 1
    resnumbers 2
    fixed   0.     0.     0.    0.    0.    0.
end structure

structure {large_conf}
    number 1
    resnumbers 2
    constrain_rotation x {np.random.randint(180)}. 90.
    constrain_rotation y {np.random.randint(180)}. 90.
    constrain_rotation z {np.random.randint(180)}. 90.
    inside box  {self.random_long_mol_subbox}
end structure
        """
        self._inp_file.write(text)

    def write_box_solvent_inp(self):
        """
        Writes the Packmol input to pack the solvent in the existing box.
        """
        text = self.INP_HEADER + f"""
output boxed.pdb

structure ./final.pdb
    number 1
    resnumbers 2
    fixed   0.     0.     0.    0.    0.    0.
end structure

        """
        for path, amount in self.input_info['solvent'].items():
            text += f"""
structure {path}
    number {amount}
    resnumbers 2
    inside box  {self.solvent_box_str}
end structure

"""
        self._inp_file.write(text)


def get_rotation_matrix(vector: np.ndarray, angle: float) -> np.ndarray:
    """
    Returns the rotation matrix of a turn around a general axis
    """
    x, y, z = vector/np.sqrt(sum(vector**2))
    cos = np.cos(angle)
    sin = np.sin(angle)

    return np.array(
        [[cos + x**2 * (1-cos), x*y*(1-cos) - z*sin, x*z*(1-cos)+y*sin],
        [y*x*(1-cos)+z*sin, cos+y**2*(1-cos), y*z*(1-cos)-x*sin],
        [z*x*(1-cos)-y*sin, z*y*(1-cos)+x*sin, cos+z**2*(1-cos)]]
        )


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
