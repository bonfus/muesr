
import numpy as np
from ase.atoms import Atoms

class Occupations(Atoms):
    """
    This class extends the Atoms object with fractional atomic occupations.
    It's just a stub for the time being.
    """

    def __init__(self,
                 symbols=None,
                 positions=None,
                 numbers=None,
                 masses=None,
                 magmoms=None,
                 scaled_positions=None,
                 cell=None,
                 pbc=None,
                 occupations=None):
        """
        occupations describe site occupancy factors for each atom in the
        same order as in positions.
        """
        Atoms.__init__(self,
                        symbols,
                        positions,
                        numbers,
                        masses,
                        magmoms,
                        scaled_positions,
                        cell,
                        pbc)
        self._occupations = occupations
        self.generate_groups()

    @classmethod
    def from_atoms(cls, atoms, occupations=None):
        if (occupations is None):
            return cls(symbols=atoms.get_chemical_symbols(),
                        positions=atoms.get_positions(),
                        cell=atoms.cell.array,
                        occupations=np.ones(len(atoms)))
        else:
            return cls(symbols=atoms.get_chemical_symbols(),
                        positions=atoms.get_positions(),
                        cell=atoms.cell.array,
                        occupations=occupations)
        cls.generate_groups()

    def set_groups(self, groups):
        self._occupations_groups = np.array(groups, dtype=np.int32)
    def generate_groups(self):
        self._occupations_groups = np.arange(len(self.scaled_positions), dtype=np.int32) + 1
        self._sites_correlation = np.zeros([len(self.scaled_positions),len(self.scaled_positions)], dtype=np.float)
    def get_occupations(self):
        return self._occupations
    def get_occupations_groups(self):
        return self._occupations_groups
    def get_sites_correlation(self):
        return self._sites_correlation

    def set_occupations(self, value):
        if len(value) == len(self):
            self._occupations = np.array(value, dtype=np.float)
        else:
            raise ValueError("Invalid size. Must be {} elements long".format(len(self)))

    def set_occupations_groups(self, value):
        """
        Specifies which sites are mutually related.
        The groups are specified with labels. Atoms with
        the same labels enter the same group.
        Labels must be specified as positive integers.
        """
        if len(value) == len(self):
            self._occupations_groups = np.array(value, dtype=np.int32)
        else:
            raise ValueError("Invalid size. Must be {} elements long".format(len(self)))

    def create_group(self, a, b):
        """
        Specifies that a and b form a single group.
        a and b indicate corresponding atoms in Atoms class  (starting from 1).
        """

        if a>len(self):
            raise ValueError("Invalid atom number form parameter a. Must be within 1 and {}.".format(len(self)))

        if b>len(self):
            raise ValueError("Invalid atom number form parameter b. Must be within 1 and {}.".format(len(self)))

        current_groups = np.unique(self._occupations_groups)
        new_group_label = np.max(current_groups) + 1
        self._occupations_groups[a-1] = new_group_label
        self._occupations_groups[b-1] = new_group_label
        if self._occupations[a-1] + self._occupations[b-1] > 1.:
            raise warnings.warn("Attention, the total fractional occupation of the newly formed group is higher than 1.")


    def set_sites_correlation(self, value):
        """
        Sets a fake energy according to sites occupations.
        It works like this:
            | a_1    | a_2   | a_3
        ----|--------|-------|------
        a_1 | T_a11  | T_a12 | T_a13
        a_2 |
        a_3 |

        For each atom, the temperature associated to the simultaneous
        occupation of the other sites is evaluated. If larger than a random
        number between 0 and 1 the occupations are re-sampled.

        """
        # This is not a good check
        if len(value) == len(self):
            self._sites_correlation = np.array(value, dtype=np.float)
        else:
            raise ValueError("Invalid size. Must be {} elements long".format(len(self)))


