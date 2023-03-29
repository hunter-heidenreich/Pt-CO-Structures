import os

from ase.build import add_adsorbate, fcc111, molecule
from ase.constraints import FixAtoms
from ase.io import write
from ase.visualize import view


d_ads = 2.5  # Angs
vacuum = 10.0  # Angs
pt_a = 4.03  # Angs

co = molecule('CO')
c_index = ([a.index for a in co if a.symbol == 'C']).pop()

for site in ['ontop', 'bridge', 'hcp', 'fcc']:
    slab = fcc111('Pt', size=(3, 3, 6), a=pt_a, vacuum=vacuum, periodic=True)
    add_adsorbate(slab, co, d_ads, position=site, offset=1, mol_index=c_index)
    constraint = FixAtoms(indices=[atom.index for atom in slab if atom.position[-1] < 13.0])
    slab.set_constraint(constraint)

    # view(slab)
    os.makedirs(f'out/pt111-p3x3-1_9-{site}', exist_ok=True)
    write(f'out/pt111-p3x3-1_9-{site}/POSCAR', slab, format='vasp')
