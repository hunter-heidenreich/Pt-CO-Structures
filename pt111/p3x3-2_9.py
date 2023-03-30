import os

from itertools import product

from ase.build import add_adsorbate, fcc111, molecule
from ase.constraints import FixAtoms
from ase.io import write
from ase.visualize import view

site2sym = {
    'ontop': 'T',
    'bridge': 'B',
    'fcc': 'F',
    'hcp': 'H'
}

d_ads = 2.5  # Angs
vacuum = 10.0  # Angs
pt_a = 4.03  # Angs

co = molecule('CO')
c_index = ([a.index for a in co if a.symbol == 'C']).pop()

sites = ['ontop', 'bridge', 'hcp', 'fcc']

seen = set()
for site1, site2 in product(*(sites, sites)):
    s_ = tuple(sorted((site1, site2)))
    if s_ in seen:
        continue
    seen.add(s_)

    slab = fcc111('Pt', size=(3, 3, 6), a=pt_a, vacuum=vacuum, periodic=True)

    add_adsorbate(slab, co, d_ads, position=site1, offset=0, mol_index=c_index)
    add_adsorbate(slab, co, d_ads, position=site2, offset=1 if site1 == site2 else (1, 0), mol_index=c_index)

    constraint = FixAtoms(indices=[atom.index for atom in slab if atom.position[-1] < 13.0])
    slab.set_constraint(constraint)

    # view(slab)
    s = site2sym[site1] + site2sym[site2]
    os.makedirs(f'out/pt111-p3x3-2_9-{s}', exist_ok=True)
    write(f'out/pt111-p3x3-2_9-{s}/POSCAR', slab, format='vasp')
