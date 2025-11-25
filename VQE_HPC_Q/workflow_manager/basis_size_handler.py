from molecule_data_extractor import molecule_data_parser as mdp
from pyscf import gto

def orchestration_handler(filepath):
    mol_spec = mdp(filepath)

    atoms = mol_spec["atom"]
    basis = mol_spec["basis"]
    charge = mol_spec["charge"]
    spin = mol_spec["spin"]
    symmetry = mol_spec["symmetry"]

    # Build molecule geometry string
    geom = "; ".join(atoms)

    # Build PySCF Mole
    if symmetry is None:
        mol = gto.M(
            atom=geom,
            basis=basis,
            charge=charge,
            spin=spin,
        )
    else:
        mol = gto.M(
            atom=geom,
            basis=basis,
            charge=charge,
            spin=spin,
            symmetry=symmetry,
        )

    mol.build()

    return mol.nao  # number of basis functions




