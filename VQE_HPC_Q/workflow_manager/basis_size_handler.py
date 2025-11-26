from molecule_data_extractor import molecule_data_parser as mdp
from pyscf import gto

def molecule_builder(filepath):
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

    return mol # number of basis functions


def HPC_local_selector(filepath,threshold = 80):
    mol = molecule_builder(filepath)
    complexity = mol.nao
    if complexity < threshold:
        print("Recommending to use Local system for building Hamiltonian")
        prefer_mode = "local" 
    else:
        print("Recommending to use HPC for building Hamiltonian")
        prefer_mode = "hpc"
    mode = input(f"Choose backend [local/hpc] (recommended: {prefer_mode}): ").strip().lower()
    if mode not in ["local", "hpc"]:
        mode = prefer_mode
    return mode

def SLURM_script_initiator1():
    return None
def Hamiltonian_builder_orchestrator(filepath):

   return None


