from rdkit import Chem
from rdkit.Chem import Descriptors


def get_descriptors(smiles):
    """
    Get a dictionary of RDKit descriptors from a SMILES string.

    Parameters
    ----------
    smiles : str
        The SMILES string of the chemical of interest

    Returns
    -------
    descriptors : dict
        A collection of molecular descriptors
    
    Notes: Developed with RDKit 2019.03.4, although doc pages listed 2019.03.1
    """

    mol = Chem.MolFromSmiles(smiles)

    descriptors = {}

    # Molecular weight
    descriptors['molwt'] = Descriptors.ExactMolWt(mol)

    return descriptors