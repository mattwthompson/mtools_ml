from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Descriptors3D, GraphDescriptors


def get_features(smiles):
    features = dict()
    features.update(get_descriptors(smiles))

    return features


def get_reference_data(smiles):
    pass

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
    mol = Chem.AddHs(mol)

    Chem.EmbedMolecule(mol, Chem.ETKDG())

    descriptors = {}

    # Starting with simple descriptors:
    # https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html

    # Molecular weight
    descriptors['molwt'] = Descriptors.ExactMolWt(mol)

    # Partial charge metrics
    descriptors['max_abs_partial_charge'] = Descriptors.MaxAbsPartialCharge(mol)
    descriptors['max_partial_charge'] = Descriptors.MaxPartialCharge(mol)
    descriptors['min_abs_partial_charge'] = Descriptors.MinAbsPartialCharge(mol)
    descriptors['min_partial_charge'] = Descriptors.MinPartialCharge(mol)

    # Basic electron counts
    descriptors['num_radical_electrons'] = Descriptors.NumRadicalElectrons(mol)
    descriptors['num_valence_electrons'] = Descriptors.NumValenceElectrons(mol)

    # 3-D descriptors
    # https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors3D.html

    # Calculating these should produce the same result, according to some basic tests
    # descriptors['asphericity'] = rdMolDescriptors.CalcAsphericity(mol)
    # descriptors['eccentricity'] = rdMolDescriptors.CalcEccentricity(mol)
    descriptors['asphericity'] = Descriptors3D.Asphericity(mol)
    descriptors['eccentricity'] = Descriptors3D.Eccentricity(mol)

    descriptors['inertial_shape_factor'] = Descriptors3D.InertialShapeFactor(mol)

    descriptors['radius_of_gyration'] = Descriptors3D.RadiusOfGyration(mol)
    descriptors['spherocity_index'] = Descriptors3D.SpherocityIndex(mol)

    # Graph descriptors
    # https://www.rdkit.org/docs/source/rdkit.Chem.GraphDescriptors.html
    descriptors['balaban_j'] = GraphDescriptors.BalabanJ(mol)
    descriptors['bertz_ct'] = GraphDescriptors.BertzCT(mol)

    descriptors['chi0'] = GraphDescriptors.Chi0(mol)
    descriptors['chi0n'] = GraphDescriptors.Chi0n(mol)
    descriptors['chi0v'] = GraphDescriptors.Chi0v(mol)
    descriptors['chi1'] = GraphDescriptors.Chi1(mol)
    descriptors['chi1n'] = GraphDescriptors.Chi1n(mol)
    descriptors['chi1v'] = GraphDescriptors.Chi1v(mol)
    descriptors['chi2n'] = GraphDescriptors.Chi2n(mol)
    descriptors['chi2v'] = GraphDescriptors.Chi2v(mol)
    descriptors['chi3n'] = GraphDescriptors.Chi3n(mol)
    descriptors['chi3v'] = GraphDescriptors.Chi3v(mol)
    descriptors['chi4n'] = GraphDescriptors.Chi4n(mol)
    descriptors['chi4v'] = GraphDescriptors.Chi4v(mol)

    descriptors['hall_kier_alpha'] = GraphDescriptors.HallKierAlpha(mol)

    descriptors['kappa1'] = GraphDescriptors.Kappa1(mol)
    descriptors['kappa2'] = GraphDescriptors.Kappa2(mol)
    descriptors['kappa3'] = GraphDescriptors.Kappa3(mol)

    # Predicted properties from Wildman and Crippen
    descriptors['log_p'] = Descriptors.MolLogP(mol)
    descriptors['refractivity'] = Descriptors.MolMR(mol)

    return descriptors