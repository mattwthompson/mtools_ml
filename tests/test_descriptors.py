from mtools_ml.descriptors import get_descriptors

def test_descriptors():
    """Just make sure it can return something."""
    descriptors = get_descriptors('CCC')

    assert isinstance(descriptors, dict)
    assert 'molwt' in descriptors.keys()