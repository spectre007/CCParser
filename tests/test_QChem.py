
from CCParser import QChem
import pytest

@pytest.fixture
def get_instance():

    created_instances = {}

    def _get_instance(class_name):
        if class_name in created_instances.keys():
            return created_instances[class_name]
        else:
            _class = getattr(QChem, class_name)
            print(_class.__name__)
            instance = _class()
            created_instances[class_name] = instance
            return instance
    return _get_instance




class TestGeneral(object):
    
    def test_basis_name(self, get_instance):
        gen = get_instance("General")
        data = ["Requested basis set is cc-pVDZ\n"]
        bas = gen.basis_name(0, data)
        assert bas == "cc-pVDZ"

