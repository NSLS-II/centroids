def test_import():
    from pycentroids import _pycentroids as _pc
    print()
    print("pycentroids c++ module path = {}".format(_pc.__file__))
    print("pycentroids c++ module version = {}".format(_pc.__version__))
