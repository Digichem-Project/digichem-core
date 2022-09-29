"""Tests for configurable objects"""

import pytest

from silico.config.configurable.base import Configurable
from silico.config.configurable.options import Options
from silico.config.configurable.option import Option
from silico.exception.configurable import Configurable_option_exception


# Setup our two test classes.
class Parent(Configurable):
    
    scf = Option(default = True, type = bool)
    
    dft = Options(
        grid = Options(
            size = Option(type = int, default = 10)
        )
    )
    
class Intermediate(Parent):
    
    dft = Options(
        functional = Option(default = "B3LYP")
    )
    
class Child(Intermediate):
    
    dft = Options(
        grid = Options(
            grid_name = Option(default = "big")
        )
    )

@pytest.fixture
def child1():
    return Child()

@pytest.fixture
def child2():
    return Child()

@pytest.fixture
def parent():
    return Parent()


def test_basic(parent):
    """Test basic access."""
    
    # Can we retrieve default values?
    assert parent.scf is True
    assert parent.dft['grid']['size'] == 10
    
    # Can we change and retrieve values?
    parent.scf = False
    parent.dft['grid']['size'] = 20
    
    assert parent.scf is False
    assert parent.dft['grid']['size'] == 20


def test_inheritance(child1, child2, parent):
    """Test the inheritance mechanism of Options objects."""
    
    # Before any change, all objects should have equivalent values of grid size.
    assert child1.dft['grid']['size'] == child2.dft['grid']['size'] and child2.dft['grid']['size'] == parent.dft['grid']['size']
    
    # Child objects should have grid name...
    assert child1.dft['grid']['grid_name'] == child2.dft['grid']['grid_name']
    
    # ...but not parent.
    with pytest.raises(Configurable_option_exception):
        parent.dft['grid']['grid_name']
        
    # Make a change to one of the children.
    child1.dft['grid']['size'] = 100
    
    # Check the change, and that the others have not changed
    assert child1.dft['grid']['size'] == 100
    assert child2.dft['grid']['size'] == 10
    assert parent.dft['grid']['size'] == 10
    
    # Check we can also change the parent property.
    parent.dft['grid']['size'] = 1
    
    # Check the change, and that the others have not changed
    assert child1.dft['grid']['size'] == 100
    assert child2.dft['grid']['size'] == 10
    assert parent.dft['grid']['size'] == 1
    
    # Check we can change an actual child property.
    child2.dft['grid']['grid_name'] = "small"
    
    # Check the change, and that the others have not changed
    assert child1.dft['grid']['grid_name'] == "big"
    assert child2.dft['grid']['grid_name'] == "small"
    
    with pytest.raises(Configurable_option_exception):
        nothing = parent.dft['grid']['grid_name']
        
    # Check we can't set a property of parent that doesn't exist.
    with pytest.raises(Configurable_option_exception):
        parent.dft['grid']['grid_name'] = "medium"
        

def test_dumping(child1):
    """Test dumping of the objects to text."""
    
    # Non-explicit dump; everything default.
    assert child1.dump(False) == {}
    
    # Explicit dump.
    assert child1.dump(True) == {
        'scf': True,
        'dft': {
            'functional': "B3LYP",
            'grid': {
                'size': 10,
                'grid_name': "big"
            }
        }
    }