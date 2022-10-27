"""Tests for methods, loaders and their resolution."""
import pytest

from silico.config.configurable.loader import split_tag_list


@pytest.mark.parametrize("tag_string, tag_list", [
        # Test single items.
        ["gaussian", ["gaussian"]],
        # Simple list.
        ["gaussian: opt", ["gaussian", "opt"]],
        # Check whitespace.
        ["gaussian :   Opt & Freq ", ["gaussian", "Opt & Freq"]],
        # Multiple segments
        ["gaussian :   Opt & Freq:DFT: 6-31G(d,p)", ["gaussian", "Opt & Freq", "DFT", "6-31G(d,p)"]],
        # Namespace.
        ["Turbomole:: Opt : Freq : PBE0PBE : 6-31G(d,p)", ["Turbomole Opt", "Turbomole Freq", "Turbomole PBE0PBE", "Turbomole 6-31G(d,p)"]],
        # Numbers should also be fine.
        ["5246", 5246]
    ])
def test_tag_list(tag_string, tag_list):
    # Split the tag_string.
    split_tag_string = split_tag_list(tag_string)
    
    # Check it's the same as the test string.
    assert split_tag_string == tag_list