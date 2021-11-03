from silico.config.configurable.option import Option
from silico.submit.base import Method_target


class Extended_basis_set(Method_target):
    """
    Top-level class for basis set targets.
    """
    
    CLASS_HANDLE = ("basis_set",)
    TYPE = Option(default = "basis_set")
    
    basis_set = Option(help = "Basis set data", required = True, type = lambda data: str(data).strip())
    