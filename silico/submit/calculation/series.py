from silico.submit.calculation import Calculation_target
from silico.config.configurable.option import Option

class Calculation_series(Calculation_target):
    """
    A pseudo calculation that automatically runs a number of calculations in series.
    """
    
    CLASS_HANDLE = ["Series"]
    
    # Configurable options.
    calculation_IDs = Option(
        "calculations",
        help = "A list of calculations to perform in series",
        #choices = lambda option, configurable: [name for calc in configurable.calculations for name in calc.NAMES],
        required = True,
        type = tuple
    )
    
    def configure(self, CONFIGS, **kwargs):
        """
        Configure this Series calculation.
        """
        self.calculations = CONFIGS
        super().configure(CONFIGS = CONFIGS, **kwargs)
        
    def expand(self):
        """
        Expand this calculation target if it represents multiple real calcs.
        
        For Calculation_series, we return the calcs we represent.
        
        :return: A list of ready-to-go calculation targets.
        """
        #self.calculations = [Calculation_target.from_name_in_configs(calculation_name, configs, silico_options = silico_options, available_basis_sets = available_basis_sets) for calculation_name in calculations]
        return [self.calculations.get_config(ID) for ID in self.calculation_IDs]