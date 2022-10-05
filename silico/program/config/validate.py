from silico.program.base import Program


class Validate_program(Program):
    """
    The Silico config validate program.
    """
    
    name = "Silico config validator"
    command = "validate"
    description = "Check the values of loaded config options."
    help = "Validate configs"
    
    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = super().arguments(sub_parsers_object)
        
        sub_parser.add_argument("type", help = "the type of configurables to validate. If none are given, all will be validated.", nargs = "*", choices = ("destinations", "programs", "calculations", "all"), default = "all")
    
        return sub_parser
    
    def main(self):
        """
        Logic for our program.
        """
        if self.args.type == "all":
            self.args.type = ("destinations", "programs", "calculations")
    
        for configurable_type in self.args.type:
            for index, configurable in enumerate(getattr(self.config, configurable_type)):
                configurable.validate()
                self.logger.info("Validated: {}) {}".format(index+1, configurable.meta['name']))
                
        self.logger.info("All configs validated")