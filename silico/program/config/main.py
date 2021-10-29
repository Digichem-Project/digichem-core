from silico.program.base import Program
import silico.program.config.validate


class Config_program(Program):
    """
    The Silico config program.
    """
    
    name = "Silico config helper"
    command = "config"
    description = "setup options for Silico"
    help = "Setup options"
    
    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = sub_parsers_object.add_parser(
            self.command,
            description = self.description,
            epilog = self.epilog,
            help = self.help
        )
        
        subparsers = sub_parser.add_subparsers(dest="prog")
    
        # Create sub parsers for each sub-program. Each will define its own parser.
        silico.program.config.validate.Validate_program.arguments(subparsers)
    
        return sub_parser
    
    def main(self):
        """
        Logic for our program.
        """
        raise NotImplementedError("Not implemented")