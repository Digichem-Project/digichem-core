# General imports.
import tabulate

# Silico imports.
from silico.program.base import Program
from silico.exception.base import Silico_exception


class Status_program(Program):
    """
    The Silico method status program.
    """
    
    name = "Method status"
    command = "status"
    description = "check status of known submission methods"
    help = "Check status"
    
    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = super().arguments(sub_parsers_object)
        
        sub_parser.add_argument("methods", help = "Selected methods to show status for", nargs = "*", default = ())
    
        return sub_parser
    
    def main(self):
        """
        Logic for our program.
        """
        # Load our calculation definitions.
        try:
            known_methods = self.config.methods
            
        except Exception:
            raise Silico_exception("Failed to load calculation methods")
        
        # Now get the one's we've been asked to list.
        if len(self.args.methods) == 0:
            # None specified, use all.
            methods = known_methods
            
        else:
            methods = type(known_methods)([known_methods.get_config(method_id) for method_id in self.args.methods])
            
        # Build a table of status to show.
        statuses = []
        for method in methods:
            try:
                status = method.status
            except NotImplementedError:
                # No status for this method.
                status = "N/A (status not available)"
            except Exception:
                status = "Error retrieving status"
                self.logger.error("Failed to fetch status information for method '{}'".format(method.NAME), exc_info = True)
                
            statuses.append((method.ID, method.NAME, status))
                
        # Print with tabulate.
        print(tabulate.tabulate(statuses, ("ID", "Name", "Status")))
        