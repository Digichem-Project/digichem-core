from silico.program.base import Program
import textwrap
from silico.exception.configurable import Long_tag_path_error


class Show_program(Program):
    """
    A program for showing/listing the available methods.
    """
    
    name = "Show details of the method library"
    command = "show"
    description = "Show the codes and unique IDs of the available methods in the library."
    help = "Show/list methods"
    
    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = super().arguments(sub_parsers_object)
        
        sub_parser.add_argument("path", help = "path to part of a method to show children off", nargs = "*")
        sub_parser.add_argument("-b", "--base", help = "the top-level type of method to start showing from", choices = ["destinations", "programs", "calculations"], default = "destinations")
        sub_parser.add_argument("-d", "--depth", help = "the maximum levels of child methods to show", type = int, default = 1)
    
        return sub_parser
    
    def main(self):
        """
        Logic for our program.
        """
        base = getattr(self.config, self.args.base)
        path = base.path_by_tags(self.args.path)
        if len(self.args.path) == 0:
            header = self.args.base.capitalize() + ":"
        else:
            header = self.get_name(base, path)
        
        body = self.get_body(base, path)
            
        print(header + textwrap.indent(body, "    "))
        
    def get_body(self, base, path, number = 1):
        """
        """
        body = ""
        children = path[-1].get_concrete_children(show_hidden = True)
        for child_path in children:
            full_child_path = path + child_path
            body += "\n" + self.get_name(base, full_child_path, True)
            
            if number < self.args.depth:
                body += textwrap.indent(self.get_body(base, full_child_path, number + 1), "    ")
        
        return body
        
    def get_name(self, base, loader_path, last_only = False):
        """
        """
        if not loader_path[-1].partial:
            code = str(base.index_of_path(loader_path))
            
        else:
            code = None
        
        if not last_only:
            name = " ".join(["'{}'".format(loader.TAG) for loader in loader_path if not loader.pseudo])
        else:
            name = "'{}'".format(loader_path[-1].TAG)
        
        if code is not None:
            return "[{}] ".format(code) + name
        else:
            return name
        
        