import textwrap

from silico.program.base import Program


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
        
        print(self.args.base.capitalize() + ":")
        # If we are starting from a node, add a bit of extra space to indent properly.
        if len(path) == 1 and path[0] == base:
            indent = 1
        else:
            indent = 2
            print("    " + self.get_title(path, False))
        
        tree = self.build_tree(path, self.args.depth + indent, cur_depth = 1 + indent)
        
        self.print_tree(tree)
        
    def get_title(self, path, final = True):
        # If our node has a single loader at the end, get a code.
        if not path[-1].partial:
            code = "[{}]".format(path[0].index_of_path(path))
        else:
            code = None
        
        if final:
            name = "'{}'".format(path[-1].TAG)
        else:
            name = " ".join(["'{}'".format(loader.TAG) for loader in path if not loader.partial])
        
        if code is not None:
            title = code + " " + name
        else:
            title = name
            
        return title
        
    def print_tree(self, tree):
        """
        """
        for node in tree:
            depth = node[0]
            path = node[1]
            children = node[2]
            
            print(textwrap.indent(self.get_title(path), "    " * (depth-1)))
            
            if len(children) > 0:
                self.print_tree(children)
        
        
    def build_tree(self, path, max_depth, *, cur_depth = 1):
        """
        """
        tree = []        
        children = path[-1].get_concrete_children(True)
        
        # Don't save the parent path if it ends with a single loader (because we are now a different type)
        if not path[-1].partial:
            prev_path = []
        else:
            prev_path = path
        
        for child_path in children:                
                
            # Get grand children if we've been asked to.
            if cur_depth < max_depth:
                grand_children = self.build_tree(prev_path + child_path, max_depth, cur_depth = cur_depth +1)
                
            else:
                grand_children = []
            
            tree.append([cur_depth, prev_path + child_path, grand_children])
        
        return tree
        
        
        
        
        
        
        
        
        
        
        
        
        