import argparse


class Extend_action(argparse.Action):
    """
    An argparse action class to replace the "extend" action for versions of python (3.6 etc) that do not support it natively.
    """

    def __call__(self, parser, namespace, values, option_string = None):
        """
        """
        try:
            attr = getattr(namespace, self.dest)
            
        except AttributeError:
            attr = []
            setattr(namespace, self.dest, attr)
            
        attr.extend(values)