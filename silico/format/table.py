# General imports.
from tabulate import tabulate

# Silico imports.
from silico.format.summary import *
from silico.format.property import *


class Table_summary_group_format(Summary_group_format):
    """
    Class for writing calculation summaries to text-table format, suitable for viewing in a terminal.
    """
    
    def join_results(self, extracted_results):
        """
        Method called to combine a list of extracted results from multiple result sets.
        """
        # Get our table data from our parent.
        # This method will also save column headers for us to self.fieldnames
        data_rows =  super().join_results(extracted_results)
            
        # Now convert to a texty-table (thanks tabulate).
        return "{}\n".format(tabulate(data_rows, headers = self.fieldnames, numalign = "center", stralign = "center"))
                    
        
    @classmethod
    def recursive_subclasses(self):
        """
        Recursively get all the subclasses of this class.
        """
        return Summary_format.recursive_subclasses()
    
class Property_table_group_format(Tabular_property_group_format):
    """
    Group format for property tables in text format (good for reading orbitals, atoms etc in text files and terminals).
    """
    
    @classmethod
    def tabulate(self, fieldnames, table_data):
        """
        Convert some raw table data into our desired format (texty-table in our case).
        
        :param fieldnames: List of fieldnames/table headers.
        :param table_data: List of table rows.
        """
        return "{}\n".format(tabulate(table_data, headers = fieldnames, numalign = "center", stralign = "center", floatfmt=".4f"))

    @classmethod
    def recursive_subclasses(self):
        """
        Recursively get all the subclasses of this class.
        """
        return Property_table_format.recursive_subclasses()
    
    
