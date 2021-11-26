# General imports.
from io import StringIO
from csv import writer

# Silico imports.
from silico.format.summary import *
from silico.format.property import *


class CSV_summary_group_format(Summary_group_format):
    """
    Class for writing calculation results to comma-separated values format (CSV).
    """
    
    def join_results(self, extracted_results):
        """
        Method called to combine a list of extracted results from multiple result sets.
        """
        # Get our table data from our parent.
        # This method will also save column headers for us to self.fieldnames
        table_data =  super().join_results(extracted_results)
                    
        # Now convert to csv.
        return CSV_property_group_format.tabulate(self.fieldnames, table_data)
        
    @classmethod
    def recursive_subclasses(self):
        """
        Recursively get all the subclasses of this class.
        
        Note that text formats don't actually extend from this group format.
        """
        return Summary_format.recursive_subclasses()
    
class CSV_property_group_format(Tabular_property_group_format):
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
        # Get a string buffer.
        buffer = StringIO()
        
        # And a CSV writer.
        csv_writer = writer(buffer)
        
        # Write headers first.
        csv_writer.writerow(fieldnames)
        # Now data.
        csv_writer.writerows(table_data)
        
        # Return the buffer.
        return buffer.getvalue()

    @classmethod
    def recursive_subclasses(self):
        """
        Recursively get all the subclasses of this class.
        """
        return Property_table_format.recursive_subclasses()
    
    
    
    