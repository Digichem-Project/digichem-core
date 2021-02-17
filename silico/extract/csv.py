from io import StringIO
from csv import writer
from silico.extract.summary import *
from silico.extract.long import *

class CSV_summary_group_extractor(Summary_group_extractor):
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
        return Long_CSV_group_extractor.tabulate(self.fieldnames, table_data)
        
    @classmethod
    def recursive_subclasses(self):
        """
        Recursively get all the subclasses of this class.
        
        Note that text extractors don't actually extend from this group extractor.
        """
        return Summary_extractor.recursive_subclasses()
    
class Long_CSV_group_extractor(Long_tabular_group_extractor):
    """
    Group extractor for long tables in text format (good for reading orbitals, atoms etc in text files and terminals).
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
        return Long_table_extractor.recursive_subclasses()
    
    
    
    