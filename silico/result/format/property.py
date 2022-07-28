from silico.result.format.table import Table_dumper_ABC

import tabulate
import io
import csv


class Table_property_dumper_ABC(Table_dumper_ABC):
    """
    ABC for result dumpers that print tables of one type of result (a table of orbitals, atoms, excited states etc).
    """
    
    def tabulate_data(self, dumped_results):
        """
        """
        # Dumped results is a list of lists of dicts.
        # Each item in the first list is one result (one .log file).
        # Each item in the second list is one filtered results.
        # Each dict is the actual data.
        tables = []
        
        for filtered_results in dumped_results:
            for filtered_dict in filtered_results:
                for name, item in filtered_dict.items():
                    # Tabulate the data.
                    table = self.flatten_list(item)
                    
                    table = self.join_lists(table)
                    
                    # Prepend the name of the dict to each row.
                    table = [dict([("{}:{}".format(name, col), data) for col, data in row.items()]) for row in table]
                    
                    tables.append(table)
                    
        return tables
    
    def flatten_list(self, nested_rows):
        # We can only flatten lists, panic if this is not a list (or a dict containing a list).
        if isinstance(nested_rows, dict) and "values" in nested_rows:
            return self.flatten_list(nested_rows['values'])
        
        elif not isinstance(nested_rows, list):
            raise ValueError("Only list properties can be flattened, not '{}'".format(type(nested_rows)))
        
        # A list of dicts (each having the same keys).
        flat_rows = []
        
        for nested_row in nested_rows:
            flat_rows.append(self.flatten(nested_row))
            
        return flat_rows


class Tabulate_property_dumper(Table_property_dumper_ABC):
    """
    Class for dumping a single result to a text table.
    """
    
    def process(self, dumped_results):
        tables = self.tabulate_data(dumped_results)
        
        return "\n".join([tabulate.tabulate(table, headers = "keys") + "\n" for table in tables])
    
class CSV_property_dumper(Table_property_dumper_ABC):
    """
    Class for dumping a single result to a csv table.
    """
    
    def process(self, dumped_results):
        # Tables to dump.
        tables = self.tabulate_data(dumped_results)
        # String buffer to write to.
        stream = io.StringIO()
        
        for table in tables:    
            headers = self.get_headers(table)
            
            writer = csv.DictWriter(stream, headers)
            
            writer.writeheader()
            for row in table:
                writer.writerow(row)
        
        return stream.getvalue()
