from silico.result.format.table import Table_dumper_ABC
from silico.result.format.filter import Result_filter

import tabulate
import io
import csv


class Table_property_dumper_ABC(Table_dumper_ABC):
    """
    ABC for result dumpers that print tables of one type of result (a table of orbitals, atoms, excited states etc).
    """
    
    @property
    def default_filters(self):
        """
        Return a list of filter objects that will be used by default (when none are specified explicitly).
        """
        return [Result_filter(filter_string, self.silico_options, allow_error = True) for filter_string in [
            "energies:cc",
            "energies:mp",
            "energies:scf",
            "atoms",
            "orbitals",
            "beta_orbitals",
            "excited_states",
            "soc",
            "vibrations",
        ]]
    
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
                    
                    # Prepend the name of the dict to each row.
                    table = [dict([("{}:{}".format(name, col), data) for col, data in row.items()]) for row in table]
                    
                    if len(table) > 0:
                        tables.append(table)
                    
        return tables


class Tabulate_property_dumper(Table_property_dumper_ABC):
    """
    Class for dumping a single result to a text table.
    """
    
    def process(self, dumped_results):
        tables = self.tabulate_data(dumped_results)
        
        return "\n".join([tabulate.tabulate(table, headers = "keys", floatfmt = self.floatfmt) + "\n" for table in tables])


class CSV_property_dumper(Table_property_dumper_ABC):
    """
    Class for dumping a single result to a csv table.
    """
    
    def process(self, dumped_results):
        # Tables to dump.
        tables = self.tabulate_data(dumped_results)
        # String buffer to write to.
        stream = io.StringIO()
        
        for index, table in enumerate(tables):
            headers = self.get_headers(table)
            
            writer = csv.DictWriter(stream, headers)
            
            writer.writeheader()
            for row in table:
                writer.writerow(row)
                
            # If there's at least one more table, separate with newline.
            if index +1 < len(tables):
                stream.write("\n")
        
        return stream.getvalue()
