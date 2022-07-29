from silico.result.format.filter import Result_filter

from silico.result.format.base import Result_dumper
import itertools
import tabulate
import csv
import io


class Table_dumper_ABC(Result_dumper):
    """
    ABC for classes that dump result sets to tabulated formats.
    """
    
    @property
    def default_filters(self):
        """
        Return a list of filter objects that will be used by default (when none are specified explicitly).
        """
        return [Result_filter(filter_string, self.silico_options, allow_error = True) for filter_string in [
            "metadata",
            "ground_state",
            "energies",
            "atoms",
            "orbitals",
            "orbitals:HOMO",
            "orbitals:LUMO",
            "beta_orbitals",
            "beta_orbitals:HOMO",
            "beta_orbitals:LUMO",
            "pdm",
            "excited_states",
            "excited_states:S(1)",
            "excited_states:T(1)",
            # NOTE: SOC by itself currently dumps nothing.
            "soc",
            "soc:S(0),T(1)",
            "soc:S(1),T(1)",
            "vibrations",
            "emission"
        ]]
    
    @classmethod
    def get_headers(self, table):
        """
        Return a list of unique keys from a list of dicts.
        """
        # Get headers.
        headers = list(itertools.chain(*[list(flat_result.keys()) for flat_result in table]))
        # Make unique.
        return list(dict.fromkeys(headers))
    
    def join_lists(self, table):
        """
        Convert any table cells that contain lists to a nicer text form.
        """
        for row in table:
            for col, cell in row.items():
                if isinstance(cell, list) or isinstance(cell, tuple):
                    row[col] = ", ".join([str(subitem) for subitem in cell])
                    
        return table
    
    def tabulate_data(self, dumped_results):
        """
        Take a list of lists of dicts of dumped results and return a list of dicts in flattened 'table' format.
        
        :returns: A tuple of the table headers and data.
        """
        # First, flatten our data ready for writing.
        table = []
        
        for dumped_result in dumped_results:
            table.append(self.flatten_results(dumped_result))
            
        # Look for any remaining lists and convert to a text form.
        table = self.join_lists(table)
        
        return table

    def flatten(self, dumped_result):
        """
        Flatten a recursive list and/or dict so it can be printed as a single, non-nested, row.
        """
        flat_result = {}
        
        for name, item in dumped_result.items():
            if (isinstance(item, list) or isinstance(item, tuple)) and (len(item) == 0 or isinstance(item[0], dict)):
                    # Don't include lists of dicts (they are too big to flatten).
                    continue
            
            if isinstance(item, dict):
                if "value" in item and "units" in item:
                    # TODO: What if units is None? Skip?
                    flat_result["{} / {}".format(name, item['units'])] = item['value']
                    item.pop("units")
                    item.pop("value")
                
                flat_child = self.flatten(item)
                flat_result.update({"{}:{}".format(name, child_name): value for child_name, value in flat_child.items()})
            
            else:
                flat_result[name] = item
                
        return flat_result
            
    def flatten_results(self, dumped_result):
        """
        Flatten a list of dicts so it can be printed as a single, non-nested, row.
        """
        flat_result = {}
        
        # Dumped result is a list of dicts, each item in the list is one filtered result.
        for filtered in dumped_result:
            flat_result.update(self.flatten(filtered))
                
        return flat_result


class Tabulate_dumper(Table_dumper_ABC):
    """
    Class for dumping a result set to a text table.
    """
    
    def process(self, dumped_results):
        flat_results = self.tabulate_data(dumped_results)
        
        return tabulate.tabulate(flat_results, headers = "keys") + "\n"


class CSV_dumper(Table_dumper_ABC):
    """
    Class for dumping a result set to a CSV table.
    """
        
    def process(self, dumped_results):
        flat_results = self.tabulate_data(dumped_results)
        headers = self.get_headers(flat_results)
        
        stream = io.StringIO()
        writer = csv.DictWriter(stream, headers)
        
        writer.writeheader()
        for flat_result in flat_results:
            writer.writerow(flat_result)
            
        return stream.getvalue()
