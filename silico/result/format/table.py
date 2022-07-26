from builtins import isinstance
import csv
import io

from silico.result.format.base import Result_dumper
import itertools
import tabulate


class Table_dumper_ABC(Result_dumper):
    
    def tabulate_data(self, dumped_results):
        """
        Take a list of lists of dicts of dumped results and return a list of dicts in flattened 'table' format.
        
        :returns: A tuple of the table headers and data.
        """
        # First, flatten our data ready for writing.
        flat_results = []
        
        for dumped_result in dumped_results:
            flat_results.append(self.flatten_results(dumped_result))
            
        # Look for any remaining lists and convert to a text form.
        for flat_result in flat_results:
            for name, item in flat_result.items():
                if isinstance(item, list) or isinstance(item, tuple):
                    flat_result[name] = ", ".join([str(subitem) for subitem in item])
                    
        # Get headers.
        headers = list(itertools.chain(*[list(flat_result.keys()) for flat_result in flat_results]))
        # Make unique.
        headers = list(dict.fromkeys(headers))
        
        return headers, flat_results


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
                    flat_result[name + " / " + item['units']] = item['value']
                    item.pop("units")
                    item.pop("value")
                
                flat_child = self.flatten(item)
                flat_result.update({"{}:{}".format(name, child_name): value for child_name, value in flat_child.items()})
            
            else:
                flat_result[name] = item
                
        return flat_result
            
    def flatten_results(self, dumped_result):
        """
        Flatten a recursive list and/or dict so it can be printed as a single, non-nested, row.
        """
        flat_result = {}
        
        for filtered in dumped_result:
            flat_result.update(self.flatten(filtered))
                
        return flat_result


class Tabulate_dumper(Table_dumper_ABC):
    """
    """
    
    def process(self, dumped_results):
#         # First, flatten our data ready for writing.
#         flat_results = []
#         
#         for dumped_result in dumped_results:
#             flat_results.append(self.flatten_results(dumped_result))
#             
#         headers = list(itertools.chain(*[list(flat_result.keys()) for flat_result in flat_results]))
#         
#         headers = list(dict.fromkeys(headers))
        
        
        
        headers, flat_results = self.tabulate_data(dumped_results)
        
        return tabulate.tabulate(flat_results, headers = "keys") + "\n"


class CSV_dumper(Table_dumper_ABC):
    """
    """
        
    def process(self, dumped_results):
#         # First, flatten our data ready for writing.
#         flat_results = []
#         
#         for dumped_result in dumped_results:
#             flat_results.append(self.flatten_results(dumped_result))
#             
#         headers = list(itertools.chain(*[list(flat_result.keys()) for flat_result in flat_results]))
#         
#         headers = list(dict.fromkeys(headers))
        
        headers, flat_results = self.tabulate_data(dumped_results)
        
        
        stream = io.StringIO()
        writer = csv.DictWriter(stream, headers)
        
        writer.writeheader()
        for flat_result in flat_results:
            writer.writerow(flat_result)
            
        return stream.getvalue()
