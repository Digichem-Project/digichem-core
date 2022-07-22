import yaml

from silico.result.format.base import Result_dumper

class Yaml_dumper(Result_dumper):
    """
    Dump a result set object to yaml text.
    """
    
    def dump(self, results):
        data = []
        
        for result in results:
            datum = []
            # If no filters, just get everything.
            if len(self.filters) == 0:
                datum.append(result.dump())
            
            for filter in self.filters:
                datum.append(filter.filter(result))
                
            data.append(datum)
                
        return self.process(data)
    
    def process(self, dumped_results):
        
        data = yaml.safe_dump_all([inner_data for data in dumped_results for inner_data in data], allow_unicode=True)
        
        # Yaml will insist on appending the document end marker if datum consists
        # only of a simple literal (eg just a string or float).
        # This is annoying.
        # TOOD: A better way of doing this?
        if data.endswith("...\n"):
            data = data[:-4]

        return data
        