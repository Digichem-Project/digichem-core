import yaml


class Yaml_dumper():
    """
    Dump a result set object to yaml text.
    """
    
    def __init__(self, filters):
        self.filters = filters
        
    def dump(self, result):
        data = []
        
        # If no filters, just get everything.
        if len(self.filters) == 0:
            data.append(self.process(result.dump()))
        
        for filter in self.filters:
            data.append(self.process(filter.filter(result)))
            
        return "...\n".join(data)
    
    def process(self, filtered_data):
        data = yaml.safe_dump(filtered_data, allow_unicode=True)
        
        # Yaml will insist on appending the document end marker if datum consists
        # only of a simple literal (eg just a string or float).
        # This is annoying.
        # TOOD: A better way of doing this?
        if data.endswith("...\n"):
            data = data[:-4]
        return data
        
        