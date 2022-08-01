from silico.result.format.filter import Result_filter

class Result_dumper():
    """
    ABC for classes that dump results to other formats
    """
    
    def __init__(self, filters, silico_options, floatfmt = ".2f"):
        self.filters = filters
        self.silico_options = silico_options
        self.floatfmt = floatfmt
    
    @property
    def default_filters(self):
        """
        Return a list of filter objects that will be used by default (when none are specified explicitly).
        """
        # This filter simply returns all data.
        return [Result_filter("", self.silico_options)]
    
    def prettify(self, nested_dict):
        """
        """
        pretty_dict = {}
        
        for name, value in nested_dict.items():
            if isinstance(value, dict):
                value = self.prettify(value)
                
            elif isinstance(value, list) and len(value) > 0 and isinstance(value[0], dict):
                value = [self.prettify(subval) for subval in value]
            
            if isinstance(name, str):
                #TODO: If we want to implement this, need to be aware that it breaks special handling
                # for dicts containing 'value' and 'units' fields.
                #name = ":".join([sub_name[0].capitalize() + sub_name[1:] for sub_name in name.split(":")])
                
                pretty_dict[name.replace("_", " ")] = value
                
            else:
                pretty_dict[name] = value
            
        return pretty_dict
        
    def dump(self, results, pretty = False):
        data = []
        
        filters = self.filters
        if len(filters) == 0:
            # No filter specified, use a default filter.
            filters = self.default_filters
        
        for result in results:
            datum = []            
            for filter in filters:
                datum.append(filter.filter(result))
                
            data.append(datum)
            
        if pretty:
            data = [[self.prettify(datum) for datum in result] for result in data]
                
        return self.process(data)
    
    def write(self, file_name, results, pretty = False):
        with open(file_name, "wt") as out_file:
            out_file.write(self.dump(results, pretty))
    
    def process(self, dumped_results):
        raise NotImplementedError("Implement in subclass")
        
        