from silico.result.format.filter import Result_filter

class Result_dumper():
    """
    ABC for classes that dump results to other formats
    """
    
    def __init__(self, filters, silico_options):
        self.filters = filters
        self.silico_options = silico_options
    
    @property
    def default_filters(self):
        """
        Return a list of filter objects that will be used by default (when none are specified explicitly).
        """
        # This filter simply returns all data.
        return [Result_filter("", self.silico_options)]
        
    def dump(self, results):
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
                
        return self.process(data)
    
    def process(self, dumped_results):
        raise NotImplementedError("Implement in subclass")
        
        