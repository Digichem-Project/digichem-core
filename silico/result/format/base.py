
class Result_dumper():
    """
    ABC for classes that dump results to other formats
    """
    
    def __init__(self, filters, silico_options):
        self.filters = filters
        self.silico_options = silico_options
        
    def dump(self, results):
        data = []
        
        for result in results:
            datum = []            
            for filter in filters:
                datum.append(filter.filter(result))
                
            data.append(datum)
                
        return self.process(data)
    
    def process(self, dumped_results):
        raise NotImplementedError("Implement in subclass")
        
        