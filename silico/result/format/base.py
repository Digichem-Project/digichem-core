
class Result_dumper():
    """
    ABC for classes that dump results to other formats
    """
    
    def __init__(self, filters):
        self.filters = filters
        
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
        raise NotImplementedError("Implement in subclass")
        
        