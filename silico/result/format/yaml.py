import yaml

from silico.result.format.base import Result_dumper

class Yaml_dumper(Result_dumper):
    """
    Dump a result set object to yaml text.
    """
    
    def process(self, dumped_results):
        
        #data = yaml.safe_dump_all([inner_data for data in dumped_results for inner_data in data], allow_unicode=True)
        # A list (per result) of dicts.
        pre_data = []
        for dumped_result in dumped_results:
            pre_data.append({})
            for dumped_dict in dumped_result:
                pre_data[-1].update(dumped_dict)
                
        data = yaml.safe_dump_all(pre_data, allow_unicode = True)
        
        # Yaml will insist on appending the document end marker if datum consists
        # only of a simple literal (eg just a string or float).
        # This is annoying.
        # TOOD: A better way of doing this?
        if data.endswith("...\n"):
            data = data[:-4]

        return data
        