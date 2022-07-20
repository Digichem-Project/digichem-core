"""
Temp file for writing results to yaml format.

This file will move once the result formatting code has been rewritten.
"""
from silico.format.base import Result_format_group

class To_yaml(Result_format_group):
    
    def __init__(self, *args, **kwargs):
        pass
    
    def extract(self, results):
        return "\n---\n".join([result.to_yaml() for result in results])