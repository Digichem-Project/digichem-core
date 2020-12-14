from tabulate import tabulate
from silico.extract.summary import *
from silico.extract.long import *
#from silico.extract.long import Long_table_group_extractor

class Table_summary_group_extractor(Summary_group_extractor):
	"""
	Class for writing calculation summaries to text-table format, suitable for viewing in a terminal.
	"""
	
	def join_results(self, extracted_results):
		"""
		Method called to combine a list of extracted results from multiple result sets.
		"""
		# Get our table data from our parent.
		# This method will also save column headers for us to self.fieldnames
		data_rows =  super().join_results(extracted_results)
			
		# Now convert to a texty-table (thanks tabulate).
		return "{}\n".format(tabulate(data_rows, headers = self.fieldnames, numalign = "center", stralign = "center"))
					
		
	@classmethod
	def recursive_subclasses(self):
		"""
		Recursively get all the subclasses of this class.
		"""
		return Summary_extractor.recursive_subclasses()
	
class Long_table_group_extractor(Long_tabular_group_extractor):
	"""
	Group extractor for long tables in text format (good for reading orbitals, atoms etc in text files and terminals).
	"""
	
	@classmethod
	def tabulate(self, fieldnames, table_data):
		"""
		Convert some raw table data into our desired format (texty-table in our case).
		
		:param fieldnames: List of fieldnames/table headers.
		:param table_data: List of table rows.
		"""
		return "{}\n".format(tabulate(table_data, headers = fieldnames, numalign = "center", stralign = "center", floatfmt=".4f"))

	@classmethod
	def recursive_subclasses(self):
		"""
		Recursively get all the subclasses of this class.
		"""
		return Long_table_extractor.recursive_subclasses()
	
	
