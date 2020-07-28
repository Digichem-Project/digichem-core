from silico.report import Report
from mako.lookup import TemplateLookup
from pathlib import Path
import pkg_resources
import silico
import shutil

class HTML_report(Report):
	"""
	A report type that produces html output.
	"""
	
	FULL_REPORT_NAME = "report.mako"
	ATOMS_REPORT_NAME = "atom_report.mako"
	
	def __init__(self, *args, **kwargs):
		"""
		Constructor for HTML_reports.
		"""
		super().__init__(*args, **kwargs)
		
		# Set our template dir. Eventually this will be changeable by the user so they can use their own templates.
		self.template_dir = self.default_template_directory
		self.src_static_dir = self.default_src_static_directory
		self.report_html_file = None
		self.report_type = None
	
	@property
	def default_template_directory(self):
		#return Path(pkg_resources.resource_filename('silico', 'data/templates'))
		return Path(silico.default_template_directory(), "report")
	
	@property
	def default_src_static_directory(self):
		return Path(pkg_resources.resource_filename('silico', 'data/static'))
	
	def _write(self, output, *, report_type = None, **kwargs):
		"""
		Write this HTML_report to file.
		
		:param report_template: Filename (relative to the report/report subdirectory of the template dir) of the report to make.
		"""
		# Call our parent first (this creates a lot of images and calls make() for us).
		super()._write(output, **kwargs)
		
		# Do no processing on output; it should be an exact path to a file to write.
		self.report_html_file = output
		self.static_dir = Path(self.report_html_file.parents[0], "static")
		
		# The type of report.
		self.report_type = report_type if report_type is not None else "full"
		if self.report_type == "full":
			report_template = self.FULL_REPORT_NAME
		elif self.report_type == "atoms":
			report_template = self.ATOMS_REPORT_NAME
		else:
			raise TypeError("Unknown report type '{}'".format(self.report_type))
		
		# Now get and load our template.
		template_body = TemplateLookup(directories = str(self.template_dir)).get_template("/report/{}".format(report_template)).render_unicode(result = self)
		
		try:
			# Write our template to file.
			with open(self.report_html_file, "wt") as report_html_file:
				report_html_file.write(template_body)
		except Exception:
			self.report_html_file = None
			raise
		
		# Finally, copy our static stuff.
		# Try and delete the existing dist directory if it already exists (so we update it).
		try:
			shutil.rmtree(str(self.static_dir))
		except FileNotFoundError:
			# This is ok.
			pass
		
		try:
			shutil.copytree(str(self.src_static_dir), str(self.static_dir))
		except FileExistsError:
			# This is maybe ok...
			pass
		
		# Sadly distutils is bugged if it gets called twice.
		#distutils.dir_util.copy_tree(str(self.src_static_dir), str(self.static_dir))
		
		# Done.
		
		
		