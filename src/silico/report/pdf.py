from weasyprint import HTML
from mako.lookup import TemplateLookup
from silico.report.html import HTML_report
import shutil

class PDF_report(HTML_report):
	"""
	A report type that produces pdf output.
	
	A lot of the code in this file comes from the weasyprint tutorial at: https://weasyprint.readthedocs.io/en/stable/tips-tricks.html
	Here are the notes from that snippet:
	Notes:
	------
	- When Weasyprint renders an html into a PDF, it goes though several intermediate steps.
	  Here, in this class, we deal mostly with a box representation: 1 `Document` have 1 `Page`
	  or more, each `Page` 1 `Box` or more. Each box can contain other box. Hence the recursive
	  method `get_element` for example.
	  For more, see:
	  https://weasyprint.readthedocs.io/en/stable/hacking.html#dive-into-the-source
	  https://weasyprint.readthedocs.io/en/stable/hacking.html#formatting-structure
	- Warning: the logic of this class relies heavily on the internal Weasyprint API. This
	  snippet was written at the time of the release 47, it might break in the future.
	- This generator draws its inspiration and, also a bit of its implementation, from this
	  discussion in the library github issues: https://github.com/Kozea/WeasyPrint/issues/92
	"""
	
	def __init__(self, *args, prog_version, **kwargs):
		"""
		Constructor for PDF report objects.
		
		:param prog_version: Version string of out program (this gets inserted into the pdf).
		"""
		super().__init__(*args, **kwargs)
		self.prog_version = prog_version
	
	def _write(self, output, **kwargs):
		"""
		Write this PDF_report to file.
		
		:param output: Filename/path to a pdf file to write.
		"""
		# Call our parent first, which creates our html report and stores its path under self.report_html_file.
		super()._write(output.with_suffix(".html"), **kwargs)
		
		# Get our pdf name. Same as 'output' in this case.
		self.pdf_file = output
		# We also need an absolute path for weasyprint.
		self.absolute_pdf_file_path = self.pdf_file.resolve()

		# Now render our finished pages.
		main_doc = HTML(self.report_html_file, base_url=str(self.absolute_pdf_file_path)).render()

		# Apply our header and footer.
		self._apply_overlay_on_main(main_doc)
		
		# Finally, write the file to disk.
		main_doc.write_pdf(self.pdf_file)

		
	def _compute_overlay_element(self, template_file, *args, **kwargs):
		"""
		
		:param element: The html text representation of the element.
		:return: A Weasyprint pre-rendered representation of an html element.
		"""
		# We recompute the header and footer for each page so we can update page numbers etc.
		html_text = TemplateLookup(directories = self.template_dir, strict_undefined = True).get_template(template_file).render_unicode(result = self, *args, **kwargs)
		
		html = HTML(
			string=html_text,
			base_url=str(self.absolute_pdf_file_path)
		)
		element_doc = html.render()
		element_page = element_doc.pages[0]
		element_body = self.get_element(element_page._page_box.all_children(), 'body')
		element_body = element_body.copy_with_children(element_body.all_children())

		return element_body

	def _apply_overlay_on_main(self, main_doc):
		"""
		Insert the header and the footer in the main document.

		Parameters
		----------
		main_doc: Document
			The top level representation for a PDF page in Weasyprint.
		header_body: BlockBox
			A representation for an html element in Weasyprint.
		footer_body: BlockBox
			A representation for an html element in Weasyprint.
		"""
		for page_number, page in enumerate(main_doc.pages):
			page_body = self.get_element(page._page_box.all_children(), 'body')
			
			# Get our header and footer.
			header_body = self._compute_overlay_element('/page/page_header.mako')
			footer_body = self._compute_overlay_element('/page/page_footer.mako', prog_version = self.prog_version, page_number = page_number +1, pages = len(main_doc.pages))

			if header_body is not None:
				page_body.children += header_body.all_children()
			if footer_body is not None:
				page_body.children += footer_body.all_children()
				
	def cleanup_intermediate_files(self):
		"""
		Remove any intermediate files that may have been created by this object.
		"""
		super().cleanup_intermediate_files()
		# Delete our HTML file.
		try:
			self.report_html_file.unlink()
		except FileNotFoundError:
			# This is ok.
			pass
		
		# And delete our static directory.
		try:
			shutil.rmtree(str(self.static_dir))
		except FileNotFoundError:
			# This is ok?
			pass

	@classmethod
	def get_element(self, boxes, element):
		"""
		Given a set of boxes representing the elements of a PDF page in a DOM-like way, find the
		box which is named `element`.

		Look at the notes of the class for more details on Weasyprint insides.
		"""
		for box in boxes:
			if box.element_tag == element:
				return box
			return self.get_element(box.all_children(), element)
