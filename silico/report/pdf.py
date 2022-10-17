import weasyprint
from weasyprint import HTML
from packaging import version
from mako.lookup import TemplateLookup
import shutil
import warnings

# Silico imports.
import silico.log
from silico.report.html import HTML_report


# Check weasy version.
if version.parse(weasyprint.VERSION) < version.parse("54.2"):
    warnings.warn("Weasyprint version < 54.2 is only partially supported. Some elements of PDF reports may display incorrectly.")


class PDF_report(HTML_report):
    """
    A report type that produces pdf output.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.prog_version = silico.version
        
    def _write(self, output, **kwargs):
        """
        Write this PDF_report to file.
        
        :param output: Filename/path to a pdf file to write to.
        """
        # Call our parent first, which creates our html report and stores its path under self.report_html_file.
        super()._write(output.with_suffix(".html"), **kwargs)
        
        # Get our pdf name. Same as 'output' in this case.
        self.pdf_file = output
        # We also need an absolute path for weasyprint.
        self.absolute_pdf_file_path = self.pdf_file.resolve()

        silico.logging.get_logger().info("Writing PDF file '{}'".format(self.pdf_file))
        
        # Now render our finished pages.
        main_doc = HTML(self.report_html_file, base_url=str(self.absolute_pdf_file_path)).render()
        
        # Finally, write the file to disk.
        main_doc.write_pdf(self.pdf_file)
        
    def cleanup(self):
        """
        Remove any intermediate files that may have been created by this object.
        """
        super().cleanup()
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


class PDF_report_legacy(HTML_report):
    """
    A report type that produces pdf output.
    
    This legacy class computes page header and footer elements manually. This is no longer necessary (and is actually unsupported) from weasyprint 55.0.
    
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
    
    def __init__(self, result, *, options, calculation = None):
        """
        Constructor for PDF report objects.
        
        :param result: A parsed result set object to render a report of.
        :param options: A silico Config dictionary which contains various options that control the appearance of this report.
        :param calculation: An optional calculation class about which this report is being written. This is used by the submission mechanism to reuse certain calculation options when, for example, calculating Turbomole cube files and Gaussian NTOs etc.
        """
        super().__init__(result, options = options, calculation = calculation)
        self.prog_version = silico.version
    
    def _write(self, output, **kwargs):
        """
        Write this PDF_report to file.
        
        :param output: Filename/path to a pdf file to write to.
        """
        # Call our parent first, which creates our html report and stores its path under self.report_html_file.
        super()._write(output.with_suffix(".html"), **kwargs)
        
        # Get our pdf name. Same as 'output' in this case.
        self.pdf_file = output
        # We also need an absolute path for weasyprint.
        self.absolute_pdf_file_path = self.pdf_file.resolve()

        silico.log.get_logger().info("Writing PDF file")
        
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
        html_text = TemplateLookup(directories = self.template_directory(self.report_style), strict_undefined = True).get_template(template_file).render_unicode(report = self, *args, **kwargs)
        
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
            
            if page_body is None:
                # Couldn't find the body element...
                continue
            
            # Should we display error colors?
            error = not self.result.metadata.success or self.result.metadata.optimisation_converged == False
            
            # Get our header and footer.
            header_body = self._compute_overlay_element('/page/header.mako', error = error)
            footer_body = self._compute_overlay_element('/page/footer.mako', prog_version = self.prog_version, page_number = page_number +1, pages = len(main_doc.pages), error = error)
            
            # Disable header and footer fot atoms mini-reports (for now).
            
            if self.report_type == "atoms":
                header_body = None
                footer_body = None
                
            if header_body is not None:
                page_body.children += header_body.all_children()
            if footer_body is not None:
                page_body.children += footer_body.all_children()
                
    def cleanup(self):
        """
        Remove any intermediate files that may have been created by this object.
        """
        super().cleanup()
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
