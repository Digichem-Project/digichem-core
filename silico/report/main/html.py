# General imports.
from mako.lookup import TemplateLookup
from pathlib import Path
import pkg_resources
import silico
import shutil

# Silico imports.
from silico.report.main import Report
from silico.misc.directory import copytree
from silico.reference.cross_reference import Captions

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
        
        self.report_html_file = None
        self.report_type = None
        self.report_style = None
        # An object used in templates for recording captions and references.
        self.captions = Captions()
    
    def template_directory(self, report_style):
        return Path(silico.default_template_directory(), "report", report_style)
    
    def src_static_directory(self, report_style):
        return Path(pkg_resources.resource_filename('silico', str(Path('data/static', report_style))))
    
    def _write(self, output, *, report_type = "full", report_style = "journal", **kwargs):
        """
        Write this HTML_report to file.
        
        :param output: Path to a html to file write to.
        """
        # Call our parent first (this creates a lot of images and calls make() for us).
        super()._write(output.parent, **kwargs)
        
        # Do no processing on output; it should be an exact path to a file to write.
        self.report_html_file = output
        self.static_dir = Path(self.report_html_file.parents[0], "static")
        
        # The type of report.
        self.report_type = report_type
        self.report_style = report_style
        if self.report_type == "full":
            report_template = self.FULL_REPORT_NAME
        elif self.report_type == "atoms":
            report_template = self.ATOMS_REPORT_NAME
        else:
            raise TypeError("Unknown report type '{}'".format(self.report_type))
        
        # Now get and load our template.
        template_body = TemplateLookup(directories = str(self.template_directory(self.report_style))).get_template("/report/{}".format(report_template)).render_unicode(report = self)
        
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
        
        # Because the in-built python copytree functions appears to be broken, we'll just use cp for now.
        # This is a work around for #21.
#         subprocess.run([
#                 "cp", "-R",
#                 str(self.src_static_dir),
#                 str(self.static_dir)
#             ],
#             universal_newlines = True,
#             check = True
#         )
        # Be aware of #21
        # We only want to copy data, nothing else.
        copytree(str(self.src_static_directory(self.report_style)), str(self.static_dir), copy_function = shutil.copyfile)
        
        # Done.
        
        
        