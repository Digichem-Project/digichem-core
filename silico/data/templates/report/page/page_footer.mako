<%page args="result, prog_version, page_number, pages"/>

<html>
	<head>
## There appears to be a bizarre bug in some part of weasyprint (for CentOS 7.7) that makes @font-face broken for header and footer.
##		<link rel="stylesheet" type="text/css" href="static/css/font.css">
		<link rel="stylesheet" type="text/css" href="static/css/page_footer.css">
	</head>
	<footer class="footer">
		##<div class="footer__name">${Path(result.metadata.name).name} - ${", ".join(result.metadata.calculations)}</div>
		<div class="footer__name">Silico ${prog_version}</div>
		<div class="footer__calcs"></div>
		%if pages > 1:
		<div class="footer__page">Page ${page_number} of ${pages}</div>
		%endif
	</footer>
</html>