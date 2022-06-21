<%page args="prog_version, page_number, pages, error"/>

<html>
    <head>
## There appears to be a bizarre bug in some part of weasyprint (for CentOS 7.7) that makes @font-face broken for header and footer.
##        <link rel="stylesheet" type="text/css" href="static/css/font.css">
        <link rel="stylesheet" type="text/css" href="static/css/page_footer.css">
    </head>
    <footer class="footer ${"footer--bad" if error else ""}">
        <div class="footer__name">Silico ${prog_version}</div>
        <div class="footer__calcs"></div>
        %if pages > 1:
        <div class="footer__page">Page ${page_number} of ${pages}</div>
        %endif
    </footer>
</html>