<%page args="report, error" />

<%!
    from pathlib import Path
%>

<html>
    <head>
## There appears to be a bizarre bug in some part of weasyprint (for CentOS 7.7) that makes @font-face broken for header and footer.
##        <link rel="stylesheet" type="text/css" href="static/css/font.css">
        <link rel="stylesheet" type="text/css" href="static/css/page_header.css">
    </head>
    <header class="header ${"header--bad" if error else ""}">
        <div class=header__title>
            ##Calculation Report
            ##${Path(result.metadata.name).name} - ${", ".join(result.metadata.calculations)}
            ##${Path(report.result.metadata.name).name} - ${report.result.title}
            ${report.result.metadata.name} - ${report.result.title}
        </div>
    </header>
</html>