<%page args="report, error" />

<html>
    <head>
        <link rel="stylesheet" type="text/css" href="static/css/page_header.css">
    </head>
    <div class="header ${"header--bad" if error else ""}">
    	<div class="report__info">Silico Calculation Report</div>
        <div class="report__title">  
            ${report.result.metadata.name} - ${report.result.title}
        </div>
    </div>
</html>