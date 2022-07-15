<%page args="report" />

<%
    error = not report.result.metadata.success or report.result.metadata.optimisation_converged == False
%>

<div class="pageHeader ${"pageHeader--bad" if error else ""}">
    <div class="pageHeader__body">
        <div class="pageHeader__info">Silico Calculation Report</div>
        <div class="pageHeader__title">  
            ${report.result.metadata.name} - ${report.result.title}
        </div>
    </div>
</div>