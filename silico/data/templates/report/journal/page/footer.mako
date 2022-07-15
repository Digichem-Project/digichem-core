<%page args="report"/>

<%
    error = not report.result.metadata.success or report.result.metadata.optimisation_converged == False
    page_number = 1
    pages = 10
%>

<div class="pageFooter ${"pageFooter--bad" if error else ""}">
    <div class="pageFooter__body">
        <div class="pageFooter__name">Silico ${report.prog_version}</div>
        <div class="pageFooter__calcs"></div>
        %if pages > 1:
        <div class="pageFooter__page">Page <span class="pageFooter__pageCounter"></span>
        of <span class="pageFooter__pageCount"></span></div>
        %endif
    </div>
</div>