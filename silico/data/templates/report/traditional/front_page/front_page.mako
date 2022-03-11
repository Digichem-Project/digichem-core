<%page args="report" />

<%!
    from pathlib import Path
%>
<div class="section section--frontPage">
    <div class="section__body section__body--frontPage">
        <div class="header">
            <img class="header__banner" src="static/image/banner.jpg">
            <img class="header__name" src="static/image/zysmancolman_group_title_en.png">
        </div>
        <div class="title title--report">
        	## We show a different title if the calculation failed.
        	%if not report.result.metadata.success:
        	<h1 class="title__superTitle title__superTitle--report title__superTitle--bad">Failed Calculation Report</h1>
        	%elif report.result.metadata.optimisation_converged == False:
        	<h1 class="title__superTitle title__superTitle--report title__superTitle--bad">Unconverged Calculation Report</h1>
        	%else:
            <h1 class="title__superTitle title__superTitle--report">Calculation Report</h1>
            %endif
            <h2 class="title__mainTitle title__mainTitle--report">${Path(report.result.metadata.name).name}</h2>
            %if len(report.result.metadata.calculations) > 0:
            <div class="title__subTitle title__subTitle--report">${report.result.title}</div>
            %endif
        </div>
        <div class="imageBlock imageBlock--multi imageBlock--frontPage">
            <div class="image__aligner image__aligner--frontPage">
            	%if report.options.report['front_page_image'] == 'skeletal' and 'skeletal' in report.images:
                <img class="image__img image__img--frontPage" src="${report.relative_image('skeletal')}">
                %elif report.options.report['front_page_image'] == 'rendered' and 'structure' in report.images:
                <img class="image__img image__img--frontPage" src="${report.relative_image('structure', 'x0y0z0')}">
        		%endif
            </div>
        </div>
    </div>
</div>