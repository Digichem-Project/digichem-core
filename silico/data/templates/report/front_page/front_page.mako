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
			<h1 class="title__superTitle title__superTitle--report">Calculation Report</h1>
			<h2 class="title__mainTitle title__mainTitle--report">${Path(report.result.metadata.name).name}</h2>
			%if len(report.result.metadata.calculations) > 0:
			<div class="title__subTitle title__subTitle--report">${report.result.title}</div>
			%endif
		</div>
		<div class="imageBlock imageBlock--multi imageBlock--frontPage">
			<div class="image__aligner image__aligner--frontPage">
				<img class="image__img image__img--frontPage" src="${report.relative_image('aligned_structure', 'x0y0z0')}">
			</div>
		</div>
	</div>
</div>