<%page args="metadata" />

<%!
	from pathlib import Path
%>

<div class="header">
	<img class="header__banner" src="static/image/banner.jpg">
	<img class="header__name" src="static/image/zysmancolman_group_title_en.png">
</div>
<div class="title">
	<h1 class="title__superTitle">Calculation Report</h1>
	<h2 class="title__name">${Path(metadata.name).name}</h2>
	%if len(result.metadata.calculations) > 0:
	<div class="title__types">${", ".join(result.metadata.calculations)}</div>
	%endif
</div>