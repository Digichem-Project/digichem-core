## -*- coding: utf-8 -*-

## Template for displaying references in RSC style.

<%page args="reference, reference_title = None"/>

<div class="reference">
	%if reference_title is not None:
	<div class="reference__title">${reference_title}</div>
	%endif
	
	%if reference.url is not None:
	<a class="reference__link" href="${reference.url}">
	%endif
	
		<div class="reference__body">
		<%page args="reference, reference_title = None"/>
		%if reference.record['ENTRYTYPE'] == 'article':
		<%include file="RSC_journal.mako" args="reference = reference, reference_title = reference_title" />
		%else:
		<%include file="RSC_general.mako" args="reference = reference, reference_title = reference_title" />
		%endif
		</div>
		
	%if reference.url is not None:
	</a>
	%endif
</div>