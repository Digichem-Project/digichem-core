### -*- coding: utf-8 -*-
##
<%page args="reference_name, report"/>
##
<%!
	import silico.reference
%>
##
<%
	reference = silico.reference.silico_references[reference_name]
%>
##
##
<tr class="reference">
	<td class="reference__number">${report.captions("citation", reference_name)}.</td>
	<td class="reference__body">
		%if reference.url is not None:
	    	<a class="reference__link" href="${reference.url}">
	    %endif    
	    %if reference.record['ENTRYTYPE'] == 'article':
	        <%include file="journal.mako" args="reference = reference" />
	    %else:
	       	<%include file="general.mako" args="reference = reference" />
        %endif
	    %if reference.url is not None:
	    	</a>
	    %endif
	</td>
</tr>