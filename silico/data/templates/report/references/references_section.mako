## -*- coding: utf-8 -*-

<%!
    import silico.reference
%>

<div class="section section--references">
    <h2 class="section__header">Bibliography</h2>
    <div class="section__body section__body--text section__body--references">
        %for num,reference in enumerate(silico.reference.silico_references.items()):
        <%include file="RSC_reference.mako" args="reference = reference[1], reference_title = '[{}]'.format(num+1)"/>
        %endfor
    </div>
</div>