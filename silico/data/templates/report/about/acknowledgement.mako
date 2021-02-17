## -*- coding: utf-8 -*-

<%page args="title, citations" />
## Citations is a list of dictionaries of the from {'name', 'number'}

<div class="acknowledgement">
    <div class="acknowledgement__title">${title}:</div>
    <div class="acknowledgement__body">
        %for index,citation in enumerate(citations):
        <%include file="/citation/citation.mako" args="name = citation['name'], number = citation['number']" />${"," if index != len(citations) -1 else ""}
        %endfor
    </div>
</div>