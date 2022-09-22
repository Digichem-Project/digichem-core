## -*- coding: utf-8 -*-

## Template to display images of generic orbitals (ie, not the FMOs)

<%page args="orbitals, report" />

## We still use flexbox for some layout here, which breaks pagination in weasyprint. For now, we fix by splitting our list into fours (which is the number of MOs per page).

## This solution was taken from https://stackoverflow.com/questions/6614891/turning-a-list-into-nested-lists-in-python/6615011
%for mo_group in [orbitals[i:i+4] for i in range(0, len(orbitals), 4)]:
<div class="section section--fullPage">
    <h2 class="section__header">${", ".join([orbital.label for orbital in mo_group])}</h2>
    <div class="section__body section__body--orbital">
        %for orbital in mo_group:
        <%include file="/density/plot.mako" args="image_name = orbital.label, caption = orbital.label + ' density (isovalue: ' + str(report.images[orbital.label].isovalue) + ')', report = report"/>
        %endfor
    </div>
</div>
%endfor