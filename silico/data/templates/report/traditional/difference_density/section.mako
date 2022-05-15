## -*- coding: utf-8 -*-

## Template to display images of difference densities of excited states

<%page args="excited_states, report" />

<%
	diffden_images = []
	# Prepare a list of difference density images to display.
	for excited_state in excited_states:
		if excited_state.state_symbol + "_difference_density" in report.images:
			diffden_images.append((excited_state, excited_state.state_symbol + "_difference_density"))
	
%>

## We still use flexbox for some layout here, which breaks pagination in weasyprint. For now, we fix by splitting our list into fours (which is the number of MOs per page).

## This solution was taken from https://stackoverflow.com/questions/6614891/turning-a-list-into-nested-lists-in-python/6615011
%for diffden_group in [diffden_images[i:i+4] for i in range(0, len(diffden_images), 4)]:
<div class="section section--fullPage">
    <h2 class="section__header">${", ".join([excited_state.state_symbol for excited_state, image_name in diffden_group])} Differential Densities</h2>
    <div class="section__body section__body--orbital">
        %for excited_state, image_name in diffden_group:
        <%include file="/density/plot.mako" args="image_name = image_name, caption = excited_state.state_symbol + ' positive (hole) (' + report.images[image_name].primary_colour + ') & negative (electron) (' + report.images[image_name].secondary_colour + ') difference density (isovalue: ' + str(report.images[image_name].isovalue) + ')', report = report"/>
        %endfor
    </div>
</div>
%endfor