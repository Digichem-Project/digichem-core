## -*- coding: utf-8 -*-

## Template to display spin density.

<%page args="result, report" />

<%!
	from silico.exception import Result_unavailable_error
%>

<%
    # Decide which spin value (up/alpha or down/beta) corresponds to the 'hole' and 'electron'.
    # We do this by comparing the alpha and beta HOMO values, the higher one is the electron.
    try:
	    alpha_energy = result.molecular_orbitals.HOMO_energy
	    beta_energy = result.beta_orbitals.HOMO_energy
	    
	    if alpha_energy > beta_energy:
	        # Alpha/up/positive is electron
	        alpha_designation = "(electron)"
	        beta_designation = "(hole)"
	    else:
	        # Beta/down/negative is electron.
	        alpha_designation = "(hole)"
	        beta_designation = "(electron)"
    except Result_unavailable_error:
    	# Couldn't get HOMO and LUMO, probably missing beta orbitals or similar.
    	alpha_designation = ""
    	beta_designation = ""
%>


<div class="section">
    <h2 class="section__header">Spin Density</h2>
    <div class="section__body section__body--orbital">
    	<%include file="/density/plot.mako" args="image_name = 'positive_spin_density', caption = 'Positive spin density ' + alpha_designation + ' (isovalue: '+ str(report.images['positive_spin_density'].isovalue) + ')', report = report"/>
    	<%include file="/density/plot.mako" args="image_name = 'negative_spin_density', caption = 'Negative spin density ' + beta_designation + ' (isovalue: '+ str(report.images['negative_spin_density'].isovalue) + ')', report = report"/>
    	<%include file="/density/plot.mako" args="image_name = 'spin_density', caption = 'Positive ' + alpha_designation + ' (' + report.images['spin_density'].primary_colour + ') & negative ' + beta_designation + ' (' + report.images['spin_density'].secondary_colour + ') spin density (isovalue: ' + str(report.images['spin_density'].isovalue) + ')', report = report"/>
    </div>
</div>