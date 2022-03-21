## -*- coding: utf-8 -*-

<%page args="report"/>

<%
    result = report.result
    # Decide which spin value (up/alpha or down/beta) corresponds to the 'hole' and 'electron'.
    # We do this by comparing the alpha and beta HOMO values, the higher one is the electron.
    try:
	    alpha_energy = result.molecular_orbitals.HOMO_energy
	    beta_energy = result.beta_orbitals.HOMO_energy
	    
	    if alpha_energy > beta_energy:
	        # Alpha/up/positive is electron
	        alpha_designation = "electron"
	        beta_designation = "hole"
	    else:
	        # Beta/down/negative is electron.
	        alpha_designation = "hole"
	        beta_designation = "electron"
    except Result_unavailable_error:
    	# Couldn't get HOMO and LUMO, probably missing beta orbitals or similar.
    	alpha_designation = ""
    	beta_designation = ""
%>

<div class="content">
	<h5>Spin Density</h5>
	The calculated difference in spin density between the alpha and beta cases is shown in figure ${report.captions("figure", 'positive_spin_density')} for the positive difference and figure ${report.captions("figure", 'negative_spin_density')} for the negative difference. A combined plot of both the positive and negative difference is shown in figure ${report.captions("figure", 'spin_density')}.
	<%include file="/geometry/image.mako" args="image_name = 'positive_spin_density', caption = 'Plot of the positive difference in spin density (alpha, {}), plotted with an isovalue of {}'.format(alpha_designation, report.images['positive_spin_density'].isovalue), report = report" />
	<%include file="/geometry/image.mako" args="image_name = 'negative_spin_density', caption = 'Plot of the negative difference in spin density (beta, {}), plotted with an isovalue of {}'.format(beta_designation, report.images['negative_spin_density'].isovalue), report = report" />
	<%include file="/geometry/image.mako" args="image_name = 'spin_density', caption = 'Plot of the positive (alpha, {}, {}) and negative (beta, {}, {}) difference in spin density, plotted with an isovalue of {}'.format(alpha_designation, report.images['spin_density'].primary_colour, beta_designation, report.images['spin_density'].secondary_colour, report.images['spin_density'].isovalue), report = report" />
</div>