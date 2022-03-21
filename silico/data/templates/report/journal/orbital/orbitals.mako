## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
    from silico.misc.text import andjoin, listjoin
    import inflect
%>

<%
    inflector = inflect.engine()

    molecular_orbitals = report.result.molecular_orbitals
    beta_orbitals = report.result.beta_orbitals
    
    # Split our list of orbitals to render images for into two, one for HOMO-n, another for LUMO+n.
    pre_HOMO_orbitals = [orbital for orbital in report.orbitals_to_render if orbital.HOMO_difference < 0]
    HOMO_LUMO_orbitals = [orbital for orbital in report.orbitals_to_render if orbital.HOMO_difference == 0 or orbital.HOMO_difference == 1]
    post_LUMO_orbitals = [orbital for orbital in report.orbitals_to_render if orbital.HOMO_difference > 1]
    
    # Setup captions for all orbital images so they are numbered correctly (we need to account for the HOMO/LUMO overlap appearing in the middle).
    for orbital in pre_HOMO_orbitals:
        report.captions("figure", orbital.label)
    for orbital in HOMO_LUMO_orbitals:
        report.captions("figure", orbital.label)
    if len(report.result.beta_orbitals) > 0:
        report.captions("figure", "alpha_HOMO_LUMO")
        report.captions("figure", "beta_HOMO_LUMO")
    else:
        report.captions("figure", "HOMO_LUMO")
    for orbital in post_LUMO_orbitals:
        report.captions("figure", orbital.label)
%>

<%def name="text_for_orbital_plots()">
    Plots of the orbital density for the ${andjoin([orbital.label for orbital in report.orbitals_to_render])} are shown in ${inflector.plural("figure", report.orbitals_to_render)} ${listjoin([report.captions("figure", orbital.label) for orbital in report.orbitals_to_render])} respectively,
    while the orbital overlap between the HOMO and LUMO is shown in 
    %if len(report.result.beta_orbitals) > 0:
    figures ${report.captions("figure", "alpha_HOMO_LUMO")} and ${report.captions("figure", "beta_HOMO_LUMO")} (alpha and beta respectively).
    %else:
    figure ${report.captions("figure", "HOMO_LUMO")}.
    %endif
</%def>

<div class="content">
    <h5>Molecular Orbitals</h5>
    %if len(beta_orbitals) == 0:
    ## Alpha only.
    In total, ${len(molecular_orbitals)} doubly occupied molecular orbitals were calculated, divided into ${len(molecular_orbitals.occupied)} occupied orbitals and ${len(molecular_orbitals.virtual)} unoccupied (or virtual) orbitals.
    The calculated energies of the <div class="result"><div class="result__title">HOMO and LUMO</div> were <div class="result__value">${"{:.2f}".format(report.result.molecular_orbitals.HOMO_energy)} and ${"{:.2f}".format(report.result.molecular_orbitals.LUMO_energy)} eV</div></div> respectively, corresponding to a <div class="result"><div class="result__title">HOMO/LUMO band gap</div> of <div class="result__value">${"{:.2f}".format(report.result.molecular_orbitals.HOMO_LUMO_energy)} eV</div></div> (figure ${report.captions("figure", 'orbital_energies')}).
    ${text_for_orbital_plots()}
    <div class="resultImage resultImage--graph resultImage--orbitalEnergies">
        <img class="resultImage__image" src="${report.relative_image('orbital_energies')}">
        <div class="resultImage__caption"><div class="caption">Figure ${report.captions("figure", 'orbital_energies')}:</div> Graph of the calculated molecular orbital energies in close proximity to the HOMO/LUMO gap. Solid lines: occupied orbitals, dashed lines: virtual orbitals.</div>
    </div>
    
    
    %else:
    ## Alpha and beta.
    In total, ${len(molecular_orbitals) + len(beta_orbitals)} singly occupied molecular orbitals were calculated, divided into ${len(molecular_orbitals.occupied)} alpha occupied orbitals, ${len(beta_orbitals.occupied)} beta occupied orbitals, ${len(molecular_orbitals.virtual)} alpha unoccupied (or virtual) orbitals and ${len(beta_orbitals.virtual)} beta unoccupied orbitals.
    The calculated energies of the <div class="result"><div class="result__title">alpha and beta HOMOs</div> were <div class="result__value">${"{:.2f}".format(report.result.molecular_orbitals.HOMO_energy)} and ${"{:.2f}".format(report.result.beta_orbitals.HOMO_energy)} eV</div></div> respectively, while the energies of the <div class="result"><div class="result__title">alpha and beta LUMOs</div> were <div class="result__value">${"{:.2f}".format(report.result.molecular_orbitals.LUMO_energy)} and ${"{:.2f}".format(report.result.beta_orbitals.LUMO_energy)} eV</div></div>. These values correspond to a calculated <div class="result"><div class="result__title">HOMO/LUMO band gap</div> of <div class="result__value">${"{:.2f}".format(report.result.molecular_orbitals.HOMO_LUMO_energy)} and ${"{:.2f}".format(report.result.beta_orbitals.HOMO_LUMO_energy)} eV</div></div> for the alpha and beta case respectively (figures ${report.captions("figure", 'orbital_energies')}).
   ${text_for_orbital_plots()}
    <div class="resultImage resultImage--graph resultImage--orbitalEnergies">
        <img class="resultImage__image resultImage__image--pair" src="${report.relative_image('alpha_orbital_energies')}"><!--
        --><img class="resultImage__image resultImage__image--pair" src="${report.relative_image('beta_orbital_energies')}">
        <div class="resultImage__caption"><div class="caption">Figure ${report.captions("figure", 'orbital_energies')}:</div> Graph of the calculated molecular orbital energies in close proximity to the HOMO/LUMO gap. Solid lines: occupied orbitals, dashed lines: virtual orbitals.</div>
    </div>
    %endif
    ##
    ## Orbital plots.
    %for orbital_list in [pre_HOMO_orbitals, HOMO_LUMO_orbitals, post_LUMO_orbitals]:
        %for orbital in orbital_list:
        <%include file="/geometry/image.mako", args="image_name = orbital.label, caption = 'Orbital density plots of the {}, plotted with isovalue: {}'.format(orbital.label, report.images[orbital.label].isovalue), report = report"/>
        %endfor
        %if orbital_list is HOMO_LUMO_orbitals:
            ## Simple brute force to determine whether we have alpha/beta or unrestricted.
            %for HOMO_LUMO_name in ["alpha_HOMO_LUMO", "beta_HOMO_LUMO", "HOMO_LUMO"]:
            ## Also the HOMO/LUMO plotted on top of each other please.
                %if HOMO_LUMO_name in report.images:
                <%
                    if HOMO_LUMO_name == "alpha_HOMO_LUMO":
                        prefix = " alpha" 
                    elif HOMO_LUMO_name == "beta_HOMO_LUMO":
                        prefix = " beta"
                    else:
                        prefix = ""
                %>
                <%include file="/geometry/image.mako", args="image_name = HOMO_LUMO_name, caption = 'Orbital density plots of the{} HOMO ({}) and LUMO ({}), plotted simultaneously with isovalue: {}'.format(prefix, report.images[HOMO_LUMO_name].primary_colour, report.images[HOMO_LUMO_name].secondary_colour, report.images[HOMO_LUMO_name].isovalue), report = report"/>
                %endif
            %endfor
        %endif
    %endfor
</div>