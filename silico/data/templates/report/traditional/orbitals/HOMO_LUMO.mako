## -*- coding: utf-8 -*-

<%!
    from silico.exception import Result_unavailable_error 
%>

<%page args="molecular_orbitals, report"/>

<%
    # First work out the type of orbitals we are dealing with (alpha, beta, or restricted).
    orbital_type = molecular_orbitals.spin_type
    
    spin_prefix = orbital_type + "_" if orbital_type != "none" else ""
    
    # Render our alternative orbital energy diagrams (we don't actually use these in the template (yet?), just write to file.
    report.relative_image(spin_prefix + 'HOMO_LUMO_energies')
    
    # Our title.
    title = "HOMO & LUMO"
    
    # Add on some more to our title if we have unrestricted orbitals.
    if orbital_type != "none":
        title += " ({})".format(orbital_type.capitalize())
        
    # Get our FMOs.
    try:
        HOMO = molecular_orbitals.get_orbital(HOMO_difference = 0)
    except Result_unavailable_error:
        HOMO = None
        
    try:
        LUMO = molecular_orbitals.get_orbital(HOMO_difference = 1)
    except Result_unavailable_error:
        LUMO = None
    
    combined_orbital_image_name = spin_prefix + "HOMO_LUMO"
%>

<div class="section">
    <h2 class="section__header">${title}</h2>
    <div class="section__body section__body--orbital">
        %if HOMO is not None:
        <div class="imageBlock imageBlock--multi imageBlock--orbital">
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(HOMO.label, 'x0y0z0')}">
                </div>
                <div class="image__caption">X/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(HOMO.label, 'x90y0z0')}">
                </div>
                <div class="image__caption">X/Z plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(HOMO.label, 'x0y90z0')}">
                </div>
                <div class="image__caption">Z/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(HOMO.label, 'x45y45z45')}">
                </div>
                <div class="image__caption">45° to axes</div>
            </div>
            <div class="imageBlock__caption">HOMO density (isovalue: ${report.images[HOMO.label].isovalue})</div>
        </div>
        %endif
        %if LUMO is not None:
        <div class="imageBlock imageBlock--multi imageBlock--orbital">
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(LUMO.label, 'x0y0z0')}">
                </div>
                <div class="image__caption">X/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(LUMO.label, 'x90y0z0')}">
                </div>
                <div class="image__caption">X/Z plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(LUMO.label, 'x0y90z0')}">
                </div>
                <div class="image__caption">Z/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(LUMO.label, 'x45y45z45')}">
                </div>
                <div class="image__caption">45° to axes</div>
            </div>
            <div class="imageBlock__caption">LUMO density (isovalue: ${report.images[LUMO.label].isovalue})</div>
        </div>
        %endif
        %if combined_orbital_image_name in report.images:
        <div class="imageBlock imageBlock--multi imageBlock--orbital">
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(combined_orbital_image_name, 'x0y0z0')}">
                </div>
                <div class="image__caption">X/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(combined_orbital_image_name, 'x90y0z0')}">
                </div>
                <div class="image__caption">X/Z plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(combined_orbital_image_name, 'x0y90z0')}">
                </div>
                <div class="image__caption">Z/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(combined_orbital_image_name, 'x45y45z45')}">
                </div>
                <div class="image__caption">45° to axes</div>
            </div>
            <div class="imageBlock__caption">HOMO (${report.images[combined_orbital_image_name].primary_colour}) & LUMO (${report.images[combined_orbital_image_name].secondary_colour}) density<br>(isovalue: ${report.images[combined_orbital_image_name].isovalue})</div>
        </div>
        %endif
        <div class="imageBlock imageBlock--multi imageBlock--orbital">
            <div class="image__aligner image__aligner--orbital_diagram">
                <img class="image_img image__img--orbital_diagram" src="${report.relative_image(spin_prefix + 'orbital_energies')}">
            </div>
        </div>
    </div>
</div>