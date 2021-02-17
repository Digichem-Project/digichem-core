## -*- coding: utf-8 -*-

## Template to display spin density.

<%page args="result, report" />

<%
    # Decide which spin value (up/alpha or down/beta) corresponds to the 'hole' and 'electron'.
    # We do this by comparing the alpha and beta HOMO values, the higher one is the electron.
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
%>


<div class="section">
    <h2 class="section__header">Spin Density</h2>
    <div class="section__body section__body--orbital">
        <div class="imageBlock imageBlock--multi imageBlock--orbital">
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('positive_spin_density', 'x0y0z0')}">
                </div>
                <div class="image__caption">X/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('positive_spin_density', 'x90y0z0')}">
                </div>
                <div class="image__caption">X/Z plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('positive_spin_density', 'x0y90z0')}">
                </div>
                <div class="image__caption">Z/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('positive_spin_density', 'x45y45z45')}">
                </div>
                <div class="image__caption">45° to axes</div>
            </div>
            <div class="imageBlock__caption">Positive spin density (${alpha_designation}) (isovalue: ${report.images['positive_spin_density'].isovalue})</div>
        </div>
        <div class="imageBlock imageBlock--multi imageBlock--orbital">
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('negative_spin_density', 'x0y0z0')}">
                </div>
                <div class="image__caption">X/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('negative_spin_density', 'x90y0z0')}">
                </div>
                <div class="image__caption">X/Z plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('negative_spin_density', 'x0y90z0')}">
                </div>
                <div class="image__caption">Z/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('negative_spin_density', 'x45y45z45')}">
                </div>
                <div class="image__caption">45° to axes</div>
            </div>
            <div class="imageBlock__caption">Negative spin density (${beta_designation}) (isovalue: ${report.images['negative_spin_density'].isovalue})</div>
        </div>
        <div class="imageBlock imageBlock--multi imageBlock--orbital">
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('spin_density', 'x0y0z0')}">
                </div>
                <div class="image__caption">X/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('spin_density', 'x90y0z0')}">
                </div>
                <div class="image__caption">X/Z plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('spin_density', 'x0y90z0')}">
                </div>
                <div class="image__caption">Z/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image('spin_density', 'x45y45z45')}">
                </div>
                <div class="image__caption">45° to axes</div>
            </div>
            <div class="imageBlock__caption">Positive (${alpha_designation}) (${report.images['spin_density'].primary_colour}) & negative (${beta_designation}) (${report.images['spin_density'].secondary_colour}) spin density (isovalue: ${report.images['spin_density'].isovalue})</div>
        </div>
    </div>
</div>