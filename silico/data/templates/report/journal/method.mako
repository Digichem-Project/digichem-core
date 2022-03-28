## -*- coding: utf-8 -*-

<%page args="report"/>

<div class="content">
    <h5>Analysis</h5>
	The report presented here was generated using the Silico software package.
	This toolset relies upon a number of third-party applications and libraries which should be cited appropriately in derivative works.
	In particular, the calculation results described within were parsed by the cclib library.<%include file="/citation.mako" args="citation = 'cclib', report = report"/>
	Scientific constants which were used, among other things, for the interconversion of scientific units were provided by SciPy.<%include file="/citation.mako" args="citation = 'SciPy', report = report"/>
	%if len(report.result.excited_states) > 0 or len(report.result.adiabatic_emission) > 0 or len(report.result.vertical_emission) > 0:
		Commission internationale de l'éclairage (CIE) coordinates, along with visual representations of the equivalent colour, were calculated using the Colour Science library.<%include file="/citation.mako" args="citation = 'ColourScience', report = report"/>
	%endif
	%if len(report.result.spin_orbit_coupling):
		Spin-orbit coupling (SOC, H<sub>SO</sub>) was calculated using a custom implementation of the PySOC program.<%include file="/citation.mako" args="citation = 'PySOC', report = report"/>
	%endif
	Three-dimensional plots of atom positions and calculated densities, including molecular orbitals,
	were rendered using Visual Molecular Dynamics (VMD)<%include file="/citation.mako" args="citation = 'VMD', report = report"/>
	and the Tachyon ray-tracer.<%include file="/citation.mako" args="citation = 'Tachyon', report = report"/>
	%if report.options['report']['front_page_image'] == "skeletal":
		The two-dimensional, skeletal-style structure image was rendered using
		%if report.options['skeletal_image']['render_backend'] == "rdkit":
			the RDKit library.<%include file="/citation.mako" args="citation = 'RDKit', report = report"/>
		%elif report.options['skeletal_image']['render_backend'] == "obabel":
			the OpenBabel library.<%include file="/citation.mako" args="citation = 'Openbabel', report = report"/>
		%else:
			an unrecognised backend.
		%endif
	%endif
	Finally, two-dimensional graphs were plotted using the MatPlotlib library,<%include file="/citation.mako" args="citation = 'Matplotlib', report = report"/>
	while this report itself was prepared using the Mako template library<%include file="/citation.mako" args="citation = 'Mako', report = report"/>
	and the Weasyprint library<%include file="/citation.mako" args="citation = 'Weasyprint', report = report"/>, the latter of which was responsible for generarion of the PDF file.
</div>