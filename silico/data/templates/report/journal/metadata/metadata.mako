## -*- coding: utf-8 -*-

<%page args="metadata, report, primary = False"/>

<%!
	from silico import misc
	from silico.result.excited_state import Energy_state
	from silico.misc.text import text_integer, andjoin, ordinal_suffix
	import inflect
%>

<%
	inflector = inflect.engine()
%>

<div class="content">
	%if primary:
	<h5>Metadata</h5>
	%endif
	%if metadata.num_calculations == 1:
	## Just one calculation to deal with.
	The calculation of the ${metadata.human_calculations_string} was performed
	%else:
	## Multiple calculations.
	This report was generated from the combined results of ${text_integer(metadata.num_calculations)} individual calculations. The individual metadatas for each separate calculation are presented in the following sections, and the overall calculation was performed
	%endif
	using
	%if metadata.package_string != "":
	the <div class="result"><div class="result__title">${metadata.package_string}</div></div> program,
	%endif
	the <div class="result"><div class="result__title">${andjoin(metadata.methods)}</div></div> ${inflector.plural("method", len(metadata.methods))}
	%if metadata.functional is not None:
	with the <div class="result"><div class="result__title">${metadata.functional}</div></div> functional
	%endif
	%if metadata.basis_set is not None:
	and the <div class="result"><div class="result__title">${metadata.basis_set}</div></div> basis set.
	%else:
	and using a basis set that could not be determined.
	%endif
	%if metadata.date is not None and metadata.duration is not None:
	It was completed on the <div class="result"><div class="result__title">${metadata.date.strftime("%d<sup>{}</sup> %B %Y".format(ordinal_suffix(metadata.date.strftime("%d"))))}</div></div>
	after a total duration of <div class="result"><div class="result__title">${misc.timedelta_to_string(metadata.duration)}</div></div>
	%elif metadata.date is not None:
	It was completed on the <div class="result"><div class="result__title">${metadata.date.strftime("%d<sup>{}</sup> %B %Y".format(ordinal_suffix(metadata.date.strftime("%d"))))}</div></div>
	%elif metadata.duration is not None:
	It was completed after a total duration of <div class="result"><div class="result__title">${misc.timedelta_to_string(metadata.duration)}</div></div>
	%endif
	%if metadata.date is not None or metadata.duration is not None:
		%if metadata.success:
		and <div class="result"><div class="result__title result__value--good">finished successfully</div></div>.
		%elif metadata.success is False:
		but <div class="result"><div class="result__title result__value--bad">did not finish successfully</div></div>.
		%else:
		but <div class="result"><div class="result__title result__value--bad">it could not be determined whether the calculation finished correctly</div></div>.
		%endif
	%else:
		The calculation
		%if metadata.success:
		<div class="result"><div class="result__title result__value--good">finished successfully</div></div>.
		%elif metadata.success is False:
		<div class="result"><div class="result__title result__value--bad">did not finish successfully</div></div>.
		%else:
		<div class="result"><div class="result__title result__value--bad">it could not be determined whether the calculation finished successfully</div></div>.
		%endif
	%endif
	%if metadata.multiplicity is not None:
	The base multiplicity of the system under study was <div class="result"><div class="result__title">${metadata.multiplicity} (${Energy_state.multiplicity_number_to_string(metadata.multiplicity)})</div></div>.
	%endif
	## TODO: Can this be None?
	%if metadata.orbital_spin_type == "restricted":
		Finally, a <div class="result"><div class="result__title">restricted wavefunction</div></div> was used, resulting in a single set of doubly occupied orbitals.
	%elif metadata.orbital_spin_type == "unrestricted":
		Finally, an <div class="result"><div class="result__title">unrestricted wavefunction</div></div> was used, resulting in two sets of singly occupied orbitals, designated as either alpha or beta, to account separately for both spin up and spin down electrons.
	%endif
	%if metadata is report.result.metadata:
	## Main metadata.
	The full calculation metadata is tabulated in table ${report.captions("table", "metadata")}.
	%endif
</div>