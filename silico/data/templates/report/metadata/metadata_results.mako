## -*- coding: utf-8 -*-

<%!
    from silico.result.excited_states import Energy_state
    from silico import misc
%>

<%page args="metadata, title = None" />


<div class="resultsContainer">
    <div class="reportHeader reportHeader--minor reportHeader--results">${"Metadata" if title is None else title}</div>
    <table class="results">
        %if metadata.date is not None:
        <tr>
            <td class="results__name">Date:</td>
            <td class="results__value">${misc.date_to_string(metadata.date)}</td>
        </tr>
        %endif
        %if metadata.duration is not None:
        <tr>
            <td class="results__name">Duration:</td>
            <td class="results__value">${misc.timedelta_to_string(metadata.duration)}</td>
        </tr>
        %endif
        <tr>
            <td class="results__name">Success:</td>
            %if metadata.success:
            <td class="results__value results__value--good">True</td>
            %elif not metadata.success:
            <td class="results__value results__value--bad">False</td>
            %else:
            <td class="results__value results__value--bad">Unknown</td>
            %endif
        </tr>
        %if metadata.optimisation_converged is not None:
        <tr>
            <td class="results__name">Converged:</td>
            %if metadata.optimisation_converged:
            <td class="results__value results__value--good">True</td>
            %else:
            <td class="results__value results__value--bad">False</td>
            %endif
        </tr>
        %endif
        %if metadata.package_string != "":
        <tr>
            <td class="results__name">Computational package:</td>
            <td class="results__value">${metadata.package_string}</td>
        </tr>
        %endif
        %if len(metadata.methods) > 0:
        <tr>
            <td class="results__name">Methods:</td>
            <td class="results__value">${", ".join(metadata.methods)}</td>
        </tr>
        %endif
        %if metadata.functional is not None:
        <tr>
            <td class="results__name">Functional:</td>
            <td class="results__value">${metadata.functional}</td>
        </tr>
        %endif
        %if metadata.basis_set is not None:
        <tr>
            <td class="results__name">Basis set:</td>
            <td class="results__value">${metadata.basis_set}</td>
        </tr>
        %endif
        %if len(metadata.calculations) > 0:
        <tr>
            <td class="results__name">Calculations:</td>
            <td class="results__value">${metadata.calculations_string}</td>
        </tr>
        %endif
        %if metadata.orbital_spin_type is not None:
        <tr>
            <td class="results__name">Orbital spin:</td>
            <td class="results__value">${metadata.orbital_spin_type}</td>
        </tr>
        %endif
        %if metadata.multiplicity is not None:
        <tr>
            <td class="results__name">Multiplicity:</td>
            <td class="results__value">${metadata.multiplicity} (${Energy_state.multiplicity_number_to_string(metadata.multiplicity)})</td>
        </tr>
        %endif
        %if metadata.temperature is not None:
        <tr>
            <td class="results__name">Calc temperature:</td>
            <td class="results__value">${metadata.temperature} K</td>
        </tr>
        %endif
        %if metadata.pressure is not None:
        <tr>
            <td class="results__name">Calc pressure:</td>
            <td class="results__value">${metadata.pressure} atm</td>
        </tr>
        %endif
        %if metadata.num_calculations != 1:
        <tr>
            <td class="results__name">No. merged calculations:</td>
            <td class="results__value">${metadata.num_calculations}</td>
        </tr>
        %endif
    </table>
</div>