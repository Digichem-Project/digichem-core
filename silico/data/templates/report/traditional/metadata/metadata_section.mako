## -*- coding: utf-8 -*-

<%!
    from silico.result.excited_state import Energy_state
    import datetime
%>

<%page args="result" />

<%
    report_date = report_date if report_date is not None else datetime.today()
%>

<div class="section section--fullPage">
    <h2 class="section__header">Metadata & Summary</h2>
    <div class="section__body section__body--metadata">
        <div class="resultsContainer">
            <div class="reportHeader reportHeader--minor reportHeader--results">Metadata</div>
            <table class="results">
                <tr>
                    <td class="results__name">Report date:</td>
                    <td class="results__value">${report_date}</td>
                </tr>
                <tr>
                    <td class="results__name">Success:</td>
                    %if result.metadata.success:
                    <td class="results__value">True</td>
                    %elif not result.metadata.success:
                    <td class="results__value results__value--bad">False</td>
                    %else:
                    <td class="results__value results__value--bad">Unknown</td>
                    %endif
                </tr>
                %if result.metadata.optimisation_converged is not None:
                <tr>
                    <td class="results__name">Converged:</td>
                    %if result.metadata.optimisation_converged:
                    <td class="results__value">True</td>
                    %else:
                    <td class="results__value results__value--bad">False</td>
                    %endif
                </tr>
                %if result.metadata.package_string != "":
                <tr>
                    <td class="results__name">Computational package:</td>
                    <td class="results__value">${result.metadata.package_string}</td>
                </tr>
                %endif
                %if len(result.metadata.methods) > 0:
                <tr>
                    <td class="results__name">Methods:</td>
                    <td class="results__value">${", ".join(result.metadata.methods)}</td>
                </tr>
                %endif
                %if result.metadata.functional is not None:
                <tr>
                    <td class="results__name">Functional:</td>
                    <td class="results__value">${result.metadata.functional}</td>
                </tr>
                %endif
                %if result.metadata.basis_set is not None:
                <tr>
                    <td class="results__name">Basis set:</td>
                    <td class="results__value">${result.metadata.basis_set}</td>
                </tr>
                %endif
                %if result.metadata.multiplicity is not None:
                <tr>
                    <td class="results__name">Multiplicity:</td>
                    <td class="results__value">${result.metadata.multiplicity} (${Energy_state.multiplicity_number_to_string(result.metadata.multiplicity)})</td>
                </tr>
                %endif
            </table>
        </div>
        <div class="resultsContainer">
            <div class="reportHeader reportHeader--minor reportHeader--results">Summary</div>
            <table class="results">
                %if len(result.SCF_energies) > 0:
                <tr>
                    <td class="results__name">SCF energy:</td>
                    <td class="results__value">${"{:0.4f}".format(result.SCF_energies.final)} eV</td>
                </tr>
                %endif
                %if len(result.MP_energies) > 0:
                <tr>
                    <td class="results__name">MP energy:</td>
                    <td class="results__value">${"{:0.4f}".format(result.MP_energies.final)} eV</td>
                </tr>
                %endif
                %if len(result.CC_energies) > 0:
                <tr>
                    <td class="results__name">CC energy:</td>
                    <td class="results__value">${"{:0.4f}".format(result.CC_energies.final)} eV</td>
                </tr>
                %endif
                <%include file="/excited_states/dest_result_rows.mako" args="excited_states = result.excited_states" />
                %if len(result.vibrations) > 0:
                <tr>
                    <td class="results__name">Negative frequencies:</td>
                    %if len(result.vibrations.negative_frequencies) == 0:
                        <td class="results__value">${len(result.vibrations.negative_frequencies)}</td>
                    %else:
                    <td class="results__value results__value--bad">${len(result.vibrations.negative_frequencies)}</td>
                    %endif
                </tr>
                %endif
            </table>
            %if result.metadata.temperature is not None or result.metadata.pressure is not None:
            <div class="reportHeader reportHeader--minor reportHeader--results">Thermodynamics</div>
            <table class="results">
                %if result.metadata.temperature is not None:
                <tr>
                    <td class="results__name">Calc temperature:</td>
                    <td class="results__value">${result.metadata.temperature} K</td>
                </tr>
                %endif
                %if result.metadata.pressure is not None:
                <tr>
                    <td class="results__name">Calc pressure:</td>
                    <td class="results__value">${result.metadata.pressure} atm</td>
                </tr>
                %endif
            </table>
            %endif
        </div>
    </div>
</div>