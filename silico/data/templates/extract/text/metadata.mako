## -*- coding: utf-8 -*-
##
<%!
    from silico.exception import Result_unavailable_error
    from silico.result.excited_states import Energy_state
    from silico import misc
%>\
##
<%page args="metadata, result_name = ''"/>\
##
<%
    # Note sure if this is even possible (I suspect we would have crashed well before; metadata is pretty vital), but no harm in checking.
    if metadata is None:
        raise Result_unavailable_error("metadata", "there is no metadata")
%>\
##
## We ignore 'result_name' because we always display it anyway.
<%include file="title.mako" args="title='Metadata', result_name=''"/>
##
##
Name: ${metadata.name}
%if metadata.user is not None:
Username: ${metadata.user}
%endif
%if metadata.date is not None:
Date: ${misc.date_to_string(metadata.date)}
%endif
%if metadata.duration is not None:
Duration: ${misc.timedelta_to_string(metadata.duration)}
%endif
Computational package: ${metadata.package_string}
Calculations: ${metadata.calculations_string}
Methods: ${metadata.methods_string}
Functional: ${metadata.functional}
Basis set: ${metadata.basis_set}
Multiplicity: ${Energy_state.multiplicity_number_to_string(metadata.multiplicity).capitalize() if metadata.multiplicity is not None else None}
Charge: ${metadata.charge}
Orbital spin: ${metadata.orbital_spin_type}
Success: ${metadata.success}
%if metadata.optimisation_converged is not None:
Converged: ${metadata.optimisation_converged}
%endif
%if metadata.temperature is not None:
Calculation temperature /K: ${"{:0.2f}".format(metadata.temperature)}
%endif
%if metadata.pressure is not None:
Calculation pressure /atm: ${"{:0.2f}".format(metadata.pressure)}
%endif
