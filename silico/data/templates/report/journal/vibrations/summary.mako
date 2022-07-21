## -*- coding: utf-8 -*-

<%page args="vibrations, report"/>

<%!
	from silico.misc.text import andjoin
%>

<%
    # The maximum number of negative frequencies to print here.
    max_neg_freq = 5
    peaks = report.images['simulated_IR_graph'].selected_peaks(0, 5)
%>

<div class="resultsTable resultsTable--summary">
    <div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "vibrations summary")}:</div>
		Summary of the properties of the calculated vibration frequencies.
	</div>
    <table class="resultsTable__table">
    	<tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">No. frequencies</td>
            <td class="resultsTable__value resultsTable__cell">${len(vibrations)}</td>
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Simulated peaks</td>
            <td class="resultsTable__value resultsTable__cell">
            	<%
            		selected_peaks = report.images['simulated_IR_graph'].selected_peaks(0, 5)
            	%>
            	${andjoin(selected_peaks)}
            	%if len(report.images['simulated_IR_graph'].peaks) > len(selected_peaks):
            	...
            	%endif
            	cm<sup>-1</sup>
            </td>
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">No. negative frequencies</td>
            %if len(vibrations.negative) == 0:
                <td class="resultsTable__value resultsTable__cell resultsTable__value--good">${len(vibrations.negative)}</td>
            %else:
            <td class="resultsTable__value resultsTable__cell resultsTable__value--bad">${len(vibrations.negative)}</td>
            %endif
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Negative frequencies</td>
            <td class="resultsTable__value resultsTable__cell">
            	%if len(vibrations.negative) > 0:
	            	${andjoin(["{:0.2f}".format(vibration.frequency) for vibration in vibrations.negative[:max_neg_freq]])}
	            	%if len(vibrations.negative) > max_neg_freq:
	            	...
	            	%endif
	            	cm<sup>-1</sup>
	            %else:
	            	N/A
	            %endif
            </td>
        </tr>
    </table>
</div>