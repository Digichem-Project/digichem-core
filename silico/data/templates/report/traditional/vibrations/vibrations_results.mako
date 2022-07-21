## -*- coding: utf-8 -*-

<%page args="vibrations"/>

<%
    # The maximum number of negative frequencies to print here.
    max_neg_freq = 4
%>

<div class="resultsContainer">
    <div class="reportHeader reportHeader--minor reportHeader--results">Vibrational Frequencies</div>
    <table class="results">
        <tr>
            <td class="results__name">Negative frequencies:</td>
            %if len(vibrations.negative) == 0:
                <td class="results__value results__value--good">${len(vibrations.negative)}</td>
            %else:
            <td class="results__value results__value--bad">${len(vibrations.negative)}</td>
            %endif
        </tr>
        ## Only show the first 5 negative frequencies here...
        %for vibration in vibrations.negative[:max_neg_freq]:
        <tr>
            <td class="results__name">Frequency:</td>
            <td class="results__value results__value--bad">${"{:0.2f}".format(vibration.frequency)} cm<sup>-1</sup></td>
        </tr>
        %endfor
        ## If we have more than 5, print ellipsis.
        %if len(vibrations.negative) > max_neg_freq:
        <tr>
            <td class="results__name"></td>
            <td class="results__value results__value--bad">More...</td>
        </tr>
        %endif
    </table>
</div>