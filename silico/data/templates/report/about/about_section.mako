## -*- coding: utf-8 -*-

<%!
    import silico
%>

<%
    acknowledgements = [
        {"title": "Extraction and processing of results", "citations": [{"name": "cclib", "numbers": 1}]},
        {"title": "Rendering of 3D images", "citations": [{"name": "VMD", "numbers": 2}, {"name": "Tachyon", "numbers": 3}]},
        {"title": "Rendering of graphs", "citations": [{"name": "Matplotlib", "numbers": 4}]},
        {"title": "Calculation of CIE colour coordinates", "citations": [{"name": "Colour Science", "numbers": 5}]},
        {"title": "Generation of reports", "citations": [{"name": "Mako", "numbers": 6}, {"name": "Weasyprint", "numbers": 7}]},
        {"title": "Scientific constants", "citations": [{"name": "SciPy", "numbers": 8}]},
        {"title": "Conversion of file formats", "citations": [{"name": "Pybel", "numbers": 9}, {"name": "Openbabel", "numbers": 10}]},
        {"title": "Calculation of spin-orbit coupling", "citations": [{"name": "PySOC", "numbers": 11}]},
        {"title": "Rendering of 2D structures", "citations": [{"name": "RDKit", "numbers": 12}]},
        {"title": "Saving of state during submission", "citations": [{"name": "Dill", "numbers": [13, 14]}]}
    ]
%>

<div class="section section--about">
    <h2 class="section__header">About</h2>
    <div class="section__body section__body--text section__body--about">
        <div class="title">
            <div class="title__superTitle">Silico Calculation Report</div>
            <div class="title__mainTitle">Part of the silico software package</div>
            <div class="title__subTitle">Version ${silico.version}</div>
            <div class="title__subTitle">${silico.last_updated.day} ${silico.last_updated.strftime("%B %Y")}</div>
        </div>
        <div class="acknowledgements">
            <div class="acknowledgements__title">
            Silico makes use of a number of 3<sup>rd</sup> party libraries and programs; please cite these appropriately in your works:
            </div>
            %for acknowledgement in acknowledgements:
            <%include file="acknowledgement.mako" args="title = acknowledgement['title'], citations = acknowledgement['citations']"/>
            %endfor
        </div>
    </div>
</div>