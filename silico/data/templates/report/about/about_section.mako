## -*- coding: utf-8 -*-

<%!
    import silico
%>

<%
    acknowledgements = [
        {"title": "Extraction and processing of results", "citations": [{"name": "cclib", "number": 1}]},
        {"title": "Generation of 3D images", "citations": [{"name": "VMD", "number": 2}, {"name": "Tachyon", "number": 3}]},
        {"title": "Generation of graphs", "citations": [{"name": "Matplotlib", "number": 4}]},
        {"title": "Calculation of CIE colour coordinates", "citations": [{"name": "Colour Science", "number": 5}]},
        {"title": "Generation of report", "citations": [{"name": "Mako", "number": 6}, {"name": "Weasyprint", "number": 7}]},
        {"title": "Scientific constants", "citations": [{"name": "SciPy", "number": 8}]},
        {"title": "Conversion of file formats", "citations": [{"name": "Pybel", "number": 9}, {"name": "Openbabel", "number": 10}]},
        {"title": "Calculation of spin-orbit coupling", "citations": [{"name": "PySOC", "number": 11}]},
    ]
%>

<div class="section section--about">
    <h2 class="section__header">About</h2>
    <div class="section__body section__body--text section__body--about">
        <div class="title">
            <div class="title__superTitle"><u>cr</u>eport: <u>C</u>alculation <u>R</u>eport Generator</div>
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