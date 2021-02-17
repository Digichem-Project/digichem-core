## -*- coding: utf-8 -*-

## This is a mini version of the report that only contains atom positions (useful for including in the ESI).

<%page args="report"/>

<!DOCTYPE html>

<html>
    <head>
        <meta charset="utf-8"/>
        <% 
        stylesheets = [
            'font.css',
            'report.css',
            'table.css',
            'image.css',
            'results.css',
            'front_page.css',
            'metadata.css',
            'geometry.css',
            'mo.css',
            'vibrations.css',
            'excited_states.css',
            'excited_states_table.css',
            'reference.css',
            'summary.css',
            'color_box.css',
            'energies.css',
            'absorptions.css',
            'about.css'
        ]
        %>
        %for stylesheet in stylesheets:
        <link rel="stylesheet" type="text/css" href="static/css/${stylesheet}">
        %endfor
        <style>
            @page {
                margin: 1cm;
                size: A4;
            }
        </style>
    </head>
    <body>
        ##<%include file="/front_page/front_page.mako" args="result = result"/>
        %if len(report.result.atoms) > 0:
            <%include file="/geometry/atom_list_section.mako" args="atoms = report.result.alignment, title = False"/>
        %endif
        ##<%include file="/about/about_section.mako"/>
        ##<%include file="/references/references_section.mako"/>
    </body>
</html>