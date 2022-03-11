## -*- coding: utf-8 -*-

<%page args="report"/>

<!DOCTYPE html>

<html>
    <head>
        <meta charset="utf-8"/>
        <% 
        stylesheets = [
        
        ]
        %>
        %for stylesheet in stylesheets:
        <link rel="stylesheet" type="text/css" href="static/css/journal/${stylesheet}">
        %endfor
    </head>
    <body>
    	<%include file="/header.mako" args="report = report"/>
    </body>
</html>