## -*- coding: utf-8 -*-
##
<%page args="atoms"/>\
##
<%
    # Build our formula string here becase we don't want any spaces between the letters, which can be tricky to do in the template.
    formula_string = ""
    
    # Go through each element and add to the string.
    for element in atoms.element_dict:
        # Add the element symbol (which is our key).
        formula_string += element
        # Add the number, unless it is 1.
        if atoms.element_dict[element] > 1:
            formula_string += "<sub>{}</sub>".format(atoms.element_dict[element])
            
    # Finally, add on our charge if we have one, but don't include the number if we are +/- 1.
    if atoms.charge == 1:
        formula_string += "<sup>+</sup>"
    elif atoms.charge == -1:
        formula_string += "<sup>-</sup>"
    elif atoms.charge != 0:
        formula_string += "<sup>{:+}</sup>".format(atoms.charge)
%>\
${formula_string}