## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
## An alternative define driver for restarted calculations.
##
## First (re)-set job name.
${calculation.safe_name(calculation.descriptive_name)}
##
## We will now be given a number of questions that we'll skip.











##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
## ! SCF                     !
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
<%include file="scf.mako" args="calculation = calculation"/>\
##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
## ! DFT                     !
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
<%include file="dft.mako" args="calculation = calculation"/>\
##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
## ! CC                      !
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
<%include file="cc.mako" args="calculation = calculation"/>\
##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
## ! prop                    !
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
<%include file="prop.mako" args="calculation = calculation"/>\
##
##
## All done.
*
##