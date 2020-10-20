## -*- coding: utf-8 -*-
##
<%page args="program"/>\
##
$uff
    ${program.calculation.maxcycle}    ${"1" if program.calculation.modus else "-1"}    ${program.calculation.nqeq}  ! maxcycle,modus,nqeq
    ${program.calculation.iterm}                                            ! iterm
    ${program.calculation.econv}    ${program.calculation.gconv}                     ! econv,gconv
    ${program.calculation.qtot}    ${program.calculation.dfac}                      ! qtot,dfac
    ${program.calculation.epssteep}    ${program.calculation.epssearch}    ${program.calculation.dqmax}   ! epssteep,epssearch,dqmax
    ${program.calculation.mxls}    ${program.calculation.dhls}    ${program.calculation.ahls}    ! mxls,dhls,ahls
    ${program.calculation.alpha}    ${program.calculation.beta}    ${program.calculation.gamma}    ! alpha,beta,gamma
    ${"T" if program.calculation.transform else "F"}    ${"T" if program.calculation.lnumhess else "F"}    ${"T" if program.calculation.lmd else "F"}    ! transform,lnumhess,lmd