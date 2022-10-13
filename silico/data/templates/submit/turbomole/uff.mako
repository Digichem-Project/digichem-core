## -*- coding: utf-8 -*-
##
<%page args="program"/>\
##
$uff
    ${program.calculation.max_cyle}    ${"1" if program.calculation.method['mm'].modus else "-1"}    ${program.calculation.method['mm'].nqeq}  ! maxcycle,modus,nqeq
    ${program.calculation.method['mm'].iterm}                                            ! iterm
    ${program.calculation.method['mm'].econv}    ${program.calculation.method['mm'].gconv}                     ! econv,gconv
    ${program.calculation.qtot}    ${program.calculation.method['mm'].dfac}                      ! qtot,dfac
    ${program.calculation.method['mm'].epssteep}    ${program.calculation.method['mm'].epssearch}    ${program.calculation.method['mm'].dqmax}   ! epssteep,epssearch,dqmax
    ${program.calculation.method['mm'].mxls}    ${program.calculation.method['mm'].dhls}    ${program.calculation.method['mm'].ahls}    ! mxls,dhls,ahls
    ${program.calculation.method['mm'].alpha}    ${program.calculation.method['mm'].beta}    ${program.calculation.method['mm'].gamma}    ! alpha,beta,gamma
    ${"T" if program.calculation.method['mm'].transform else "F"}    ${"T" if program.calculation.method['mm'].lnumhess else "F"}    ${"T" if program.calculation.method['mm'].lmd else "F"}    ! transform,lnumhess,lmd