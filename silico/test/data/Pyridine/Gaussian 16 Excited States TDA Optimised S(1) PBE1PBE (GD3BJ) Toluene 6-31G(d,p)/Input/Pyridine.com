%Chk="Pyridine.chk"
%Rwf="Pyridine.rwf"
%NProcShared=4
%Mem=10GB
#p Opt TDA=(Singlet, nstates=10, root=1) PBE1PBE/6-31G(d,p) SCRF=(Solvent=Toluene) EmpiricalDispersion=(GD3BJ) Symmetry=(Tight) Population=(Regular) Density=(Current) 6D 10F GFInput

Pyridine_Excited_States_TDA_Optimised_S_1__PBE1PBE__GD3BJ__Toluene_6_31G_d_p_

0, 1
C           1.13865        -0.71996        -0.00002
C           1.19540         0.67092         0.00002
C           0.00004         1.38179         0.00002
C          -1.19537         0.67098         0.00002
C          -1.13870        -0.71990         0.00003
N          -0.00003        -1.41721        -0.00006
H           2.05706        -1.30465        -0.00004
H           2.15406         1.17940         0.00004
H           0.00006         2.46782         0.00005
H          -2.15399         1.17953         0.00002
H          -2.05712        -1.30454        -0.00004

