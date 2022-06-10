%Chk="Naphthalene.chk"
%Rwf="Naphthalene.rwf"
%NProcShared=4
%Mem=10GB
#p Opt TDA=(Singlet, nstates=10, root=1) PBE1PBE/6-31G(d,p) SCRF=(Solvent=Toluene) EmpiricalDispersion=(GD3BJ) Symmetry=(Tight) Population=(Regular) Density=(Current) 6D 10F GFInput

Naphthalene_Excited_States_TDA_Optimised_S_1__PBE1PBE__GD3BJ__Toluene_6_31G_d_p_

0, 1
C           0.00000         1.24046         1.39914
C           0.00000         2.42600         0.70664
C           0.00000         2.42600        -0.70664
C           0.00000         1.24046        -1.39914
C          -0.00000         0.00000        -0.71423
C          -0.00000         0.00000         0.71423
C          -0.00000        -1.24046         1.39914
C          -0.00000        -1.24046        -1.39914
C          -0.00000        -2.42600        -0.70664
C          -0.00000        -2.42600         0.70664
H           0.00000         1.23670         2.48620
H           0.00000         3.36970         1.24397
H           0.00000         3.36970        -1.24397
H           0.00000         1.23670        -2.48620
H          -0.00000        -1.23670         2.48620
H          -0.00000        -1.23670        -2.48620
H          -0.00000        -3.36970        -1.24397
H          -0.00000        -3.36970         1.24397

