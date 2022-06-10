%Chk="Naphthalene.chk"
%Rwf="Naphthalene.rwf"
%NProcShared=4
%Mem=10GB
#p TDA=(50-50, nstates=10) PBE1PBE/6-31G(d,p) SCRF=(Solvent=Toluene) EmpiricalDispersion=(GD3BJ) Symmetry=(Tight) Population=(Regular) Density=(Current) 6D 10F GFInput

Naphthalene_Excited_States_TDA_10_Singlets_10_Triplets_PBE1PBE__GD3BJ__Toluene_6_31G_d_p_

0, 1
C          -1.24046        -1.39914        -0.00000
C          -2.42600        -0.70664        -0.00000
C          -2.42600         0.70664         0.00000
C          -1.24046         1.39914         0.00000
C           0.00000         0.71423        -0.00000
C           0.00000        -0.71423        -0.00000
C           1.24046        -1.39914         0.00000
C           1.24046         1.39914        -0.00000
C           2.42600         0.70664        -0.00000
C           2.42600        -0.70664         0.00000
H          -1.23670        -2.48620        -0.00000
H          -3.36970        -1.24397         0.00000
H          -3.36970         1.24397        -0.00000
H          -1.23670         2.48620         0.00000
H           1.23670        -2.48620         0.00000
H           1.23670         2.48620        -0.00000
H           3.36970         1.24397         0.00000
H           3.36970        -1.24397        -0.00000

