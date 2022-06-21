%Chk="Pyridine.chk"
%Rwf="Pyridine.rwf"
%NProcShared=4
%Mem=10GB
#p TDA=(50-50, nstates=10) PBE1PBE/6-31G(d,p) SCRF=(Solvent=Toluene) EmpiricalDispersion=(GD3BJ) Symmetry=(Tight) Population=(Regular) Density=(Current) 6D 10F GFInput

Pyridine_Excited_States_TDA_10_Singlets_10_Triplets_PBE1PBE__GD3BJ__Toluene_6_31G_d_p_

0, 1
C          -1.13866        -0.71995         0.00002
C          -1.19540         0.67093        -0.00002
C          -0.00003         1.38179        -0.00002
C           1.19537         0.67097        -0.00002
C           1.13869        -0.71991        -0.00003
N           0.00002        -1.41721         0.00006
H          -2.05707        -1.30463         0.00004
H          -2.15405         1.17942        -0.00004
H          -0.00004         2.46782        -0.00005
H           2.15400         1.17951        -0.00002
H           2.05711        -1.30456         0.00004

