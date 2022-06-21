%Chk="Pyridine.chk"
%Rwf="Pyridine.rwf"
%NProcShared=4
%Mem=10GB
#p Opt Freq PBE1PBE/6-31G(d,p) SCRF=(Solvent=Toluene) EmpiricalDispersion=(GD3BJ) Symmetry=(Tight) Population=(Regular) Density=(Current)

Pyridine_Optimisation_Frequencies_PBE1PBE__GD3BJ__Toluene_6_31G_d_p_

0, 1
C          -1.20591         0.69678         0.00001
C          -1.20638        -0.69598         0.00011
C          -0.00049        -1.39268        -0.00011
C           1.20594        -0.69675         0.00000
C           1.20638         0.69595         0.00012
N           0.00046         1.39270        -0.00012
H          -2.14642         1.24015        -0.00003
H          -2.14726        -1.23872        -0.00013
H          -0.00077        -2.47889         0.00014
H           2.14640        -1.24021        -0.00003
H           2.14724         1.23872        -0.00014

