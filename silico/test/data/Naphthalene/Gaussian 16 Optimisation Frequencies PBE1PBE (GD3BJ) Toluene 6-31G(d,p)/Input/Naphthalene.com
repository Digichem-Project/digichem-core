%Chk="Naphthalene.chk"
%Rwf="Naphthalene.rwf"
%NProcShared=4
%Mem=10GB
#p Opt Freq PBE1PBE/6-31G(d,p) SCRF=(Solvent=Toluene) EmpiricalDispersion=(GD3BJ) Symmetry=(Tight) Population=(Regular) Density=(Current)

Naphthalene_Optimisation_Frequencies_PBE1PBE__GD3BJ__Toluene_6_31G_d_p_

0, 1
C           1.22392        -1.39264        -0.00465
C           2.43349        -0.69589        -0.00607
C           2.43463         0.69593        -0.00478
C           1.22621         1.39466        -0.00238
C           0.00050         0.70909        -0.00158
C          -0.00066        -0.70506        -0.00249
C          -1.22636        -1.39063        -0.00154
C          -1.22407         1.39667        -0.00023
C          -2.43364         0.69992         0.00005
C          -2.43478        -0.69190        -0.00045
H           1.24325        -2.48004        -0.00559
H           3.37369        -1.24098        -0.00825
H           3.37573         1.23948        -0.00580
H           1.24733         2.48203        -0.00140
H          -1.24748        -2.47800        -0.00191
H          -1.24341         2.48407         0.00026
H          -3.37385         1.24501         0.00058
H          -3.37588        -1.23545        -0.00015

