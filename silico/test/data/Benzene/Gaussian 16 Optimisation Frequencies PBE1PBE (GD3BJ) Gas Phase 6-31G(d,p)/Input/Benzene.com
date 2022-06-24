%Chk="Benzene.chk"
%Rwf="Benzene.rwf"
%NProcShared=4
%Mem=10GB
#p Opt Freq PBE1PBE/6-31G(d,p) EmpiricalDispersion=(GD3BJ) Symmetry=(Tight) Population=(Regular) Density=(Current)

Benzene_Optimisation_Frequencies_PBE1PBE__GD3BJ__Gas_Phase_6_31G_d_p_

-1, 2
C           1.38269        -0.22179         0.00558
C           0.50635        -1.30704        -0.00813
C          -0.87140        -1.09058        -0.01457
C          -1.37323         0.21092        -0.00446
C          -0.49682         1.29597         0.01059
C           0.88099         1.07954         0.01412
H           2.45619        -0.39043         0.00959
H           0.89725        -2.32088        -0.01414
H          -1.55415        -1.93593        -0.02710
H          -2.44657         0.37942        -0.00832
H          -0.88770         2.31003         0.01913
H           1.56377         1.92487         0.02374

