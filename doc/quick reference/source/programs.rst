Choosing a Program
==================

In the Zysman-Colman group, the characterisation of organic and organometallic emitters or photocatalysts is typically performed using Gaussian.
When choosing between Gaussian 09 and Gaussian 16, the latter, more up-to-date program should generally be preferred.
However, it is important to note that calculations performed with different version of Gaussian cannot be directly compared (because differences between the two programs result in differences in the computed data). Hence when carrying out computations as part of a series, all molecules under investigation should be performed using the same program to ensure the results are comparable.


Turbomole and Multi-Resonant Thermally Activated Delayed Fluorescence
---------------------------------------------------------------------

The excited states of multi-resonant thermally activated delayed fluorescence (MR-TADF) emitters cannot be accurately modelled using conventional DFT.
Instead, these emitters require a class of 'higher order' computational methods to be modelled effectively. 
called coupled cluster (CC) methods
(although see SCS-ADC(2)). Coupled cluster methods come in various flavours, some of which are
implemented by both Gaussian and Turbomole, such as CCD, CCSD and CCSD(T), but the
majority of these methods are extremely expensive (very slow). To overcome this, the approximate
coupled-cluster singles and doubles model (CC2) uses a number of approximations to drastically
increase the speed of computation, which makes it viable for (semi-) rapid characterisation. This
CC2 method is only available in Turbomole, and so all MR-TADF type emitters must use
Turbomole for the calculation of singlet and triplet energies.
Although faster than most CC methods, CC2 is still extremely slow compared to DFT. Therefore, it
is recommended to first perform a DFT optimisation so the starting geometry can be brought closer
to the final geometry, reducing the overall computation time. It is recommended that this DFT step
be performed in Gaussian, so that calculated energies can be compared to other, non MR-TADF
emitters.