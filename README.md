# Dipole moment & wavefunction in strong-field approximation

by Wei-Chun Chu (2015)

This code is to calculate dipole moment and wave function of simple atomic systems driven by an arbitrary linearly polarized electric field. The model is based on the semi-classical strong-field approximation (SFA) by Lewenstein et. al. [M. Lewenstein, Ph. Balcou, M. Yu. Ivanov, Anne L'Huillier, and P. B. Corkum, Phys. Rev. A 49, 2117 (1994)]. The original paper has a few errors on the plus/minus sign. Read [W.-C. Chu, J. C. Travers, and P. St.J. Russell, New J. Phys. 18, 023018 (2016)].

Simply compile sfaq.tex with latex to read the short description of the algorithm used in the code.

All codes are written in FORTRAN 90 and OpenMP. Compile all executables with _make executable_. Example configuration files are in the folder "examples".
