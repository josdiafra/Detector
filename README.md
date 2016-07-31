# josdiafra.github.io

Geant4 is a public toolkit for High Energy Physics (HEP) experiments using Object-Oriented environment and written in C++. GEANT4 is a Monte Carlo code not only for HEP but cosmic rays physics, space science and medical applications.

Any Geant4 application must consist of three class: DetectorConstruction, PrimaryGenerator and PhysicsList. On the one hand, DetectorConstruction class, allows us define the geometry and materials. Then, we use the PrimaryGenerator class to define the primary particles. In our case this particles will be americium and the differents uranium radioisotopes after.
Finally, PhysicsList allows us describe any physics processes. In this class we implant two physicsList package: EmStandardPhysics and RadioactiveDecay.

We use the remaining classes to obtain information about the energy deposited by the interesting particle or its track.
