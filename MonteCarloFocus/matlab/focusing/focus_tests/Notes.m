%% Notes ont the ouputs of MCX:
% The term "flux" has many meanings in physics and is indeed very confusing
% (I am sorry for choosing that ambiguous name in the documentation).
% 
% What MCX/MMC really outputs is the "fluence rate Phi(r,t)" (i.e. the first output of mcxlab.m).
% It has a unit of "W/m^2", referring to the expected number of photons, or
% Joule of energy passing by position r per unit time at time t, regardless of direction.
% This is also the TPSF generated at each voxel.
% 
% If you multiply the fluence rate by the time-gate width (cfg.tstep) and sum all
% the time gates, that yields "fluence" - Phi(r) (i.e. the so called CW fluence), 
% which has a unit of J/m^2. 
% 
% MCX can also output "current density" (J(r,t), unit W/m^2, same as Phi(r,t)) - 
% referring to the expected number of photons or Joule of energy flowing 
% through a unit area pointing towards a particular direction per unit time.
% The current density can be calculated at the boundary of the domain by two means:
% 
% 1. using the detected photon partial path output (i.e. the second output of mcxlab.m),
% one can compute the total energy E received by a detector, then one can 
% divide E by the area/aperture of the detector to obtain the J(r) at a detector
% (E should be calculated as a function of t by using the time-of-fly of detected
% photons, the E(t)/A gives J(r,t); if you integrate all time gates, the total E/A
% gives the current I(r), instead of the current density).
% 
% 2. use the new -X 1 or --saveref/cfg.issaveref option in mcx to enable the
% diffuse reflectance recordings on the boundary. the diffuse reflectance 
% is represented by the current density J(r) flowing outward from the domain.
% 
% The current density has, as mentioned, the same unit as fluence rate,
% but the difference is that J(r,t) is a vector, and Phi(r,t) is a scalar. Both measuring
% the energy flow across a small area (the are has direction in the case of J) per unit 
% time. 
% 
% You can find more rigorous definitions of these quantities in Lihong Wang's
% Biomedical Optics book, Chapter 5.
% 
% Here comes the confusion, because the current density sometimes is 
% also called radiative flux or flux (unit W/m^2) in some fields of physics, 
% so I labeled the output of mcx as flux in the documentation, but it was 
% indeed an incorrect choice of terminology - a flux is typically a vector,
% and a fluence/fluence rate is a scalar. 
% 
% In short, the first output of mcxlab is the fluence rate Phi(r,t), the second
% output detps can derive the current density J(r,t) at the detectors.
% 
% After normalization, the fluence rate output should match the fluence 
% rate derived from the RTE (in the case of large photons), and should also
% match nicely with the diffusion approximation in the diffusive regime
% when the DA is valid (i.e. high scattering). So, my effort of validation
% focuses primarily on comparing the MCX fluence/fluence rate output 
% with the diffusion equation's Green's function (assuming unitary energy 
% source) in the diffusive regime.
% 
% having these output definitions in mind, I think the observations you made
% seem to make sense to me.
% 
% if your medium has no scattering (like air or water), then, if you launch
% 1 Joule of photons at your source, any area/voxel in front of the source
% should receive exactly 1 J / area fluence, because a photon only passes by 
% a voxel one-way one time, and never returns. 
% 
% However, for a scattering medium, a photon can pass by a voxel wall 
% many times due to scattering. Thus, the total net energy flowing through 
% any voxel boundary is more than what it is supposed to receive when
% there is no scattering.

%SOURCE: https://groups.google.com/forum/#!topic/mcx-users/N55J79w9FYY

% output='energy' records the total absorbed energy in each voxel, if you launch 1 joule of energy, how much of those ended up inside each voxel. If you add the total energy output together, this number should roughly equal to the percentage absorption reported at the end of the log.
% 
% output='fluence' is the time resolved fluence (i.e. integrated of each time gate) at each voxel. it has a unit of 1/mm^2 or joule/mm^2. the relationship between fluence output and the energy deposit output is that 
% 
% energy-deposit=fluence*mua*voxel_volume = (fluence_rate*time_gate)*mua*voxel_volume


%% 