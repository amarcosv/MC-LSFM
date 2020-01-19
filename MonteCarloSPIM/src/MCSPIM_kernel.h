#pragma once

typedef struct MCSPIMStats {

    float simPhotons;		/**< Total simulated photons*/
    float detPhotons;		/**< Total detected photons*/
    unsigned int* focusedPhotons;	/**< Photons used to focus at each focusing plane */

    float mcxTime;		/**< Total MCX simulation time*/
    float focusTime;		/**< Total focusing time*/
    float simTime;		/**< Total simulation time*/
};

int simulationLauncher(int argc, char* argv[]);
