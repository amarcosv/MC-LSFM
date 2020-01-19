#pragma once
#include "focus_utils.h"
#include "cuda_runtime.h"


typedef float3 photonPos;
typedef float3 photonDir;
typedef float3 photon0Pos;


/**
*@brief Photon struct containing position, direction and weight
*/
typedef struct  __align__(8) photon {
	float w;			/**< photon weight*/
	photonPos pos;		/**< position of the photon - (x,y,z)*/
	photonDir dir;		/**< directional vector of the photon*/	
	photon0Pos Opos;	/**< launching position of the photon - (x0,y0,z0)*/
}photon;

/**
*@brief Struct containing the kernel focusing configuration
*/
typedef struct __align__(16) focusConfig {
	float f1;			/**< Lens 1 focal distance [grid Units]*/
	float lens1R;		/**< Lens 1 radius [grid Units]*/
	float NA1;			/**< Lens 1 Numerical aperture*/
	float f2;			/**< Lens 2 focal distance [grid Units]*/
	float lens2R;		/**< Lens 3 focal distance [grid Units]*/
	float NA2;			/**< Lens 1 Numerical aperture*/

	ulong2 sensorSize;	/**< Sensor dimensions in px*/
	float sensorFactor;	/**< Conversion factor from grid units to detector pixels */
	float focusZplane;	/**< Z plane to focus [grid Units]*/
	int focusZindex;

	float* illVolume;	/**< Pointer to illumination array. */
	float mua;			/**< Absorption coefficient of the medium [mm-1] */

	ulong3 mcx_simVolume;		/**< MCX volume size in grid units*/
	float mcx_pxSize;			/**< Grid unit size, MCX voxel size in mm/px*/

	uint3 ill_simVolume;		/**< Illumination volume size in px*/
	float ill_pxSize;			/**< Illumination volume unit size, MCX voxel size in mm/px*/

	unsigned long int Nphotons;	/**< Number of photons in dataFile*/
	unsigned int photonSize;	/**< Number of fields for a photon in dataFile*/

	int blockSize; /**< Number threads per block*/		
	int gridSize;  /**< Number of blocks in the grid*/
	int photonsPerThread;	/**< Number of photons to be processed per thread*/

}focusConfig;


int  list_gpu(SPIMConfig*cfg, SPIMGPUInfo **info);
//int focus_launcher(Config* cfg, GPUInfo* gpu, float *photonsData, float **imageVolume);
int focus_launcher(SPIMConfig* cfg, SPIMGPUInfo* gpu, float** imageVolume);
void memoryCalculator(focusConfig* hfcfg, int *nChunks, unsigned int *photonsPerChunk, unsigned int *remaindingPhotons);
void volumeAdd(SPIMConfig* cfg, SPIMGPUInfo* gpu, float* imageVolume);
int focus_mcx(SPIMConfig* cfg, SPIMGPUInfo* gpu, photon* photons, float** imageVolume);
void update_outputVolume(SPIMConfig* cfg, SPIMGPUInfo* gpu, float** outputVolume, float* imagevolume);
