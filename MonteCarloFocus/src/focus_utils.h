#pragma once
#include <stdio.h>
#include "cuda_runtime.h"


#define MAX_SESSION_LENGTH  256 
#define MAX_DEVICE  256 


/**
* @brief Configuration structure for Monte Carlo focusing
*/
typedef struct SPIMConfig {

	char* workingDir;		/**< Working directory path*/
	char* fileBaseName;		/**< Files basename*/
	char* configFileName;		/**< Configuration file name*/
	char* illFileName;		/**< Illumination file name*/
	char* outputFileName;		/**< Output file base name*/
	unsigned short int nFiles;	/**< Number of photon data files to process*/
	
	float mcx_pxSize;		/**< Grid unit size, MCX voxel size in mm/px*/
	uint3 mcx_simVolume;		/**< MCX volume size in px*/

	float ill_pxSize;		/**< Illumination volume unit size, MCX voxel size in mm/px*/
	uint3 ill_simVolume;		/**< Illumination volume size in px*/

	float f1;			/**< Lens 1 focal distance [mm]*/
	float f2;			/**< Lens 2 focal distance [mm]*/
	float NA1;			/**< Lens 1 numerical aperture*/
	float NA2;			/**< Lens 2 numerical aperture*/

	unsigned int Nphotons;		/**< Number of photons*/

	uint3 sensorSize;		/**< Camera sensor X,Y size in pixels*/
	float sensorpxSize;		/**< Camera sensor pixel size in mm/px*/
	
	float2 zScan;			/**< Z focus range [x,y]=[zStart,zStop] [mm]
								* z=0 at the base of the sample, zmax =mcx_simVolume z length * mcx_pxSize */
	unsigned int nzPlanes;		/**< Number of z positions to scan*/
	unsigned short int spimVol;	/**< Flag to process photons as a SPIM volume or OPT projection. 1 for SPIM, 0 for OPT projection */

	float* photonsData;		/**< Pointer to photons data array. [px py pz vx vy vz p0x p0y p0z w], with Nphotons values per field*/
	float* illVolume;		/**< Pointer to illumination array. */

	float mua;			/**< Absorption coefficient of the medium [mm-1] */

	/*For focus statistics*/
	unsigned int* focusedPhotons;		/**< Photons used to focus at each focusing plane */

	unsigned int photonSize;	/**< Number of fields for a photon in dataFile*/
	int gpuid;			/**<the ID of the GPU to use, starting from 1, 0 for auto*/
	char isgpuinfo;			/**<1 to print gpu info when attach, 0 do not print*/
	char deviceid[MAX_DEVICE];	/**<a 0-1 mask for all the GPUs, a mask of 1 means this GPU will be used*/
}SPIMConfig;

/**
* @brief Data structure for GPU information
*/
typedef struct SPIMGPUInfo {
	char name[MAX_SESSION_LENGTH];/**< name of the GPU */
	int id;                       /**< global index of the GPU, starting from 0 */
	int devcount;                 /**< total GPU count */
	int major;                    /**< major version of the CUDA device */
	int minor;                    /**< minor version of the CUDA device */
	size_t globalmem;             /**< size of the global memory in the GPU */
	size_t constmem;              /**< size of the constant memory in the GPU */
	size_t sharedmem;             /**< size of the shared memory in the GPU */
	int regcount;                 /**< size of the register file in the GPU */
	int clock;                    /**< clock in Hz of the GPU processor */
	int sm;                       /**< number of multi processors */
	int core;                     /**< number of stream processors */
	int autoblock;                /**< optimized number of blocks to launch */
	unsigned long autothread;     /**< optimized number of threads to launch */
	int maxgate;                  /**< max number of time gates that can be saved in one call */
	int maxmpthread;              /**< maximum thread number per multi-processor */



} SPIMGPUInfo;

#ifdef	__cplusplus
extern "C" {
#endif
void initcfg(SPIMConfig*cfg);
int importPhotonsData(const char *dataFile, SPIMConfig*cfg, float** photonsData);
int importIlluminationData(SPIMConfig*cfg, float** photonsData);
int exportVolume(float* imageVolume, SPIMConfig* cfg);
int readMchFile(const char* dataFile, SPIMConfig* cfg, float** photonsData);
int parseCommandLine(int argc, char* argv[], SPIMConfig* cfg);
int readConfigFile(SPIMConfig* cfg);
void buildMchPath(SPIMConfig* cfg, char* tmp, int fidx);
#ifdef	__cplusplus
}
#endif