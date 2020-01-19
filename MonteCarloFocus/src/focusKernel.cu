#include "focusKernel.h"
#include "focus_utils.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_occupancy.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CUDA_ERROR_CHECK
#ifdef CUDA_DEBUG
#define PHOTON_DEBUG(x) print_photon(x);
#define DEBUG(x) printf x;
#else
#define PHOTON_DEBUG(x) 
#define DEBUG(x) 
#endif


#define PRINT(fp,...) fprintf(fp,__VA_ARGS__)

__constant__ focusConfig fcfg[1];
const int LOAD_BUFFER_SIZE = 2097152;
//const int LOAD_BUFFER_SIZE = 3;
/**
*@brief Error Handler for Cuda Api calls
*/
#define CUDA_ASSERT( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
inline void __cudaSafeCall(cudaError err, const char *file, const int line)
{
#ifdef CUDA_ERROR_CHECK
	if (cudaSuccess != err)
	{
		fprintf(stderr, "cudaSafeCall() failed at %s:%i : %s\n",
			file, line, cudaGetErrorString(err));
		exit(-1);
	}
#endif

	return;
}

/**
*@brief Error Handler for Kernel calls
*/
#define CUDA_ERROR_CHECK()    __cudaCheckError( __FILE__, __LINE__ )
inline void __cudaCheckError(const char *file, const int line)
{
#ifdef CUDA_ERROR_CHECK
	cudaError err = cudaGetLastError();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
			file, line, cudaGetErrorString(err));
		exit(-1);
	}

	// More careful checking. However, this will affect performance.
	// Comment away if needed.
	err = cudaDeviceSynchronize();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
			file, line, cudaGetErrorString(err));
		exit(-1);
	}
#endif

	return;
}

/**
*@brief Copy Array to GPU for Kernel calls
*arguments: (**device,*host,number of elements, data type)
*/
#define CUDA_LOAD(device,host,elements,type) CopyArrayToGPU((void**)device,(void*)host,elements, sizeof(type))
int CopyArrayToGPU(void** DeviceArray, void* HostArray, unsigned long int NumElements, size_t typelength)
{
    unsigned long int bytes = typelength * NumElements;

    // Allocate memory on the GPU for array
    if (cudaMalloc(DeviceArray, bytes) != cudaSuccess)
    {
	printf("CopyArrayToGPU(): Couldn't allocate mem for array on GPU.");
	return 1;
    }

    // Copy the contents of the host array to the GPU
    if (cudaMemcpy(*DeviceArray, HostArray, bytes, cudaMemcpyHostToDevice) != cudaSuccess)
    {
	printf("CopyArrayToGPU(): Couldn't copy host array to GPU.");
	cudaFree(*DeviceArray);
	return 1;
    }

    return 0;
}


/**
*@brief Copy GPU Array to Host for Kernel calls
*arguments: (**host,*device,elements, type)
*/
#define CUDA_UNLOAD(host,device,elements,type) CopyArrayToHost((void**)host,(void*)device,elements, sizeof(type))
int CopyArrayToHost(void **HostArray, void *DeviceArray, unsigned long int NumElements, size_t typelength)
{
	unsigned long int bytes = typelength * NumElements;

	// Copy the contents of the host array to the GPU
	if (cudaMemcpy(*HostArray, DeviceArray, bytes, cudaMemcpyDeviceToHost) != cudaSuccess)
	{
		printf("CopyArrayToGPU(): Couldn't copy device array to host.");
		cudaFree(*HostArray);
		return 1;
	}

	return 0;
}


/**
*@brief Build focus settings structure to be copied in constant memory
*/
void createfocusConfig(SPIMConfig* cfg, SPIMGPUInfo* gpu, focusConfig* fcfg) {
	fcfg->f1 = cfg->f1/cfg->mcx_pxSize; 
	fcfg->f2 = cfg->f2 /cfg->mcx_pxSize;
	fcfg->NA1 = cfg->NA1;
	fcfg->NA2 = cfg->NA2;
		//fcfg->lens1R = cfg->mcx_simVolume.x / 2.f ;
	//fcfg->lens2R = fcfg->lens1R;
	fcfg->lens1R = fcfg->f1 * tanf(asinf(0.99f * fcfg->NA1));
	fcfg->lens2R = fcfg->f2 * tanf(asinf(0.99f * fcfg->NA2));
	fcfg->sensorFactor = cfg->sensorpxSize/cfg->mcx_pxSize;
	fcfg->sensorSize.x = cfg->sensorSize.x;
	fcfg->sensorSize.y = cfg->sensorSize.y;
	fcfg->focusZplane=cfg->zScan.x / cfg->mcx_pxSize; /** Convert to Grid units*/
	fcfg->mcx_simVolume.x = cfg->mcx_simVolume.x;
	fcfg->mcx_simVolume.y = cfg->mcx_simVolume.y;
	fcfg->mcx_simVolume.z = cfg->mcx_simVolume.z;
	fcfg->mcx_pxSize = cfg->mcx_pxSize;
	fcfg->ill_simVolume.x = cfg->ill_simVolume.x;
	fcfg->ill_simVolume.y = cfg->ill_simVolume.y;
	fcfg->ill_simVolume.z = cfg->ill_simVolume.z;
	fcfg->ill_pxSize = cfg->ill_pxSize;
	fcfg->Nphotons = cfg->Nphotons;
	fcfg->blockSize = gpu->autoblock;   /**< The launch configurator returned block size*/
	fcfg->gridSize = gpu->autothread / fcfg->blockSize;    /**< The actual grid size needed, based on input size*/
	fcfg->photonsPerThread = (unsigned long int)ceil((float)cfg->Nphotons / (float)gpu->autothread);	
	fcfg->photonSize = cfg->photonSize;
	fcfg->mua = cfg->mua;

	if (cfg->illVolume) {
		float *gillVolume;
		//if (cfg->ill_simVolume.z!=1)
		//CUDA_LOAD(&gillVolume, cfg->illVolume, cfg->ill_simVolume.x * cfg->ill_simVolume.y * (2 * cfg->ill_simVolume.z + 1));
		//else
		    CUDA_LOAD(&gillVolume, cfg->illVolume, cfg->ill_simVolume.x * cfg->ill_simVolume.y * cfg->ill_simVolume.z,float);
		fcfg->illVolume = gillVolume;
		//fcfg->illVolume = NULL;
	}


}
#define FOCUS_AT(x,y) set_focusDistance<<<1,1>>>(x, y); cudaDeviceSynchronize();
/**
*@brief Sets the focal plane z distance for the focusing kernel
*
*@param[in] focusZplane:  focal plane in grid units
*@param[in] focusZinfex:  index at output volume to store image
*/
__global__ void set_focusDistance(float focusZplane, int focusZindex) {
	fcfg->focusZplane = focusZplane;
	fcfg->focusZindex = focusZindex;
}

#define SET_PHOTONS_DATA(x,y,z) set_photons<<<1,1>>>(x, y, z); cudaDeviceSynchronize();
/**
*@brief Sets the focal plane z distance for the focusing kernel [d in grid units]
*/
__global__ void set_photons(unsigned int photons, float mcx_pxSize, int photonSize) {
	fcfg->Nphotons = photons;
	fcfg->mcx_pxSize = mcx_pxSize;
	fcfg->photonSize = photonSize;
}

/**
*@brief prints photon position information
*/
__global__ void print_photonPosition(photon *photons, unsigned long index) {
	printf("Photon %d: pos[%f, %f, %f]\n", index, photons[index].pos.x, photons[index].pos.y, photons[index].pos.z);
}

__global__ void print_focusConfigGPU() {
	printf("Focus Configuration:\n\tFocus1=%f\n\tFocus2=%f\n\tNphotons=%d\n", fcfg->f1, fcfg->f2, fcfg->Nphotons);
	printf("\tmua=%f\n\tFocus2=%d\n\tPhperthread=%d\n", fcfg->mua, fcfg->photonSize, fcfg->photonsPerThread);
	return;
}



#define PHOTON_PRINT(ph,id) print_photonData<<<1,1>>>(ph,id); cudaDeviceSynchronize();
/**
*@brief prints photon position information
*/
__global__ void print_photonData(photon *photons,  unsigned long index) {
	printf("Photon %d:\n", index);
	printf("\t pos: [%f, %f, %f]\n", photons[index].pos.x, photons[index].pos.y, photons[index].pos.z);
	printf("\t dir: [%f, %f, %f]\n", photons[index].dir.x, photons[index].dir.y, photons[index].dir.z);
	printf("\t pos0: [%f, %f, %f]\n", photons[index].Opos.x, photons[index].Opos.y, photons[index].Opos.z);
	printf("\t w:	 %e\n",  photons[index].w);
	
}

/**
*@brief prints photon position information
*/
__device__ void print_photon(photon *photons) {
	printf("\t pos: [%f, %f, %f]\n", photons->pos.x, photons->pos.y, photons->pos.z);
	printf("\t dir: [%f, %f, %f]\n", photons->dir.x, photons->dir.y, photons->dir.z);
	printf("\t w:	 %e\n", photons->w);

}

/**
*@brief Calculates photon weight and shifts photon coordinates to have the (0,0) at the center of the volume 
*
* @param[in, out] photons:  photon data structure
*
*/
__global__ void prepare_MCXphotonData(photon* photons) {
    size_t thread = blockIdx.x * blockDim.x + threadIdx.x;

    if (thread < fcfg->Nphotons) {
	photons[thread].pos.x = photons[thread].pos.x - fcfg->mcx_simVolume.x / 2.0f;
	photons[thread].pos.y = photons[thread].pos.y - fcfg->mcx_simVolume.y / 2.0f;
	photons[thread].Opos.x = photons[thread].Opos.x - fcfg->mcx_simVolume.x / 2.0f;
	photons[thread].Opos.y = photons[thread].Opos.y - fcfg->mcx_simVolume.y / 2.0f;
	photons[thread].w = expf(-photons[thread].w  * fcfg->mua); /** Both are in grid units!*/
    }

}


/**
*@brief Sets photon data in photons array from the raw photonsArray
*
*@param[in,out] photons: structure of photon data
*@param[in]	photonsData: gpu buffer with raw photons info [id ppath(M) p(3) v(3) p0(3) path]
*@param[in]	Nphotons:	 number of photons to be processed
*@param[in]	loadOffset:	 offset
*
*/
__global__ void set_photonData(photon* photons, float* photonsData, unsigned int Nphotons, int offset) {
	size_t thread = blockIdx.x * blockDim.x + threadIdx.x;
	int limit = thread * fcfg->photonsPerThread + fcfg->photonsPerThread;
	if (limit > Nphotons) limit = Nphotons ;
	photons = photons + offset;

	for (int i = thread * fcfg->photonsPerThread; i < limit; i++) {
		//printf("Thread# %d: processing photonData %d into photon %d\n", (int)thread, i, i * fcfg->photonSize);
		photons[i].pos.x = photonsData[i * fcfg->photonSize + 3] -fcfg->mcx_simVolume.x / 2.0f;
		photons[i].pos.y = photonsData[i * fcfg->photonSize + 4] -fcfg->mcx_simVolume.y / 2.0f;
		photons[i].pos.z = photonsData[i * fcfg->photonSize + 5];
		photons[i].dir.x = photonsData[i * fcfg->photonSize + 6];
		photons[i].dir.y = photonsData[i * fcfg->photonSize + 7];
		photons[i].dir.z = photonsData[i * fcfg->photonSize + 8];
		photons[i].Opos.x = photonsData[i * fcfg->photonSize + 9];
		photons[i].Opos.y = photonsData[i * fcfg->photonSize + 10];
		photons[i].Opos.z = photonsData[i * fcfg->photonSize + 11];
		photons[i].w = expf(-photonsData[i * fcfg->photonSize + 2]*fcfg->mcx_pxSize * fcfg->mua);
		//photons[i].w = photonsData[i * fcfg->photonSize + 2] * 0.025f;
	}
	return;
}

/**
*@brief Sets position element in photons array from the raw photonsArray
*/
__global__ void set_photonPosition(photon *photons, float *photonsData, unsigned long int Nphotons) {
	size_t thread = blockIdx.x * blockDim.x + threadIdx.x;	
	size_t limit = thread * fcfg->photonsPerThread + fcfg->photonsPerThread;
	if (limit > Nphotons) limit = Nphotons; /* Prevent from accesing array out of bounds photons*/
	for (int i = thread * fcfg->photonsPerThread; i < limit; i++) {
		photons[i].pos.x = photonsData[i] - fcfg->mcx_simVolume.x / 2.0f;
		photons[i].pos.y = photonsData[i + Nphotons] - fcfg->mcx_simVolume.y / 2.0f;
		photons[i].pos.z = photonsData[i + 2 * Nphotons];
	}
	return;
}

/**
*@brief Sets position element in photons array from the raw photonsArray
*/
/**__global__ void set_photonPosition(photon *photons,float *photonsData, unsigned long int Nphotons){
	size_t thread = blockIdx.x * blockDim.x + threadIdx.x;
	if (thread < Nphotons) {
		photons[thread].pos.x = photonsData[thread] - fcfg->mcx_simVolume.x/2.0f;
		photons[thread].pos.y = photonsData[thread + Nphotons] - fcfg->mcx_simVolume.y / 2.0f;
		photons[thread].pos.z = photonsData[thread + 2 * Nphotons] ;
	}
	return;
}*/
/**
*@brief Sets direction element in photons array from the raw photonsArray
*/
__global__ void set_photonDirection(photon *photons, float *photonsData, unsigned long int Nphotons) {
	size_t thread = blockIdx.x * blockDim.x + threadIdx.x;	
	size_t limit = thread * fcfg->photonsPerThread + fcfg->photonsPerThread;
	if (limit > Nphotons) limit = Nphotons;
	for (int i = thread * fcfg->photonsPerThread; i < limit; i++) {
		photons[i].dir.x = photonsData[i];
		photons[i].dir.y = photonsData[i + 1 * Nphotons];
		photons[i].dir.z = photonsData[i + 2 * Nphotons];
	}
	return;
}

/**
*@brief Sets launch position in photons array from the raw photonsArray
*/
__global__ void set_photonOrigin(photon0Pos *photonsOrigin, float *photonsData, unsigned long int Nphotons) {
	size_t thread = blockIdx.x * blockDim.x + threadIdx.x;	
	size_t limit = thread * fcfg->photonsPerThread + fcfg->photonsPerThread;
	if (limit > Nphotons) limit = Nphotons;
	for (int i = thread * fcfg->photonsPerThread; i < limit; i++) {
		photonsOrigin[i].x = __float2int_rn(photonsData[i]);
		photonsOrigin[i].y = __float2int_rn(photonsData[i + 1 * Nphotons]);
		photonsOrigin[i].z = __float2int_rn(photonsData[i + 2 * Nphotons]);
	}
	return;
}

/**
*@brief Sets photonWeight value in photons array from the raw photonsArray
*/
__global__ void set_photonWeight(photon *photons, float *photonsData, unsigned long int Nphotons) {
	size_t thread = blockIdx.x * blockDim.x + threadIdx.x;	
	size_t limit = thread * fcfg->photonsPerThread + fcfg->photonsPerThread;
	if (limit > Nphotons) limit = Nphotons;
	for (int i = thread * fcfg->photonsPerThread; i < limit; i++) {
		photons[i].w = photonsData[i];
	}
	return;
}



/**
*@brief Builds the photonsStruct using the raw data from the .mch binary file.
*/
int load_photonStruct(focusConfig* hfcfg, SPIMGPUInfo* gpu, unsigned long int Nphotons, float* photonsData, photon** photons) {
	float* gphotonsData = NULL;	/**< Device buffer to store temporally raw photon data*/
	/*Prepare GPU settings*/
	int blockSize = gpu->autoblock;   /**< The launch configurator returned block size*/
	int gridSize = gpu->autothread/gpu->autoblock;    /**< The actual grid size needed, based on input size*/
	int nChunks = Nphotons / LOAD_BUFFER_SIZE;
	int remaining= Nphotons % LOAD_BUFFER_SIZE;
	int loadBuffer = LOAD_BUFFER_SIZE;
	if (remaining > 0) ++nChunks;

	printf("Running memory allocation kernel for %d photons.\nAllocation in %d steps of %d photons per chunk\n", Nphotons, nChunks, LOAD_BUFFER_SIZE);	
	/*Allocate photons structures*/
	CUDA_ASSERT(cudaMalloc((void**)photons, Nphotons * sizeof(photon)));
	//CUDA_ASSERT(cudaMalloc((void**)photonsOrigin, Nphotons * sizeof(photon0Pos)));
	CUDA_ASSERT(cudaMalloc(&gphotonsData, LOAD_BUFFER_SIZE * hfcfg->photonSize * sizeof(float)));
	//blockSize = 64;
		//gridSize = 2560;

	int offset = 0;
	int step = 0;
	while (offset < Nphotons) {
		if (offset + LOAD_BUFFER_SIZE > Nphotons) loadBuffer = Nphotons - offset;
		printf("\tStep %d of %d: Loading photons  %d to %d.\n", step+1,nChunks,offset, offset+ loadBuffer);
		CUDA_ASSERT(cudaMemcpy(gphotonsData, photonsData + offset * hfcfg->photonSize, loadBuffer * hfcfg->photonSize * sizeof(float), cudaMemcpyHostToDevice));
		set_photonData << <gridSize, blockSize >> > (*photons, gphotonsData, loadBuffer, offset);
		cudaDeviceSynchronize();
		CUDA_ERROR_CHECK();
		offset += LOAD_BUFFER_SIZE;
		step++;
	}

	CUDA_ASSERT(cudaFree(gphotonsData));
	return 0;
}







/**
*@brief Builds the photonsStruct using the raw data from the binary file. 
*/
int build_photonStruct(unsigned long int Nphotons, float *photonsData, photon **photons, photon0Pos **photonsOrigin){
	float* gphotonsData = NULL;	/**< Device buffer to store temporally raw photon data*/ 
	/*Prepare GPU settings*/
	int blockSize=128;   /**< The launch configurator returned block size*/ 
	int gridSize=12288/blockSize;    /**< The actual grid size needed, based on input size*/ 
	
	//gridSize = (Nphotons + blockSize - 1) / blockSize;
	printf("Running memory allocation kernel for %d photons. %d blocks of %d threads per block\n", Nphotons, gridSize, blockSize);

	/*Allocate photons structure*/
	CUDA_ASSERT(cudaMalloc((void **)photons, Nphotons*sizeof(photon)));
	CUDA_ASSERT(cudaMalloc((void **)photonsOrigin, Nphotons*sizeof(photon0Pos)));


	/*Assign photons position*/
	CUDA_LOAD(&gphotonsData, photonsData, 3 * Nphotons, float);
	set_photonPosition << <gridSize, blockSize >> > (*photons, gphotonsData,Nphotons);
	cudaDeviceSynchronize();
	CUDA_ERROR_CHECK();

	/*Assign photons direction*/
	CUDA_LOAD(&gphotonsData, (photonsData + 3 * Nphotons), 3 * Nphotons, float);
	set_photonDirection << <gridSize, blockSize >> > (*photons, gphotonsData, Nphotons);
	cudaDeviceSynchronize();
	CUDA_ERROR_CHECK();

	/*Assign photons launch position*/
	CUDA_LOAD(&gphotonsData, (photonsData + 6 * Nphotons), 3 * Nphotons,float);
	set_photonOrigin << <gridSize, blockSize >> > (*photonsOrigin, gphotonsData, Nphotons);
	cudaDeviceSynchronize();
	CUDA_ERROR_CHECK();

	CUDA_ASSERT(cudaFree(gphotonsData));
	/*Assign photons launch position*/
	CUDA_LOAD(&gphotonsData, (photonsData + 9 * Nphotons), Nphotons, float);
	set_photonWeight << <gridSize, blockSize >> > (*photons, gphotonsData, Nphotons);
	cudaDeviceSynchronize();
	CUDA_ERROR_CHECK();


	CUDA_ASSERT(cudaFree(gphotonsData));
	return 0;
}

/**
*@brief Build camera detector buffer.
*/
int build_cameraSensor(SPIMConfig *cfg, float **cameraSensor) {
	CUDA_ASSERT(cudaMalloc((void **)cameraSensor, (cfg->sensorSize.x * cfg->sensorSize.y) * cfg->nzPlanes *  sizeof(float)));
	CUDA_ASSERT(cudaMemset(*cameraSensor, 0, (cfg->sensorSize.x * cfg->sensorSize.y) * cfg->nzPlanes * sizeof(float)));
	return 0;
}

/**
*@brief Propagates photon a distance z
*/
__device__ void photonProp(photon* photon, float z) {
	float t = fdividef(z , photon->dir.z);
	photon->pos.x += t * photon->dir.x;
	photon->pos.y += t * photon->dir.y;
}

/**
*@brief Refract photon according to the lens focal distance
*/
__device__ void photonRefract(photon* photon, float f) {

	photonDir dirP;		/** photon direction in photon's coordinate system*/
	photonDir dirPr;	/** rotated photon direction in photon's coordinate system*/
	float2 P;
	float h;	
	float sintheta;
	float costheta;
	/**Normalize photon base vector*/	
	h = rsqrtf(photon->pos.x*photon->pos.x + photon->pos.y*photon->pos.y);
	P.x = photon->pos.x*h;
	P.y = photon->pos.y*h;	
	sincosf(1.0f/(f*h),&sintheta,&costheta);	
	/** Change photon coordinate system*/
	dirP.x = photon->dir.x*P.x + photon->dir.y*P.y;
	dirP.y = photon->dir.y*P.x - photon->dir.x*P.y;
	dirP.z = photon->dir.z;	
	/** Rotate photon*/
	dirPr.x = dirP.x*costheta + dirP.z*sintheta;
	dirPr.y = dirP.y;
	dirPr.z = -dirP.x*sintheta + dirP.z*costheta;	
	/**Revert to system, coordinate system*/
	photon->dir.x = P.x*dirPr.x - P.y*dirPr.y;
	photon->dir.y = P.y*dirPr.x + P.x*dirPr.y;
	photon->dir.z = dirPr.z;

}

/**
*@brief Check if a photon passes through a puppil
*True if photon passes through the lens 
*/
__device__ int checkPupil(photon* photon, float r) {
	if ((photon->pos.x*photon->pos.x + photon->pos.y*photon->pos.y) > (r*r))
		return 0;
	return 1;
}

/**
*@brief Check if a photon passes through a puppil
*True if photon passes through the lens
*/
__device__ int checkLensAcceptance(photon* photon, float NA) {
	float mod = sqrtf(photon->dir.x* photon->dir.x+ photon->dir.y* photon->dir.y+ photon->dir.z* photon->dir.z);
	float sintheta = sinf(acosf(photon->dir.z/mod));
	if (sintheta > NA)
		return 0;
	return 1;
}

/**
*@brief Detect a photon hitting the camera sensor
*/
__device__ int detectPhoton(photon* photon, float* cameraSensor) {
	uint2 cameraPixel;
	float2 sensorLandPos;
	int weightIndex;
	int weightOffset;
	/** Convert to detector units*/
	sensorLandPos.x = (photon->pos.x )/fcfg->sensorFactor + fdividef(fcfg->sensorSize.x, 2.0f);
	sensorLandPos.y = (photon->pos.y )/fcfg->sensorFactor + fdividef(fcfg->sensorSize.y, 2.0f);
	DEBUG(("\tLanding at [%f %f]", sensorLandPos.x, sensorLandPos.y));
	/** Check detector out of bounds*/
	if ((sensorLandPos.x < 0.f )|| (sensorLandPos.y < 0.f))
		return 0;
	cameraPixel.x = __float2uint_rn(sensorLandPos.x);
	cameraPixel.y = __float2uint_rn(sensorLandPos.y);
	if ((cameraPixel.x > (fcfg->sensorSize.x-(unsigned int)1)) || (cameraPixel.y > (fcfg->sensorSize.y-(unsigned int)1)))
		return 0;
	
	/** Add weight*/
	DEBUG(("\t\tDetected at [%d %d] IDX: %d\n", cameraPixel.x, cameraPixel.y, cameraPixel.y*fcfg->sensorSize.x + cameraPixel.x));
	//unsigned long volumeIndex = cameraPixel.y * fcfg->sensorSize.x + cameraPixel.x + fcfg->focusZindex * fcfg->sensorSize.x * fcfg->sensorSize.y;


		//***********************
		//weightIndex = photon->Opos.z * fcfg->mcx_simVolume.x * fcfg->mcx_simVolume.y + photon->Opos.y * fcfg->mcx_simVolume.x + photon->Opos.x;
		//weightOffset = ((fcfg->mcx_simVolume.z-1) - __float2uint_rn(fcfg->focusZplane)) * fcfg->mcx_simVolume.x * fcfg->mcx_simVolume.y;
		//atomicAdd(&cameraSensor[cameraPixel.x + cameraPixel.y * fcfg->sensorSize.x + fcfg->focusZindex * fcfg->sensorSize.x * fcfg->sensorSize.y], photon->w * fcfg->illVolume[weightIndex + weightOffset]);

		//***********************

		unsigned int ill_focusZplane = __float2uint_rn(fcfg->focusZplane * fdividef(fcfg->mcx_pxSize,fcfg->ill_pxSize));
		weightOffset = ((fcfg->ill_simVolume.z /2) - ill_focusZplane) * fcfg->ill_simVolume.x * fcfg->ill_simVolume.y;

		uint3 ill_0pos;
		ill_0pos.x = __float2uint_rn(photon->Opos.x * fdividef(fcfg->mcx_pxSize, fcfg->ill_pxSize));		
		ill_0pos.z = __float2uint_rn(photon->Opos.z * fdividef(fcfg->mcx_pxSize, fcfg->ill_pxSize));

		if (fcfg->ill_simVolume.y != 1)
		    ill_0pos.y = __float2uint_rn(photon->Opos.y * fdividef(fcfg->mcx_pxSize, fcfg->ill_pxSize));
		else 
		    ill_0pos.y = 0;

		//weightIndex = photon->Opos.z * fcfg->ill_simVolume.x * fcfg->ill_simVolume.y + photon->Opos.y * fcfg->ill_simVolume.x + photon->Opos.x;
		weightIndex = ill_0pos.z * fcfg->ill_simVolume.x * fcfg->ill_simVolume.y + ill_0pos.y * fcfg->ill_simVolume.x + ill_0pos.x;
		//printf("Photon launched from[%f %f %f] -> [%u %u %u]\n Focus at %f Offset= %d Index= %d ILLIDX= %d\n", photon->Opos.x, photon->Opos.y, photon->Opos.z, ill_0pos.x, ill_0pos.y, ill_0pos.z, fcfg->focusZplane, weightOffset,weightIndex, weightIndex + weightOffset);
		//DEBUG(("\t\tPhoton launched from [%f %f %f] IDX: %d\n", photon->Opos.x, photon->Opos.y, photon->Opos.z,cameraPixel.y * fcfg->sensorSize.x + cameraPixel.x));
		if ((weightIndex + weightOffset) >= 0 && (weightIndex + weightOffset) < (fcfg->ill_simVolume.x * fcfg->ill_simVolume.y * fcfg->ill_simVolume.z)) {
		    atomicAdd(&cameraSensor[cameraPixel.x + cameraPixel.y * fcfg->sensorSize.x + fcfg->focusZindex * fcfg->sensorSize.x * fcfg->sensorSize.y], photon->w * fcfg->illVolume[weightIndex + weightOffset]);
		    return 1;
		}
		//***********************
	return 0;
}
/**
*@brief Focusing Kernel for 4F system
*/
__global__ void focus4F(photon *photons, float* cameraSensor) {
	size_t thread = blockIdx.x * blockDim.x + threadIdx.x;
	
	size_t limit = thread * fcfg->photonsPerThread + fcfg->photonsPerThread;
	
	if (limit > fcfg->Nphotons) limit = fcfg->Nphotons;
	for (int i = thread * fcfg->photonsPerThread; i < limit; i++) {
		//printf("focusing Photon # [%d]", i);
		photon threadphoton = photons[i];
		DEBUG(("focusing at %f  px\n", fcfg->focusZplane));
		DEBUG(("Photon #%d\t", thread));
		PHOTON_DEBUG(&threadphoton);
		/**Propagate to Lens 1 */
		//photonProp(&threadphoton, fcfg->f1 - fcfg->focusZplane);
		photonProp(&threadphoton, fcfg->f1 - (fcfg->mcx_simVolume.z - fcfg->focusZplane));
		DEBUG(("Photon #%d\t", thread));
		PHOTON_DEBUG(&threadphoton);
		/** Lens 1 refraction*/
		if (checkPupil(&threadphoton, fcfg->lens1R) == 0) continue;
		if(checkLensAcceptance(&threadphoton, fcfg->NA1)==0) continue;
		photonRefract(&threadphoton, -fcfg->f1);
		DEBUG(("Photon #%d\t", thread));
		PHOTON_DEBUG(&threadphoton);
		/**Propagate to Lens 2 */
		photonProp(&threadphoton, fcfg->f1 + fcfg->f2);
		DEBUG(("Photon #%d\t", thread));
		PHOTON_DEBUG(&threadphoton);
		/** Lens 2 refraction*/
		if (checkPupil(&threadphoton, fcfg->lens2R) == 0) continue;
		photonRefract(&threadphoton, -fcfg->f2);
		DEBUG(("Photon #%d\t", thread));
		PHOTON_DEBUG(&threadphoton);
		/** Propagate to detector */
		photonProp(&threadphoton, fcfg->f2);
		DEBUG(("Photon #%d\t", thread));
		PHOTON_DEBUG(&threadphoton);
		/** Detect photon*/
		DEBUG(("Photon #%d\t", thread));
		detectPhoton(&threadphoton, cameraSensor);

	}
}

/**
*@brief Focusing Kernel for 4F system
*/
__global__ void focus4F_auto(photon* photons, float* cameraSensor,unsigned int* gfocused) {
    size_t thread = blockIdx.x * blockDim.x + threadIdx.x;
    int res;
   // size_t limit = thread * fcfg->photonsPerThread + fcfg->photonsPerThread;
    if (thread < fcfg->Nphotons) {	
	    //printf("focusing Photon # [%d]", i);
	    photon threadphoton = photons[thread];
	    DEBUG(("focusing at %f  px\n", fcfg->focusZplane));
	    DEBUG(("Photon #%d\t", thread));
	    PHOTON_DEBUG(&threadphoton);
	    /**Propagate to Lens 1 */
	    //photonProp(&threadphoton, fcfg->f1 - fcfg->focusZplane);
	    photonProp(&threadphoton, fcfg->f1 - (fcfg->mcx_simVolume.z - fcfg->focusZplane));
	    DEBUG(("Photon #%d\t", thread));
	    PHOTON_DEBUG(&threadphoton);
	    /** Lens 1 refraction*/
	    if (checkPupil(&threadphoton, fcfg->lens1R) == 0) return;
	    if (checkLensAcceptance(&threadphoton, fcfg->NA1) == 0) return;
	    photonRefract(&threadphoton, -fcfg->f1);
	    DEBUG(("Photon #%d\t", thread));
	    PHOTON_DEBUG(&threadphoton);
	    /**Propagate to Lens 2 */
	    photonProp(&threadphoton, fcfg->f1 + fcfg->f2);
	    DEBUG(("Photon #%d\t", thread));
	    PHOTON_DEBUG(&threadphoton);
	    /** Lens 2 refraction*/
	    if (checkPupil(&threadphoton, fcfg->lens2R) == 0) return;
	    photonRefract(&threadphoton, -fcfg->f2);
	    DEBUG(("Photon #%d\t", thread));
	    PHOTON_DEBUG(&threadphoton);
	    /** Propagate to detector */
	    photonProp(&threadphoton, fcfg->f2);
	    DEBUG(("Photon #%d\t", thread));
	    PHOTON_DEBUG(&threadphoton);
	    /** Detect photon*/
	    DEBUG(("Photon #%d\t", thread));
	    res= detectPhoton(&threadphoton, cameraSensor);
	    atomicAdd(gfocused+fcfg->focusZindex, (unsigned int)res);
    }
}

__global__ void addingKernel(float * imageVolume, unsigned int zPlanes) {
	int thread = blockIdx.x * blockDim.x + threadIdx.x;
	if (thread < fcfg->sensorSize.x * fcfg->sensorSize.y) {
		for (int i = 1; i < zPlanes; i++) 
			imageVolume[thread] = imageVolume[thread] + imageVolume[thread + fcfg->sensorSize.x * fcfg->sensorSize.y * i];
	}
}

void volumeAdd(SPIMConfig*cfg, SPIMGPUInfo* gpu, float* imageVolume) {
	int blockSize = 1024;   /**< The launch configurator returned block size*/
	int gridSize = (gpu->autothread / gpu->autoblock);    /**< The actual grid size needed, based on input size*/
	gridSize = (cfg->sensorSize.x * cfg->sensorSize.y) / blockSize+1;
	printf("Generating OPT projection...\n");
	addingKernel << <gridSize, blockSize >> > (imageVolume,  cfg->nzPlanes);
	cudaDeviceSynchronize();
}


__global__ void output_addingKernel(float* outputVolume, float* imageVolume,int npixels) {
    int thread = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread < npixels) {
	outputVolume[thread] = outputVolume[thread] + imageVolume[thread];
    }
}


/**
*@brief Initializes focus config structure from MCX parameters
*/
void update_outputVolume(SPIMConfig* cfg, SPIMGPUInfo* gpu, float** outputVolume, float* imagevolume) {
    int blockSize = 1024;   /**< The launch configurator returned block size*/
    int gridSize = (cfg->sensorSize.x * cfg->sensorSize.y * cfg->nzPlanes) / blockSize + 1;
    float* goutputVolume;
    float* gimageVolume;
    printf("Updating output volume..\n");
    CUDA_LOAD(&goutputVolume, *outputVolume,cfg->sensorSize.x * cfg->sensorSize.y * cfg->nzPlanes,float);
    CUDA_LOAD(&gimageVolume, imagevolume, cfg->sensorSize.x * cfg->sensorSize.y * cfg->nzPlanes,float);

    output_addingKernel << <gridSize, blockSize >> > (goutputVolume, gimageVolume, cfg->sensorSize.x * cfg->sensorSize.y * cfg->nzPlanes);
    CUDA_ASSERT(cudaDeviceSynchronize());
    CUDA_UNLOAD(outputVolume, goutputVolume, cfg->sensorSize.x * cfg->sensorSize.y * cfg->nzPlanes,float);
    CUDA_ASSERT(cudaFree(goutputVolume));
    CUDA_ASSERT(cudaFree(gimageVolume));

}





__global__ void illTest() {
	printf("This is a test of illumination volume\n");
	for (int i = 0; i < 10; i++) {
		printf("\t Pos %d: %f\n", i, fcfg->illVolume[i]);
	}
}

/**
*@brief check card memory and estimate tunning configuration 
*/
 void memoryCalculator(focusConfig* hfcfg, int *nChunks, unsigned int *photonsPerChunk, unsigned int *remaindingPhotons) {
	size_t cudaMemory;			/**  Cuda total device memory*/
	size_t cudaFreeMemory;		/**  Cuda device free memory*/
	size_t cudaFreememory_ph;	/**  Cuda device free memory in photons*/
	

	/*Check  memory needs*/
	CUDA_ASSERT(cudaMemGetInfo(&cudaFreeMemory, &cudaMemory));
	printf("Photons struct size in GPU = %zu KB = %.2f GB\n", hfcfg->Nphotons * sizeof(photon), ((float)hfcfg->Nphotons * sizeof(photon)) / 1073741824.0f);
	printf("Free memory space in GPU = %d KB = %.2f GB of %.2f GB\n", cudaFreeMemory, cudaFreeMemory / 1073741824.0f, cudaMemory / 1073741824.0f);
	cudaFreememory_ph = (cudaFreeMemory / sizeof(photon)) * 0.8;

	/*Estimate memory needs*/
	*nChunks = (int)hfcfg->Nphotons / cudaFreememory_ph;
	*remaindingPhotons = (unsigned int)hfcfg->Nphotons % cudaFreememory_ph;
	if (remaindingPhotons > 0) ++(*nChunks);
	*photonsPerChunk = cudaFreememory_ph;

	printf("Processing %d photons in %u steps of %d photons. Remainding photons = %u\n",hfcfg->Nphotons,*nChunks,*photonsPerChunk, *remaindingPhotons);
	
}
 /**
*@brief Print CPU focusing Kernel
*/
 void print_focusConfig(SPIMConfig* cfg) {
	 printf("Focus Configuration:\n\tf1 = %f [mm]\n\tf2 = %f [mm]\n", cfg->f1, cfg->f2);
	 printf("\tmua = %f [cm^-1]\n\t Number of .mch files = %z\n", cfg->mua, cfg->nFiles);
	 return;
 }


int focus_launcher(SPIMConfig* cfg, SPIMGPUInfo* gpu, float** imageVolume) {
	focusConfig hfcfg;
	photon* photons;
	photon0Pos* photonsOrigin;
	float* photonsData = NULL;
	int nChunks;
	unsigned int photonsPerChunk, remainingPhotons;
	size_t cudaMemory;
	size_t cudaFreeMemory;
	float* gimageVolume;
	float zFocus;
	/* Prepare GPU settings*/
	int blockSize = gpu->autoblock;   /**< The launch configurator returned block size*/
	int gridSize = gpu->autothread / blockSize;    /**< The actual grid size needed, based on input size*/
	/* Prepare variables to import data files*/
	//int fidx = 0;
	char tmp[200];
	const char* dataFile;
	int Nphotons;

	CUDA_ASSERT(cudaMemGetInfo(&cudaFreeMemory, &cudaMemory));
	printf("Free memory space in GPU = %d KB = %.2f GB of %.2f GB\n", cudaFreeMemory, cudaFreeMemory / 1073741824.0f, cudaMemory / 1073741824.0f);

	/*Prepare camera buffer*/
	build_cameraSensor(cfg, &gimageVolume);

	/*Generate kernel focusing settings*/
	createfocusConfig(cfg, gpu, &hfcfg);
	hfcfg.photonsPerThread = (unsigned long int)ceil((float)LOAD_BUFFER_SIZE / (float)gpu->autothread);
	CUDA_ASSERT(cudaMemcpyToSymbol(fcfg, &hfcfg, sizeof(focusConfig)));

	/*print_focusConfig << <1, 1 >> > ();
	cudaDeviceSynchronize();*/
	   
	int photonOffset;

	for (int fidx = 0; fidx < cfg->nFiles; fidx++) {
		printf("Processing file %d of %d \n", fidx + 1, cfg->nFiles);
		/*Load photons data from file*/
		buildMchPath(cfg, tmp, fidx);
		dataFile = &(*tmp);
		Nphotons = readMchFile(dataFile, cfg, &photonsData);	
		//Nphotons = 100;
		if (Nphotons > 0) {
			printf("Data was loaded succesfully!!\n");
			cfg->Nphotons += (unsigned int)Nphotons;

			/*Update focusing kernel structure*/
			hfcfg.Nphotons = (unsigned long int)Nphotons;
			hfcfg.mcx_pxSize = cfg->mcx_pxSize;
			hfcfg.photonSize = cfg->photonSize;					
			SET_PHOTONS_DATA(hfcfg.Nphotons, hfcfg.mcx_pxSize, hfcfg.photonSize);

			/*Calculate memory settings*/
			memoryCalculator(&hfcfg, &nChunks, &photonsPerChunk, &remainingPhotons);
	
			/* Process photons struct in chunks if necessary*/
			photonOffset = 0;
			for (int chunk = 0; chunk < nChunks; chunk++) {
				photonOffset = photonsPerChunk * chunk * cfg->photonSize;
				if (chunk == nChunks - 1)
					load_photonStruct(&hfcfg, gpu, remainingPhotons, (photonsData + photonOffset), &photons);
				else
					load_photonStruct(&hfcfg, gpu, photonsPerChunk, (photonsData + photonOffset), &photons);
				printf("Running focusing kernel\n");
				for (int z = 0; z < cfg->nzPlanes; z++) {
					zFocus = ((cfg->zScan.y - cfg->zScan.x) / cfg->nzPlanes) * z + cfg->zScan.x;
					FOCUS_AT(zFocus / cfg->mcx_pxSize, z);
					//printf("Focus distance in grid units = %f\n", zFocus / hfcfg.mcx_pxSize);
					printf("\tFocusing at plane %d of %d. z= %f mm. \n", z + 1, cfg->nzPlanes, zFocus);
					focus4F << <gridSize, blockSize >> > (photons, gimageVolume);
					cudaDeviceSynchronize();
				}

				CUDA_ASSERT(cudaFree(photons));
			}
			free(photonsData);
			photonsData = NULL;
		}
		else
			printf("Error loading data\n");

	}
	printf("Focusing step finished\n");
	printf("Bringing images back to the host memory...\n");
	//build_photonStruct(cfg->Nphotons, photonsData, &photons, &photonsOrigin);	

	//printf("Running kernel for %d photons. %d blocks of %d threads per block\n", cfg->Nphotons, gridSize, blockSize);
	if (cfg->spimVol == 0) {
		volumeAdd(cfg, gpu, gimageVolume);		
		*imageVolume = (float*)calloc(cfg->sensorSize.x * cfg->sensorSize.y * 1, sizeof(float));
		CUDA_UNLOAD(imageVolume, gimageVolume, cfg->sensorSize.x * cfg->sensorSize.y * 1,float);
	}
	else {	
		*imageVolume = (float*)calloc(cfg->sensorSize.x * cfg->sensorSize.y * cfg->nzPlanes, sizeof(float));
		CUDA_UNLOAD(imageVolume, gimageVolume, cfg->sensorSize.x * cfg->sensorSize.y * cfg->nzPlanes,float);
	}
	
	printf("NA was %f\n", hfcfg.lens1R / (hfcfg.f1));

	//CUDA_ASSERT(cudaFree(photons));
	return 0;
}

int focus_mcx(SPIMConfig* cfg, SPIMGPUInfo* gpu, photon* photons, float** imageVolume) {
    focusConfig hfcfg;
    //photon* photons;
    //photon0Pos* photonsOrigin;
    //float* photonsData = NULL;

    int nChunks;
    unsigned int photonsPerChunk, remainingPhotons;
    size_t cudaMemory;
    size_t cudaFreeMemory;
    float* gimageVolume;
    float zFocus;
    /* Prepare GPU settings*/
    int blockSize = 1024;				/**< The launch configurator returned block size*/
    int gridSize = (cfg->Nphotons) / blockSize + 1;	/**< The actual grid size needed, based on input size*/
    /* Prepare variables to import data files*/
    //int fidx = 0;
    char tmp[200];
    const char* dataFile;
    int Nphotons;
    unsigned int *gfocused;
    int focused;

    CUDA_ASSERT(cudaMemGetInfo(&cudaFreeMemory, &cudaMemory));
    printf("Free memory space in GPU = %u KB = %.2f GB of %.2f GB\n", cudaFreeMemory, cudaFreeMemory / 1073741824.0f, cudaMemory / 1073741824.0f);

    /*Prepare camera buffer(instead of bulding a new sensor, we load buffer from host into GPU)*/
    CUDA_LOAD(&gimageVolume, *imageVolume, cfg->sensorSize.x * cfg->sensorSize.y * cfg->sensorSize.z,float);

   // build_cameraSensor(cfg, &gimageVolume);

    /*Generate kernel focusing settings*/
    createfocusConfig(cfg, gpu, &hfcfg);
    CUDA_ASSERT(cudaMemcpyToSymbol(fcfg, &hfcfg, sizeof(focusConfig)));

    /*Initialize photon stats*/
    CUDA_LOAD(&gfocused, cfg->focusedPhotons, cfg->sensorSize.z, unsigned int);

    //CUDA_ASSERT(cudaMalloc((void**)& gfocused, sizeof(int)));
    //CUDA_ASSERT(cudaMemset(gfocused, 0, sizeof(int)));
	  

    /*Update focusing kernel structure*/
    hfcfg.Nphotons = (unsigned long int)cfg->Nphotons;
    hfcfg.mcx_pxSize = cfg->mcx_pxSize;
    hfcfg.photonSize = cfg->photonSize;
    SET_PHOTONS_DATA(hfcfg.Nphotons, hfcfg.mcx_pxSize, hfcfg.photonSize);

    //print_photonData << < 1,1 >> > (photons, 30);
    prepare_MCXphotonData << <gridSize, blockSize >> > (photons);
    cudaDeviceSynchronize();
    //print_photonData << < 1, 1 >> > (photons, 30);


    /*Run focusing kernel for each zPlane*/
    for (int z = 0; z < cfg->nzPlanes; z++) {
	zFocus = ((cfg->zScan.y - cfg->zScan.x) / cfg->nzPlanes) * z + cfg->zScan.x;
	FOCUS_AT(zFocus / cfg->mcx_pxSize, z * cfg->spimVol);/** SpimVol=1 for SPIM, 0 for OPT*/
	//printf("Focus distance in grid units = %f\n", zFocus / hfcfg.mcx_pxSize);
	printf("\tFocusing at plane %d of %d. z= %f mm.\n", z + 1, cfg->nzPlanes, zFocus);
	focus4F_auto << <gridSize, blockSize >> > (photons, gimageVolume,gfocused);
	CUDA_ERROR_CHECK(cudaDeviceSynchronize());
	CUDA_UNLOAD(&(cfg->focusedPhotons), gfocused, cfg->sensorSize.z, unsigned int);
	//CUDA_ASSERT(cudaMemcpy(&focused, gfocused, sizeof(int), cudaMemcpyDeviceToHost));
	//printf("%e of %e photons used\n",(float)focused,(float)cfg->Nphotons);

    }
    
    printf("Focusing step finished\n");
    printf("Bringing images back to the host memory...\n");
    //build_photonStruct(cfg->Nphotons, photonsData, &photons, &photonsOrigin);	
    CUDA_UNLOAD(&(cfg->focusedPhotons),gfocused, cfg->sensorSize.z, unsigned int);

    //printf("Running kernel for %d photons. %d blocks of %d threads per block\n", cfg->Nphotons, gridSize, blockSize);
    CUDA_UNLOAD(imageVolume, gimageVolume, cfg->sensorSize.x * cfg->sensorSize.y * cfg->sensorSize.z,float);
    /*if (cfg->spimVol == 0) {
	//volumeAdd(cfg, gpu, gimageVolume);
	//*imageVolume = (float*)calloc(cfg->sensorSize.x * cfg->sensorSize.y * 1, sizeof(float));
	CUDA_UNLOAD(imageVolume, gimageVolume, cfg->sensorSize.x * cfg->sensorSize.y * 1);
    }
    else {
	//*imageVolume = (float*)calloc(cfg->sensorSize.x * cfg->sensorSize.y * cfg->nzPlanes, sizeof(float));
	CUDA_UNLOAD(imageVolume, gimageVolume, cfg->sensorSize.x * cfg->sensorSize.y * cfg->nzPlanes);
    }*/
    CUDA_ASSERT(cudaFree(gfocused));
    CUDA_ASSERT(cudaFree(gimageVolume));
    CUDA_ASSERT(cudaFree(fcfg->illVolume));

    return 1;
}

/**
 * @brief Utility function to calculate the GPU stream processors (cores) per SM
 *
 * Obtain GPU core number per MP, this replaces
 * ConvertSMVer2Cores() in libcudautils to avoid
 * extra dependency.
 *
 * @param[in] v1: the major version of an NVIDIA GPU
 * @param[in] v2: the minor version of an NVIDIA GPU
 */
int corecount(int v1, int v2) {
	int v = v1 * 10 + v2;
	if (v < 20)      return 8;
	else if (v < 21) return 32;
	else if (v < 30) return 48;
	else if (v < 50) return 192;
	else if (v < 60) return 128;
	else          return 64;
}

/**
 * @brief Utility function to calculate the maximum blocks per SMX
 *
 * @param[in] v1: the major version of an NVIDIA GPU
 * @param[in] v2: the minor version of an NVIDIA GPU
 */
int smxblock(int v1, int v2) {
	int v = v1 * 10 + v2;
	if (v < 30)      return 8;
	else if (v < 50) return 16;
	else          return 32;
}

/**
 * @brief Utility function to query GPU info and set active GPU
 *
 * This function query and list all available GPUs on the system and print
 * their parameters. This is used when -L or -I is used.
 *
 * @param[in,out] cfg: the simulation configuration structure
 * @param[out] info: the GPU information structure
 */

int list_gpu(SPIMConfig*cfg, SPIMGPUInfo**info) {


	int dev;
	int deviceCount, activedev = 0;

	cudaGetDeviceCount(&deviceCount);
	if (deviceCount == 0) {
		PRINT(stderr, "No CUDA-capable GPU device found\n");
		return 0;
	}
	*info = (SPIMGPUInfo*)calloc(deviceCount, sizeof(SPIMGPUInfo));
	if (cfg->gpuid && cfg->gpuid > deviceCount) {
		PRINT(stderr, "Specified GPU ID is out of range\n");
		return 0;
	}
	// scan from the first device
	for (dev = 0; dev < deviceCount; dev++) {
		cudaDeviceProp dp;
		CUDA_ASSERT(cudaGetDeviceProperties(&dp, dev));

		if (cfg->isgpuinfo == 3)
			activedev++;
		else if (cfg->deviceid[dev] == '1') {
			cfg->deviceid[dev] = '\0';
			cfg->deviceid[activedev] = dev + 1;
			activedev++;
		}
		strncpy((*info)[dev].name, dp.name, MAX_SESSION_LENGTH);
		(*info)[dev].id = dev + 1;
		(*info)[dev].devcount = deviceCount;
		(*info)[dev].major = dp.major;
		(*info)[dev].minor = dp.minor;
		(*info)[dev].globalmem = dp.totalGlobalMem;
		(*info)[dev].constmem = dp.totalConstMem;
		(*info)[dev].sharedmem = dp.sharedMemPerBlock;
		(*info)[dev].regcount = dp.regsPerBlock;
		(*info)[dev].clock = dp.clockRate;
		(*info)[dev].sm = dp.multiProcessorCount;
		(*info)[dev].core = dp.multiProcessorCount*corecount(dp.major, dp.minor);
		(*info)[dev].maxmpthread = dp.maxThreadsPerMultiProcessor;
		//(*info)[dev].maxgate = cfg->maxgate;
		(*info)[dev].autoblock = (*info)[dev].maxmpthread / smxblock(dp.major, dp.minor);
		(*info)[dev].autothread = (*info)[dev].autoblock * smxblock(dp.major, dp.minor) * (*info)[dev].sm;

		if (strncmp(dp.name, "Device Emulation", 16)) {
			if (cfg->isgpuinfo) {
				PRINT(stdout, "=============================   GPU Infomation  ================================\n");
				PRINT(stdout, "Device %d of %d:\t\t%s\n", (*info)[dev].id, (*info)[dev].devcount, (*info)[dev].name);
				PRINT(stdout, "Compute Capability:\t%u.%u\n", (*info)[dev].major, (*info)[dev].minor);
				PRINT(stdout, "Global Memory:\t\t%u B\nConstant Memory:\t%u B\n"
					"Shared Memory:\t\t%u B\nRegisters:\t\t%u\nClock Speed:\t\t%.2f GHz\n",
					(unsigned int)(*info)[dev].globalmem, (unsigned int)(*info)[dev].constmem,
					(unsigned int)(*info)[dev].sharedmem, (unsigned int)(*info)[dev].regcount, (*info)[dev].clock*1e-6f);
#if CUDART_VERSION >= 2000
				PRINT(stdout, "Number of MPs:\t\t%u\nNumber of Cores:\t%u\n",
					(*info)[dev].sm, (*info)[dev].core);
				PRINT(stdout, "Optimized blocks to launch:\t%u\nOptimized threads to launch:\t%u\n",
					(*info)[dev].autoblock, (*info)[dev].autothread);
#endif
				PRINT(stdout, "SMX count:\t\t%u\n", (*info)[dev].sm);
				PRINT(stdout, "================================================================================\n\n");
			}
		}
	}
	//if (cfg->isgpuinfo == 2 && cfg->parentid == mpStandalone) { //list GPU info only
	//	exit(0);
	//}
	if (activedev < MAX_DEVICE)
		cfg->deviceid[activedev] = '\0';
	return activedev;
}
