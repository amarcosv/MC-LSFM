#include "MCSPIM_kernel.h"
#include "mcx_utils.h"
#include "mcx_core.h"
#include "focus_utils.h"
#include "focusKernel.h"
#include "string.h"
#include <time.h>
//--session line --root D:\Results\MonteCarlo_LSFM\Spiral\test_set\line -f D:\Results\MonteCarlo_LSFM\Spiral\test_set\line\vesselTest.json --outputformat mc2 --gpu 1 --photon 1470000000 --normalize 1 --save2pt 1 --reflect 0 --savedet 1 --unitinmm 0.10 --srcfrom0 1 --seed 1648335518 --saveseed 0 --specular 0 --array 0 --dumpmask 0 --repeat 1 -w XVWL --maxdetphoton 100000000 --bc aaaaaa
#define MU "\u03BC"

#define CUDA_ERROR_CHECK

#define CUDA_ASSERT( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
inline void __cudaSafeCall(cudaError err, const char* file, const int line)
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
inline void __cudaCheckError(const char* file, const int line)
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
*@brief prints photon position information
*/
__device__ void print_photon_id(photon* photons,uint id) {
    printf("Printing photon #%d\n", id);
    printf("\t pos: [%f, %f, %f]\n", photons[id].pos.x, photons[id].pos.y, photons[id].pos.z);
    printf("\t dir: [%f, %f, %f]\n", photons[id].dir.x, photons[id].dir.y, photons[id].dir.z);
    printf("\t pos0: [%f, %f, %f]\n", photons[id].Opos.x, photons[id].Opos.y, photons[id].Opos.z);
    printf("\t w:	 %e\n", photons[id].w);

}

/**
*@brief prints photon position information
*/
__global__ void print_photon_id_cpu(photon* photons, uint id) {
    printf("Printing photon #%d\n", id);
    printf("\t pos: [%f, %f, %f]\n", photons[id].pos.x, photons[id].pos.y, photons[id].pos.z);
    printf("\t dir: [%f, %f, %f]\n", photons[id].dir.x, photons[id].dir.y, photons[id].dir.z);
    printf("\t pos0: [%f, %f, %f]\n", photons[id].Opos.x, photons[id].Opos.y, photons[id].Opos.z);
    printf("\t w:	 %e\n", photons[id].w);

}

/**
*@brief prints photon position information
*/
__device__ void print_photonPosition(float* photons, int index) {
    printf("Photon %d: pos[%f]\n", index, photons[index]);
    printf("Photon %d: pos[%f]\n", index+1, photons[index+1]);
    printf("Photon %d: pos[%f]\n", index+2, photons[index+2]);
    printf("Photon %d: pos[%f]\n", index+3, photons[index+3]);
    printf("Photon %d: pos[%f]\n", index+4, photons[index+4]);
    printf("Photon %d: pos[%f]\n", index+5, photons[index+5]);
    printf("Photon %d: pos[%f]\n", index+7, photons[index+7]);
    printf("Photon %d: pos[%f]\n", index+8, photons[index+8]);
    printf("Photon %d: pos[%f]\n", index+9, photons[index+9]);
    printf("Photon %d: pos[%f]\n", index + 6, photons[index + 6]);
}

/**
*@brief prints photon position information
*/
__global__ void print_photonPosition_cpu(float* photons, int index) {
    printf("Photon %d: pos[%f]\n", index, photons[index]);
    printf("Photon %d: pos[%f]\n", index + 1, photons[index + 1]);
    printf("Photon %d: pos[%f]\n", index + 2, photons[index + 2]);
    printf("Photon %d: pos[%f]\n", index + 3, photons[index + 3]);
    printf("Photon %d: pos[%f]\n", index + 4, photons[index + 4]);
    printf("Photon %d: pos[%f]\n", index + 5, photons[index + 5]);
    printf("Photon %d: pos[%f]\n", index + 7, photons[index + 7]);
    printf("Photon %d: pos[%f]\n", index + 8, photons[index + 8]);
    printf("Photon %d: pos[%f]\n", index + 9, photons[index + 9]);
    printf("Photon %d: pos[%f]\n", index + 6, photons[index + 6]);
}

/**
*@brief Initializes focus config structure from MCX parameters
*/
void merge_ConfigFiles(MCXConfig* mcxcfg, SPIMConfig* spimcfg) {
    spimcfg->workingDir = (char*)malloc(strlen(mcxcfg->rootpath));
    strcpy(spimcfg->workingDir, mcxcfg->rootpath);
    spimcfg->fileBaseName = (char*)malloc(strlen(mcxcfg->session));
    strcpy(spimcfg->fileBaseName, mcxcfg->session);
    spimcfg->mcx_pxSize = mcxcfg->unitinmm;
    spimcfg->mcx_simVolume = mcxcfg->dim;
    spimcfg->mua = mcxcfg->prop[1].mua; /** Comes in grid units*/
}

void initstats(MCSPIMStats* stats) {
    stats->simPhotons = 0;
    stats->detPhotons = 0;
    stats->mcxTime = 0;
    stats->focusTime = 0;
    stats->simTime = 0;
 }

void print_runningstats(MCSPIMStats* stats, SPIMConfig* spimconfig) {

    printf("**********************************************************************************\n\n");
    for (int i = 0; i < spimconfig->sensorSize.z; i++)
	printf("\t z Plane %d of %d. z= %f mm . Focused photons was %e\n", i + 1, spimconfig->sensorSize.z, ((spimconfig->zScan.y - spimconfig->zScan.x) / spimconfig->nzPlanes) * i + spimconfig->zScan.x, (float)stats->focusedPhotons[i]);
    printf("Total simulated photons was %e\n", stats->simPhotons);
    printf("Total detected photons was %e\n", stats->detPhotons);
    printf("Execution time: \n");
    printf("\tTotal: %f secs = %d hh %d mm %d ss\n", stats->simTime, (int)stats->simTime / 3600, ((int)stats->simTime % 3600) / 60, ((int)stats->simTime % 3600) % 60);
    printf("\tMCX:   %f secs \n", stats->mcxTime);
    printf("\tFocus: %f secs \n", stats->focusTime);
    printf("**********************************************************************************\n\n");
}

int save_runingLog(MCSPIMStats* stats, SPIMConfig* spimconfig, MCXConfig* mcxconfig) {
	FILE* fileID;
	size_t written;
	char idx[5];
	int fidx = 0;
	const char* dataFile;	
	char tmp[200];
	strcpy(tmp, spimconfig->workingDir);
	//strcat(tmp, "\\");
	strcat(tmp, spimconfig->outputFileName);

	strcat(tmp, ".txt");
	dataFile = &(*tmp);

	/** Get current time and date*/
	time_t t;   
	time(&t);

	printf("Writting results in file: %s\n", tmp);

	fileID = fopen(dataFile, "wb");
	if (NULL == fileID)
	    return -1; 
	fprintf(fileID,"***********************\tMC SPIM Simulation\t***********************\n");
	fprintf(fileID, "%s\n", ctime(&t));
	fprintf(fileID, "***********************\tData files and directories\t***********************\n");
	fprintf(fileID, "Output file\n");
	fprintf(fileID, "\tName = %s.mcspim\n", spimconfig->outputFileName + 1);
	fprintf(fileID, "\tDimensions [x y z] = [%u %u %u] px \n", spimconfig->sensorSize.x, spimconfig->sensorSize.x, spimconfig->sensorSize.z);
	if (spimconfig->spimVol == 1)
	    fprintf(fileID, "\tVoxel size = [%f %f %f] mm\n", (spimconfig->f2 / spimconfig->f1) * spimconfig->sensorpxSize, (spimconfig->f2 / spimconfig->f1) * spimconfig->sensorpxSize, (spimconfig->zScan.y - spimconfig->zScan.x) / (float)spimconfig->nzPlanes);
	else
	    fprintf(fileID, "\Pixel size = %f mm\n", (spimconfig->f2 / spimconfig->f1)*spimconfig->sensorpxSize); 	
	fprintf(fileID, "Data directories\n");
	fprintf(fileID, "\tBasename = %s\n", spimconfig->fileBaseName);
	fprintf(fileID, "\tWorking directory = %s\n", spimconfig->workingDir);
	fprintf(fileID, "\tConfig file = %s\n", spimconfig->configFileName+1);
	fprintf(fileID, "\tIllumination file name = %s\n", spimconfig->illFileName+1);

	fprintf(fileID, "\n***********************\tSimulation parameters\t***********************\n");
	fprintf(fileID, "Number of simulations = %u\n", spimconfig->nFiles);
	fprintf(fileID, "Optical properties:\n");
	fprintf(fileID, "\t%sa = %f 1/cm\n", MU, 10.f * mcxconfig->prop[1].mua / mcxconfig->unitinmm);
	fprintf(fileID, "\t%ss = %f 1/cm\n", MU, 10.f * mcxconfig->prop[1].mus / mcxconfig->unitinmm);
	fprintf(fileID, "\tg = %f\n", mcxconfig->prop[1].g);
	fprintf(fileID, "\tn = %f\n", mcxconfig->prop[1].n);
	fprintf(fileID, "Fluorescence:\n");
	fprintf(fileID, "\tSource type = %d\n", mcxconfig->srctype);
	if(mcxconfig->srctype==16)
	fprintf(fileID, "\tSource file = %s\n", mcxconfig->fluoname);
	fprintf(fileID, "\tNumber of fluorophores = %u\n",(unsigned int) mcxconfig->srcparam1.x);

	fprintf(fileID, "MCX configuration:\n");
	fprintf(fileID, "\tVolume size = [%u %u %u] px\n", spimconfig->mcx_simVolume.x, spimconfig->mcx_simVolume.y, spimconfig->mcx_simVolume.z);
	fprintf(fileID, "\tPixel size = %f\n", spimconfig->mcx_pxSize);
	fprintf(fileID, "Illumination configuration:\n");
	fprintf(fileID, "\tVolume size = [%u %u %u] px\n", spimconfig->ill_simVolume.x, spimconfig->ill_simVolume.y, spimconfig->ill_simVolume.z);
	fprintf(fileID, "\tPixel size = %f\n", spimconfig->ill_pxSize);
	fprintf(fileID, "\tDetailed illumination = %s\n", (spimconfig->ill_simVolume.y == 1) ? "Yes" : "No");
		
	fprintf(fileID, "Camera configuration:\n");
	fprintf(fileID, "\tSensor size  = [%u %u] px = [%f %f] mm\n", spimconfig->sensorSize.x, spimconfig->sensorSize.y, ((float)spimconfig->sensorSize.x) * spimconfig->sensorpxSize, ((float)spimconfig->sensorSize.y) * spimconfig->sensorpxSize);
	fprintf(fileID, "\tSensor pixel size = %f mm\n", spimconfig->sensorpxSize);
	fprintf(fileID, "\tFocal length [f1 f2] = [%.2f %.2f] mm \n", spimconfig->f1, spimconfig->f2);
	fprintf(fileID, "\tNumerical aperture [NA1 NA2] = [%.2f %.2f] \n", spimconfig->NA1, spimconfig->NA2);	
	fprintf(fileID, "SPIM Acquisition:\n");
	fprintf(fileID, "\tZ scanning range [init end zstep] = [%.2f %.2f %.4f] mm\n", spimconfig->zScan.x, spimconfig->zScan.y, (spimconfig->zScan.y- spimconfig->zScan.x)/(float)spimconfig->nzPlanes);
	fprintf(fileID, "\tZ planes  = %u\n", spimconfig->nzPlanes);
	fprintf(fileID, "\tImage mode  = %s\n", ((spimconfig->spimVol==1) ? "SPIM" : "OPT"));

	fprintf(fileID, "\n***********************\tRunning stats\t***********************\n");
	fprintf(fileID, "Simulated photons: \n");
	fprintf(fileID, "\tTotal = %e\n", stats->simPhotons);
	fprintf(fileID, "\tDetected = %e\n", stats->detPhotons);
	fprintf(fileID, "\tDetection rate = %.2f%%\n", stats->detPhotons/stats->simPhotons*100);
	fprintf(fileID, "Execution time: \n");
	fprintf(fileID, "\tTotal: %f secs = %d hours %d mins %d secs\n", stats->simTime, (int)stats->simTime / 3600, ((int)stats->simTime % 3600) / 60, ((int)stats->simTime % 3600) % 60);
	fprintf(fileID, "\tMCX:   %f secs \n", stats->mcxTime);
	fprintf(fileID, "\tFocus: %f secs \n", stats->focusTime);
	fprintf(fileID, "Focused photons per plane : \n");
	for (int i = 0; i < spimconfig->sensorSize.z; i++)
	    fprintf(fileID, "\t z Plane %d of %d. z= %f mm . Focused photons was %e\n", i + 1, spimconfig->sensorSize.z, ((spimconfig->zScan.y - spimconfig->zScan.x) / spimconfig->nzPlanes) * i + spimconfig->zScan.x, (float)stats->focusedPhotons[i]);

	fclose(fileID);
	printf("File written succesfully!\n");
	return 0;

}

/**
*@brief Launches mcx simulation and MCfocusing code
*/
int simulationLauncher(int argc, char* argv[]) {

    /*! structure to store all simulation parameters
    */
    SPIMConfig spimconfig;	    /**< spimconfig: structure to store all focusing algorithm parameters */
    MCXConfig  mcxconfig;           /**< mcxconfig: structure to store all simulation parameters */
    GPUInfo* gpuinfo = NULL;        /**< gpuinfo: structure to store GPU information */
    SPIMGPUInfo* spimgpuinfo = NULL;    /**< Structure for GPU configuration */
    MCSPIMStats stats;		    /**< stats: structure with all execution stats*/

    float* outputVolume;	    /**< outputVolume: buffer to store output volume */
    float* imageVolume;		    /**< imageVolume: buffer to store spim volume */
    float* illVolume = NULL;	    /**< illVolume: buffer to store illuminationVolume */
    unsigned int activedev = 0;     /**< activedev: count of total active GPUs to be used */
    
    size_t cudaMemory;			/**<  Cuda total device memory*/
    size_t cudaFreeMemory;		/**<  Cuda device free memory*/
    size_t imageVolumeSize;		/**<  Output volume size*/
    size_t illVolumeSize;		/**<  Illumination volume size*/
    size_t cudaFreememory_ph;	/**<  Cuda device free memory in photons*/
  
    size_t Nphotons;			/**<  Photons to run at each simulation*/
    float* gPdet;			/**< GPU pointer to MCX detected photons array*/
    photon* detphotons;			/**< GPU pointer to casted detected photons */

    /** Variables for running statistics*/
    float simPhotons = 0;		/**< Total simulated photons*/
    float detPhotons = 0;		/**< Total detected photons*/   
    float mcxTime = 0;			/**< Total MCX simulation time*/   
    float focusTime = 0;		/**< Total focusing time*/  
    float simTime = 0;			/**< Total simulation time*/
    clock_t begin, end;
    clock_t simStart, simEnd;

    simStart = clock();

    /**
       To start an MCX simulation, we first create a simulation configuration and
       set all elements to its default settings.
     */
    mcx_initcfg(&mcxconfig);
    initcfg(&spimconfig);
    initstats(&stats);

    /**
       Then, we parse the full command line parameters and set user specified settings
     */
    mcx_parsecmd(argc, argv, &mcxconfig);
    merge_ConfigFiles(&mcxconfig, &spimconfig);
    readConfigFile(&spimconfig);

    /** The next step, we identify gpu number and query all GPU info */
    if (!(activedev = mcx_list_gpu(&mcxconfig, &gpuinfo))) {
	mcx_error(-1, "No GPU device found\n", __FILE__, __LINE__);
    }
    if (!(list_gpu(&spimconfig, &spimgpuinfo))) {
	printf("no active GPU device found\n");
    }
    /**
	Specify here compulsary settings for a MCX dedicated MCSPIM simulation
    */
    Nphotons = mcxconfig.nphoton;
    mcxconfig.savedetflag = 180;
    //mcxconfig.nphoton= 1000;
    /**
	Prepare output data here
    */
    outputVolume = (float *)calloc(spimconfig.sensorSize.x * spimconfig.sensorSize.y * spimconfig.nzPlanes,sizeof(float));
    spimconfig.focusedPhotons = (unsigned int*)calloc(spimconfig.sensorSize.z, sizeof(unsigned int));
    stats.focusedPhotons = spimconfig.focusedPhotons;

    /**Perform memory precalculations to ensure enough memory  */
    CUDA_ASSERT(cudaMemGetInfo(&cudaFreeMemory, &cudaMemory));
    imageVolumeSize = spimconfig.sensorSize.x * spimconfig.sensorSize.y * spimconfig.sensorSize.z * sizeof(float);
    illVolumeSize = spimconfig.ill_simVolume.x * spimconfig.ill_simVolume.y * spimconfig.ill_simVolume.z*sizeof(float);
    cudaFreememory_ph = ((cudaFreeMemory - (imageVolumeSize + illVolumeSize)) / sizeof(photon)) * 0.97;

    imageVolume=(float*)calloc(spimconfig.sensorSize.x * spimconfig.sensorSize.y * spimconfig.sensorSize.z, sizeof(float));
    /**
       This line runs the main MCX simulation for each GPU inside each thread
     */
  
    
    printf("Loading illumination file in memory.......\n");

    //	printf("Error loading data\n");
    int result;
    result = importIlluminationData(&spimconfig, &illVolume);
    if (result == 0)
	printf("Illumination was loaded succesfully!!\n");
    else {
	printf("ERROR: Illumination file couldn't be loaded.\n");
	printf("Exiting the program....\n");
	return 1;
    }

    for (int i = 0; i < spimconfig.nFiles; i++) {
	printf("Running simulation %d of %d...\n", i + 1,spimconfig.nFiles);
	CUDA_ASSERT(cudaMemGetInfo(&cudaFreeMemory, &cudaMemory));
	printf("Initially free memory space in GPU = %u KB = %.2f GB of %.2f GB\n", cudaFreeMemory, cudaFreeMemory / 1073741824.0f, cudaMemory / 1073741824.0f);
	mcxconfig.nphoton= Nphotons ;
	mcxconfig.maxdetphoton = (unsigned int)cudaFreememory_ph;

	begin = clock();
	mcx_run_simulation(&mcxconfig, gpuinfo, &gPdet);
	end = clock();

	if (mcxconfig.detectedcount > mcxconfig.maxdetphoton)
	    spimconfig.Nphotons = mcxconfig.maxdetphoton;
	else
	    spimconfig.Nphotons = mcxconfig.detectedcount;
	
	stats.detPhotons = stats.detPhotons + (float)spimconfig.Nphotons;
	stats.simPhotons = stats.simPhotons + (float)mcxconfig.nphoton;
	stats.mcxTime = stats.mcxTime + (float)(end - begin) / CLOCKS_PER_SEC;

	detphotons = (photon*)gPdet;
	printf("MCX Simulation ended after running %e photons\n", (float)mcxconfig.nphoton);
	printf("\Detected %e photons of a maximum of %e\n", (float)spimconfig.Nphotons, (float)mcxconfig.maxdetphoton);
	printf("\nUsing %f %% of photon memory\n", ((float)spimconfig.Nphotons / (float)mcxconfig.maxdetphoton) * 100.0f);

	/*Prepare camera buffer*/
        //build_cameraSensor(cfg, &gimageVolume);

	begin = clock();
	focus_mcx(&spimconfig, spimgpuinfo, detphotons, &imageVolume);
	end = clock();

	stats.focusTime= stats.focusTime + (float)(end - begin) / CLOCKS_PER_SEC;

	detphotons = NULL;
	CUDA_ASSERT(cudaFree(gPdet));

	//update_outputVolume(&spimconfig, spimgpuinfo, &outputVolume, imageVolume); We reuse our buffer all the time

	exportVolume(imageVolume, &spimconfig);
	if (result == 0)
	    printf("Ouput was exported succesfully!!\n");
	else
	    printf("Error writting data file\n");
	/**
	   Once simulation is complete, we clean up the allocated memory in config and gpuinfo, and exit
	 */
	mcxconfig.seed = mcxconfig.seed + 1;
	CUDA_ASSERT(cudaDeviceReset());

    }
    simEnd = clock();
    stats.simTime = (float)(simEnd - simStart) / CLOCKS_PER_SEC;

    save_runingLog(&stats, &spimconfig, &mcxconfig);
    print_runningstats(&stats, &spimconfig);



    mcx_cleargpuinfo(&gpuinfo);
    mcx_clearcfg(&mcxconfig);
    return 0;


}