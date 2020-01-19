#include "focus_utils.h""
#include "focusKernel.h"

#include <stdlib.h>
#include <cuda.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <cstring>
#include <time.h>

#ifdef  MCML_MEX
#include <matrix.h>
#include <mex.h>
/**<  Macro to read the 1st scalar cfg member */
#define GET_INPUT_FIELD(cfg,field,type)  if(strcmp(name,#field)==0) {double *val=mxGetPr(item);cfg->field=(type)val[0];printf("mcml.%s=%g;\n",#field,(float)(cfg->field));}
#define GET_FLOAT_INPUT_FIELD(cfg,field)  if(strcmp(name,#field)==0) {mxSingle *val=mxGetSingles(item);cfg->field=(float)val[0];printf("mcml.%s=%g;\n",#field,(cfg->field));}
#define GET_UINT8_INPUT_FIELD(cfg,field)  if(strcmp(name,#field)==0) {mxUint8 *val=mxGetUint8s(item);cfg->field=(unsigned short)val[0];printf("mcml.%s=%g;\n",#field,(cfg->field));}
#define GET_UINT32_INPUT_FIELD(cfg,field)  if(strcmp(name,#field)==0) {mxUint32 *val=mxGetUint32s(item);cfg->field=(unsigned long)val[0];printf("mcml.%s=%g;\n",#field,(cfg->field));}
#define MCML_PRINT(x) mexPrintf x
#else
#define MCML_PRINT(x) printf x
#endif


using namespace std;




int main(int argc, char* argv[]) {
	SPIMConfig cfg;				/**< Structure for simulation configuration */
	SPIMGPUInfo* gpuinfo = NULL;/**< Structure for GPU configuration */
	float* photonsData = NULL;
	float* illVolume = NULL;
	float* imageVolume;
	const char* illuminationFile = "D:/Results/MonteCarlo_LSFM/Spiral/mus2_mua0025_spiral/illuminationVolume.ill";
	//const char* dataFile = "D:/Results/MonteCarlo_LSFM/Spiral/mus2_mua0025_spiral/spiral_out_1.mch";
	printf("Initializing GPU....\n");
		initcfg(&cfg);
	if (!(list_gpu(&cfg, &gpuinfo))) {
		printf("no active GPU device found\n");
	}
	printf("Command line input commands:\n");
	parseCommandLine(argc, argv, &cfg);

	if (cfg.configFileName != NULL) {
		printf("Reading configuration file...\n");
		readConfigFile(&cfg);
	}

	/**/
	//int result= importPhotonsData(dataFile, &cfg, &photonsData);
	printf("Loading illumination file in memory.......\n");

	//	printf("Error loading data\n");
	int result;
	result = importIlluminationData(&cfg, &illVolume);
	if (result == 0)
		printf("Illumination was loaded succesfully!!\n");
	else
		printf("Error loading illumination data\n");

	printf("Running focusing algorithm.......\n");
	result = focus_launcher(&cfg, gpuinfo, &imageVolume);
	if (result == 0)
		printf("Focus performed succesfully!!\n");
	else
		printf("Error during the kernel execution\n");

	//dataFile = "D:/Results/MonteCarlo_LSFM/Spiral/mus2_mua0025_spiral/spiral_out_1.mcspim";
	//dataFile= "C:/Users/Asier\ Marcos/Documents/GitHub/MonteCarloFocus/data/data2fluoro_out.mc";

	exportVolume(imageVolume, &cfg);
	if (result == 0)
		printf("Ouput was exported succesfully!!\n");
	else
		printf("Error writting data file\n");
		
	return 0;
}



