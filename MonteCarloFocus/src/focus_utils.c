#include "focusKernel.h"
#include "focus_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
* Code equivalence for config file and command line entry
*/
const char* configFields[] = { "WORKINGPATH", "BASENAME", "NFILES", "MCX_SIMVOLUME", "F1", "F2", "NA1", "NA2","SENSORSIZE", "SENSORPXSIZE", "ZSCAN", "SPIMVOLUME", "DETAILEDILL", "ILLNAME"};
const char* commandFields[] = { "-w", "-b", "-n", "-m" ,"-f1" ,"-f2" , "-na1", "-na2", "-s", "-p", "-z", "-v", "-d", "-i"};
//const char* configFields[] = { "WORKINGPATH", "BASENAME", "NFILES", "MCX_PXSIZE", "MCX_SIMVOLUME", "F1", "F2", "NA1", "NA2","SENSORSIZE", "SENSORPXSIZE", "ZSCAN", "SPIMVOLUME", "DETAILEDILL" };
//const char* commandFields[] = { "-d", "-b", "-n", "-mp", "-mv" ,"-f1" ,"-f2" , "-na1", "-na2", "-s", "-p", "-z", "-v", "-i" };


/**
*@brief Initialization of config structure with default values.
*/
void initcfg(SPIMConfig*cfg) {
	cfg->Nphotons = 100;
	cfg->f1 = 100.0f;
	cfg->f2 = 100.0f;
	cfg->NA1 = 0.9f;
	cfg->NA2 = 0.9f;

	cfg->mcx_pxSize = 0.025f;
	cfg->mcx_simVolume.x = 300;
	cfg->mcx_simVolume.y = 400;
	cfg->mcx_simVolume.z = 300;

	cfg->sensorpxSize = 0.025f;
	cfg->sensorSize.x = 200;
	cfg->sensorSize.y = 200;

	cfg->ill_simVolume = cfg->mcx_simVolume;
	cfg->ill_pxSize = cfg->mcx_pxSize;
	
	cfg->photonSize = 10;

	cfg->nzPlanes = 200;
	cfg->zScan.x= 0.0f;
	cfg->zScan.y = 5.0f;

	cfg->mua = 0.25f;

	cfg->configFileName = NULL;
	cfg->outputFileName = NULL;
	cfg->illFileName = "\\illuminationVolume.ill";

	cfg->spimVol = 1;
}

/**
*@brief Build path for .mch file given a file index
*/
void buildMchPath(SPIMConfig* cfg, char* tmp, int fidx) {
	//char tmp[200];
	char idx[3];

	strcpy(tmp, cfg->workingDir);
	strcat(tmp, "\\");
	strcat(tmp, cfg->fileBaseName);
	strcat(tmp, "_out_");
	itoa(fidx, idx, 10);
	strcat(tmp, idx);
	strcat(tmp, ".mch");
	//dataFile = &(*tmp);
}


/**
*@brief Parse command line arguments
*/

int parseInput(char* argv[], SPIMConfig* cfg) {
	char* strBuffer;
	int tempLength;
	switch (argv[0][1]) {
	case 'w':
		strBuffer= malloc(strlen(argv[1]));		
		strcpy(strBuffer, argv[1]);
		int i = 2;
		while (argv[i][0] != '-') {
			tempLength = strlen(strBuffer);
			strBuffer = realloc(strBuffer,strlen(strBuffer)+1+strlen(argv[i]));
			*(strBuffer + tempLength) = ' ';
			strcpy(strBuffer+tempLength+1, argv[i]);
			i++;
		}
		printf("\tWorking directory path is %s\n", strBuffer);
		cfg->workingDir = malloc(strlen(strBuffer));		
		strcpy(cfg->workingDir, strBuffer);
		return 1;

	case 'c':
		printf("\tConfiguration file name is %s\n", argv[1]);
		strBuffer = malloc(strlen(argv[1])+1);
		*strBuffer = '\\';
		strcpy(strBuffer+1, argv[1]);
		cfg->configFileName = malloc(strlen(strBuffer));
		strcpy(cfg->configFileName, strBuffer);
		return 1;

	case 'b':
		printf("\tFile basename is %s\n", argv[1]);
		cfg->fileBaseName = malloc(strlen(argv[1]));
		//cfg->fileBaseName = argv[1];
		strcpy(cfg->fileBaseName, argv[1]);
		return 1;
		
	case 'n':
		if (argv[0][2] == 'a') {
			if (argv[0][3] == '1') {
				printf("\tNumerical aperture lens 1 = %s\n", argv[1]);
				cfg->NA1 = atof(argv[1]);
			}
			if (argv[0][3] == '2') {
				printf("\tNumerical aperture lens 2 = %s\n", argv[1]);
				cfg->NA2 = atof(argv[1]);
			}

		}
		else {
			printf("\tNumber of .mch files = %s\n", argv[1]);
			cfg->nFiles = atoi(argv[1]);
		}
		return 1;	

	case 'f':
		if (argv[0][2] == '1') {
			printf("\tFocal distance f1= %s [mm]\n", argv[1]);
			cfg->f1 = atof(argv[1]);
		}
		else if (argv[0][2] == '2') {
			printf("\tFocal distance f2= %s [mm]\n", argv[1]);
			cfg->f2 = atof(argv[1]);
		}
		return 1;
		
	case 'm':
		//if (argv[0][2] == 'v') {
			cfg->mcx_simVolume.x = atoi(argv[1]);
			cfg->mcx_simVolume.y = atoi(argv[2]);
			cfg->mcx_simVolume.z = atoi(argv[3]);
			cfg->mcx_pxSize = atof(argv[4]);
			printf("\tMCX simulation Volume dimensions = [%s %s %s]\n", argv[1], argv[2], argv[3]);
			printf("\tMCX simulation pixel size = %s (mm)\n", argv[4]);
			if (cfg->ill_simVolume.y != 1) {
			    cfg->ill_pxSize = cfg->mcx_pxSize;
			    cfg->ill_simVolume = cfg->mcx_simVolume;
			    cfg->ill_simVolume.z = 2 * cfg->ill_simVolume.z + 1;
			}
			return 4;
		//}
		/**if (argv[0][2] == 'p') {
			cfg->mcx_pxSize = atof(argv[1]);
			printf("\tMCX simulation pixel size = %s (mm)\n", argv[1]);
			if (cfg->ill_simVolume.y != 1)
			    cfg->ill_pxSize = cfg->mcx_pxSize;
			return 1;
		}*/
	case 's':
		printf("\tSensor dimensions = [%s %s] px\n", argv[1], argv[2]);
		cfg->sensorSize.x = atoi(argv[1]);
		cfg->sensorSize.y = atoi(argv[2]);
		return 2;

	case 'p':
		printf("\tSensor pixel size = %s (mm)\n", argv[1]);
		cfg->sensorpxSize = atof(argv[1]);
		return 1;

	case 'z':
		printf("\tZ planes to scan = %s. Start = %s (mm) Stop = %s (mm)\n", argv[1], argv[2], argv[3]);
		cfg->nzPlanes = atoi(argv[1]);
		cfg->zScan.x = atof(argv[2]);
		cfg->zScan.y = atof(argv[3]);
		return 3;
	
	case 'v':		
		cfg->spimVol= atoi(argv[1]);
		if (cfg->spimVol==1)
			printf("\tOutput image mode = SPIM \n");
		else
			printf("\tOutput image mode = OPT \n");
		return 1;

	case 'd':
	        printf("\tUsing detailed illumination  [xSize, zSize, pxSize] = [%s %s %s] \n", argv[1], argv[2], argv[3]);
		cfg->ill_simVolume.x = atoi(argv[1]);
		cfg->ill_simVolume.z = atoi(argv[2]);
		cfg->ill_simVolume.y = (unsigned int)1;
		cfg->ill_pxSize = atof(argv[3]);
		return 3;
	
	case 'i':
		printf("\tIllumination file name is %s\n", argv[1]);

		strBuffer = malloc(strlen(argv[1]) + 1);
		*strBuffer = '\\';
		strcpy(strBuffer + 1, argv[1]);
		cfg->illFileName = malloc(strlen(strBuffer));
		strcpy(cfg->illFileName, strBuffer);
		return 1;
	}
}

/**
*@brief Parse command line arguments
*/
int parseCommandLine(int argc, char* argv[], SPIMConfig*cfg){
	int i = 1;
	char* arg;
	while (i < argc-1) {
		//printf("arg %d: contains  string %s\n", i, argv[i]);
		//arg = argv[i];
		if (argv[i][0] == '-') {
			parseInput(&argv[i], cfg);
	
		}
		i++;
	}
	return 0;
}

/**
*@brief Match a field from the config file to a command
*/

int findCommand(char* field) {
	int nfields = sizeof(commandFields) / sizeof(commandFields[0]);
	for (int i = 0; i < nfields; i++)
		if (strcmp(field, configFields[i])==0)
			return i;
	return -1;
}

/**
*@brief Read config file from the working directory
*/
int readConfigFile(SPIMConfig* cfg) {
	FILE* fileID;
	char* configFile;
	char* config = "\\focusConfig.txt";
	if (cfg->configFileName != NULL) {
	    configFile = malloc(strlen(cfg->workingDir) + strlen(cfg->configFileName));
	    strcpy(configFile, cfg->workingDir);
	    strcat(configFile, cfg->configFileName);
	}
	else {
	    configFile = malloc(strlen(cfg->workingDir) + strlen(config));
	    strcpy(configFile, cfg->workingDir);
	    strcat(configFile, config);
	    cfg->configFileName = malloc(strlen(config));
	    strcpy(cfg->configFileName, config);
	}

	printf("Config file is: %s\n", configFile);

	fileID = fopen(configFile, "rb");
	if (NULL == fileID)
		return -1;
	char field[200];
	int fieldn;
	char* argv[5];
	char content[100];
	char* contents;
	int i=0;


	while (feof(fileID)==0) {	
		if (fscanf(fileID, "#%s", &field)) {
			fieldn = findCommand(field);
			if (fieldn > -1) {
				argv[i] = commandFields[fieldn];
				i++;
				fgets(&content, 100, fileID);
				contents = strtok(content, "\t\r\n ");
				while (contents != NULL) {
					argv[i] = contents;
					contents = strtok(NULL, "\t\r\n ");
					i++;
				}
				parseInput(argv, cfg);
			}
			
		}
		else
			fgets(field, 20, fileID);
		i = 0;
	}	
	fclose(fileID);

	/** Build output base filename*/
	char idx[5];
	const char* dataFile;

	//strcpy(field, cfg->workingDir);
	strcpy(&field, "\\");
	//strcat(field, "\\");
	strcat(field, cfg->fileBaseName);
	strcat(field, "_");
	itoa(cfg->sensorSize.x, idx, 10);
	strcat(field, idx);
	strcat(field, "_");
	itoa(cfg->sensorSize.y, idx, 10);
	strcat(field, idx);

	if (cfg->spimVol == 1) {
	    strcat(field, "_");
	    itoa(cfg->nzPlanes, idx, 10);
	    strcat(field, idx);
	}

	cfg->outputFileName = (char*)malloc(strlen(field));
	strcpy(cfg->outputFileName, field);
	cfg->sensorSize.z = (cfg->spimVol == 1) ? cfg->nzPlanes : 1;


	return 0;
}

/**
*@brief Read data from binary file. First 4 bytes contain the number of photons.
* Data organized as [px py pz vx vy vz p0x p0y p0z w], with Nphotons values per field, 32 bit float
*/
int importPhotonsData(const char *dataFile, SPIMConfig*cfg, float** photonsData) {
	FILE* fileID;
	size_t read;	
	size_t bufferElements;
	float *number;
	
	fileID= fopen(dataFile, "rb");
	if (NULL == fileID)
		return -1;
	
	read=fread(&cfg->Nphotons, sizeof(unsigned long int), 1, fileID);
	if (read == 0) return -2;
	printf("Nphotons = %d\n", cfg->Nphotons);
	bufferElements=cfg->photonSize*cfg->Nphotons;
	printf("Photons data buffer size = %d Bytes = %f GBytes\n", bufferElements, ((float)bufferElements)/1073741824.0f);
	*photonsData= (float*)calloc(bufferElements, sizeof(float));

	read=fread(*photonsData, sizeof(float), bufferElements, fileID);
	
	if (read == 0) return -2;
	fclose(fileID);
	cfg->photonsData = photonsData;
	return 0;
}
/**
*@brief Read data from.mch binary file. 
* Data organized as [px py pz vx vy vz p0x p0y p0z w], with Nphotons values per field, 32 bit float
*/
int readMchFile(const char* dataFile, SPIMConfig* cfg, float** photonsData) {
	FILE* fileID;
	size_t read;
	size_t bufferElements;
	float* number;
	int colcount; 
	unsigned int Nphotons = 0;
	printf("Loading .mch file: %s\n", dataFile);
	fileID = fopen(dataFile, "rb");
	if (NULL == fileID)
		return -1;
	fseek(fileID, 16, SEEK_CUR); /*Move pointer to colcount*/
	//fread(&colcount, sizeof(int), 1, fileID);/*########## Colcount can be what we call here photonSizeee#####################*/
	fread(&cfg->photonSize, sizeof(int), 1, fileID);/*########## Colcount can be what we call here photonSizeee#####################*/

	fseek(fileID, 8, SEEK_CUR); /*Move pointer to savedphotons*/	
	//fread(&cfg->Nphotons, sizeof(int), 1, fileID);
	fread(&Nphotons, sizeof(unsigned int), 1, fileID);
	fread(&cfg->mcx_pxSize, sizeof(float), 1, fileID);
	fseek(fileID, 28, SEEK_CUR); /*Move pointer to data*/


	/*Here now we need to read the data*/

	//read = fread(&cfg->Nphotons, sizeof(unsigned long int), 1, fileID);
	//if (read == 0) return -2;
	printf("Loading  %d photons\n", Nphotons);
	bufferElements = (unsigned long)cfg->photonSize * Nphotons;
	printf("Photons data buffer sizen in RAM = %d Bytes = %f GBytes\n", bufferElements, (4.0f*(float)bufferElements) / 1073741824.0f);
	*photonsData = (float*)calloc(bufferElements, sizeof(float));

	read = fread(*photonsData, sizeof(float), bufferElements, fileID);

	if (read == 0) return -2;
	fclose(fileID);
	//cfg->photonsData = photonsData;
	return Nphotons;
}

///**
//*@brief Read data from.mch binary file.
//* Data organized as [px py pz vx vy vz p0x p0y p0z w], with Nphotons values per field, 32 bit float
//*/
//int readMchHeader(const char* dataFile, focusConfig* fcfg) {
//	FILE* fileID;
//	size_t read;
//	size_t bufferElements;
//	float* number;
//	int colcount;
//
//	fileID = fopen(dataFile, "rb");
//	if (NULL == fileID)
//		return -1;
//	fseek(fileID, 16, SEEK_CUR); /*Move pointer to colcount*/
//	//fread(&colcount, sizeof(int), 1, fileID);/*########## Colcount can be what we call here photonSizeee#####################*/
//	fread(&fcfg->photonSize, sizeof(int), 1, fileID);/*########## Colcount can be what we call here photonSizeee#####################*/
//
//	fseek(fileID, 8, SEEK_CUR); /*Move pointer to savedphotons*/
//	//fread(&cfg->Nphotons, sizeof(int), 1, fileID);
//	fread(&fcfg->Nphotons, sizeof(unsigned int), 1, fileID);
//	fread(&fcfg->mcx_pxSize, sizeof(float), 1, fileID);
//	fclose(fileID);	
//	return 0;
//}

/**
*@brief Read illumination weights matrix data from binary file. 
*/
int importIlluminationData(SPIMConfig*cfg, float** illVolume) {
	FILE* fileID;
	size_t read;
	size_t bufferElements;
	float *number;
	const char* dataFile;
	char tmp[300];	
	//strcpy(tmp, cfg->workingDir);
	//strcat(tmp, "\\illuminationVolume.ill");
	//dataFile = &(*tmp);

	dataFile = malloc(strlen(cfg->workingDir) + strlen(cfg->illFileName));
	strcpy(dataFile, cfg->workingDir);
	strcat(dataFile, cfg->illFileName);

	printf("Illumination file path: %s\n",dataFile);
	fileID = fopen(dataFile, "rb");
	if (NULL == fileID)
		return -1;
	//if (cfg->ill_simVolume.z=!1)
	  //  bufferElements = cfg->ill_simVolume.x * cfg->ill_simVolume.y * (2 * cfg->ill_simVolume.z + 1);
	//else
	    bufferElements = cfg->ill_simVolume.x * cfg->ill_simVolume.y * cfg->ill_simVolume.z ;
	*illVolume = (float*)calloc(bufferElements, sizeof(float));

	read = fread(*illVolume, sizeof(float), bufferElements, fileID);

	if (read == 0) return -2;
	fclose(fileID);
	cfg->illVolume = *illVolume;
	return 0;
}

/**
*@brief Save in binary format the SPIM volume
*/
int exportVolume(float* imageVolume, SPIMConfig* cfg) {
	FILE* fileID;
	size_t written;
	char idx[5];
	int fidx = 0;
	const char* dataFile;
	char tmp[200];
	strcpy(tmp, cfg->workingDir);
	strcat(tmp, "\\");
	strcat(tmp, cfg->fileBaseName);
	strcat(tmp, "_");
	itoa(cfg->sensorSize.x, idx, 10);
	strcat(tmp, idx);
	strcat(tmp, "_");
	itoa(cfg->sensorSize.y, idx, 10);
	strcat(tmp, idx);
	
	if (cfg->spimVol == 1) {
		strcat(tmp, "_");
		itoa(cfg->nzPlanes, idx, 10);
		strcat(tmp, idx);
	}	
	strcat(tmp, ".mcspim");

	dataFile = &(*tmp);

	printf("Writting results in file: %s\n", tmp);
	fileID = fopen(dataFile, "wb");
	if (NULL == fileID)
		return -1;
	//if (cfg->spimVol == 1) 
	//	written = fwrite(imageVolume, sizeof(float), cfg->sensorSize.x * cfg->sensorSize.y * cfg->nzPlanes, fileID);
	//else
	written = fwrite(imageVolume, sizeof(float), cfg->sensorSize.x * cfg->sensorSize.y * cfg->sensorSize.z, fileID);

	if (written == 0) return -2;
	fclose(fileID);
	printf("File written succesfully!\n");
	return 0;

}