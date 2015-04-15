/*
	To create an EAM_Dynamo setfl model, where the user will give all the info.  

	Contributor: Mingjian Wen
*/

#include<stdio.h>
#include<string.h>


int main() {

		/* local variables */
	char SysCmd[256];
	char tmpstring[256];
	char flag_make[256];
	int i,j;
  FILE *pFile;

	char model_name[256];
  char model_driver_name[256]; /* "EAM_CubicCompleteSpline__MD_000000111111_000" */
  
	int Nspecies = 1;

	int NumRho = 15;
	double RhoStep = 0.1;
	int NumR	 = 15;
	double RStep = 0.1; 
  double cutoff = 	0.0;

/*************************************************************************
* read in data from stdin 
**************************************************************************/

	printf("\nThis util is used to generated KIM-compliant EAM Dynamo model" 
				" the parameter file is in Dynamo `setfl' format. You can find an explanation of"
 				" the format at: <http://enpub.fulton.asu.edu/cms/potentials/submain/format.htm>."
				" If you want to generate initial potential for potfit-kim, the number"
				" of species types, the number of rho points and the number of r points"
				" should be exactly the same as those in potfit input files. Others"
				" will be handled by potfit-kim automically.\n\n");


	printf("Enter model name (e.g. ``EAM_Dynamo_tmp_model''."
				 " Type `default' or `d' to use this one):\n");
  fgets(tmpstring, 80, stdin);
	sscanf(tmpstring, "%s", model_name);
	if(! strcmp(model_name, "default") || ! strcmp(model_name, "d") )
		strcpy(model_name, "EAM_Dynamo_tmp_model");


	printf("\nEnter model driver name (e.g. ``EAM_Dynamo__MD_120291908751_001''."
				 " Type `default' or `d' to use this one):\n");
  fgets(tmpstring, 80, stdin);
	sscanf(tmpstring, "%s", model_driver_name);
	if(! strcmp(model_driver_name, "default") || ! strcmp(model_driver_name, "d") )
		strcpy(model_driver_name, "EAM_Dynamo__MD_120291908751_001");

	printf("\nEnter number of species types. (e.g. 1).\n");
  fgets(tmpstring, 80, stdin);
	sscanf(tmpstring,"%d", &Nspecies);


	char species[Nspecies][3];
	int atomic_number[Nspecies];
	double atomic_mass[Nspecies];
	double lattice_const[Nspecies];
	char lattice_type[Nspecies][15];


	printf("\nEnter number of rho points, rho data step, number of r points, r step, "
				 "and cutoff radius. (e.g. ``15  0.1  15  0.1  8.0'').\n");
  fgets(tmpstring, 80, stdin);
	sscanf(tmpstring,"%d %lf %d %lf %lf", &NumRho, &RhoStep, &NumR, &RStep, &cutoff);


	for(i = 0; i < Nspecies; i++) {
		printf("\nEnter element type, atomic number, atomic mass, lattice constant, and "
					" lattice type of species %d. (e.g. ``Al  13  26.982  4.032  fcc''."
					" Type `default' or `d' to use ``El  0  0.0  0.0  lattice_type''):\n", i+1);
  	fgets(tmpstring, 80, stdin);
		if(! strncmp(tmpstring, "default", 7) || ! strncmp(tmpstring, "d", 1) ) {
			strcpy(species[i], "El");
			atomic_number[i] = 0;
			atomic_mass[i] = 0.0;	
			lattice_const[i] = 0.0;
			strcpy(lattice_type[i], "lattice_type");	
		}
		else	
		sscanf(tmpstring,"%s %d %lf %lf %s", species[i], &atomic_number[i],
						  &atomic_mass[i],	&lattice_const[i], lattice_type[i]);
	}

	printf("\nDo you want to make the model right now? Type `yes' or 'no'.\n");
  fgets(tmpstring, 80, stdin);
	sscanf(tmpstring,"%s", flag_make);

	

/*************************************************************************
* Create model dir
**************************************************************************/

	/* create model directoy */
/*NOTE(change) it would be better check whether the dir exits or not */
	strcpy(SysCmd, "mkdir ");
	strcat(SysCmd, model_name);
	system(SysCmd);


/*************************************************************************
* Create Makefile 
**************************************************************************/
	sprintf(tmpstring, "%s/Makefile", model_name);
	pFile = fopen(tmpstring,"w");
	if(pFile == NULL)
		error(1, "Error creating file: KIM Makefile. %s %d\n", __FILE__, __LINE__);
	else {
		fprintf(pFile,
			"#\n"
			"# Copyright (c) 2013--2014, Regents of the University of Minnesota.\n"
			"# All rights reserved.\n"
			"#\n"
			"# Contributors:\n"
			"#    Ryan S. Elliott\n"
			"#    Ellad B. Tadmor\n"
			"#    Mingjian Wen\n"
			"#\n\n\n" );	

		fprintf(pFile,"# load all basic KIM make configuration\n" );	
		fprintf(pFile,"ifeq ($(wildcard ../Makefile.KIM_Config),)\n"
	  	"  $(error ../Makefile.KIM_Config does not exist. Something is wrong with"
			" your KIM API package setup)\n"
			"endif\n"
			"include ../Makefile.KIM_Config\n\n" );	

		fprintf(pFile,"# set model driver specific details\n" );
		fprintf(pFile,"MODEL_DRIVER_NAME    := %s\n", model_driver_name);	
		fprintf(pFile,"MODEL_NAME           := %s\n", model_name);	

		for (i = 0; i < Nspecies; i++) {
			sprintf(tmpstring, "SPECIES_%03d_NAME", i+1);		
			fprintf(pFile,"%s     := %s\n", tmpstring, species[i]);	
		}
		sprintf(tmpstring, "%s.params", model_name);			
		fprintf(pFile,"PARAM_FILE_001_NAME  := %s\n\n", tmpstring);	
	
		fprintf(pFile,
			"# APPEND to compiler option flag lists\n"
			"#FFLAGS   +=\n"
			"#CFLAGS   +=\n"
			"#CXXFLAGS +=\n"
			"#LDFLAGS  +=\n"
			"#LDLIBS   +=\n\n");	

		fprintf(pFile, "# load remaining KIM make configuration\n"
		"include $(KIM_DIR)/$(builddir)/Makefile.ParameterizedModel\n");	
	}
	fflush(pFile);
	fclose(pFile);



/*************************************************************************** 
* Create parameter file EAM
/***************************************************************************/ 

	sprintf(tmpstring, "%s/%s.params", model_name, model_name);	
	pFile = fopen(tmpstring,"w");
	if(pFile == NULL)
		error(1, "Error creating file: KIM Makefile. %s %d\n", __FILE__, __LINE__);
	else {
  	/* 3 lines of comments */
		fprintf(pFile, "# EAM potential in Dynamo `setfl' format, generated by potfit-kim.\n");		
		fprintf(pFile, "# For more info about KIM, visit: https://openkim.org\n");				
		fprintf(pFile, "# For more info about potfit, visit: http://potfit.sourceforge.net\n");		
		
    /* number of atoms types and atom species*/
		fprintf(pFile, "    %d", Nspecies);	
		for (i = 0; i < Nspecies; i++) {
			fprintf(pFile, "  %s", species[i]);
		}
		fprintf(pFile, "\n");

		/* number of data points, delta  rcut */
		fprintf(pFile, "  %d %22.15e %d %22.15e %22.15e\n", NumRho, RhoStep, 
					  NumR, RStep, cutoff);

		/* embedding data, density data and their header */
		for (i = 0; i < Nspecies; i++) {

			/* the header */
			fprintf(pFile, "  %d  %f  %f  %s\n", atomic_number[i], atomic_mass[i],
																					 lattice_const[i], lattice_type[i]);

			/* fill embedding data with 0.0 */
			for (j = 0; j < NumRho; j++) {
				fprintf(pFile, "  %2.1f", 0.0);
			}
			fprintf(pFile, "\n");
			/* fill density data with 0.0 */
			for (j = 0; j < NumR; j++) {
				fprintf(pFile, "  %2.1f", 0.0);
			}
			fprintf(pFile, "\n");
		}

		/* pair data */
		for (i = 0; i < Nspecies*(Nspecies + 1)/2; i++) {	
			/* fill pair data with 0.0 */
			for (j = 0; j < NumR; j++) {
				fprintf(pFile, "  %2.1f", 0.0);
			}
			fprintf(pFile, "\n");
		}
	}	
	fflush(pFile);
	fclose(pFile);


/***************************************************************************
* System call to make the model 
***************************************************************************/
	/* create Makefile.KIM_Config and make model */
 	if(!strncmp(flag_make, "yes", 3)) {
		system("kim-api1-build-config --makefile-kim-config > Makefile.KIM_Config");
		printf("\n\nMaking KIM model ... starting\n");
	  sprintf(SysCmd,"cd %s; make", model_name);
		system(SysCmd);
	}

	return 0;
}
