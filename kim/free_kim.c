/***************************************************************************
*
* KIM-potfit
*
*	kim_free.c 
* 
* Defining the functions that free memory realted to KIM. e.g., free 
* neighbor list, free KIM object and destroy KIM model etc.
*
* Contributor: Mingjian Wen
*
***************************************************************************/


#include "../potfit.h"				/* to use `nconf' etc.	*/
#include "free_kim.h"

/* added */
/***************************************************************************
* 
* Free memory of KIM objects and neighbor lists 
*
***************************************************************************/

void free_KIM(void)
{
	/* local variables */
	int status;
	int i,j;
	NeighObjectType* NeighObject;
	void* pkim;

	for (i = 0; i < nconf; i++) {		
		pkim = pkimObj[i];
		
		/* call model destroy */
		status = KIM_API_model_destroy(pkim);
    if (KIM_STATUS_OK > status) 
			KIM_API_report_error(__LINE__, __FILE__,"destroy", status);

		/* free memory of neighbor list */
		/* first get neighbor object from KIM */
	  NeighObject = (NeighObjectType*) KIM_API_get_data(pkim, "neighObject", &status);
  	if (KIM_STATUS_OK > status) 
			KIM_API_report_error(__LINE__, __FILE__,"get_data", status);
		/* then free the stuff in neighbor list */	
		free(NeighObject->NNeighbors);	
		free(NeighObject->neighborList);	
		free(NeighObject->RijList);	
		free(NeighObject->BeginIdx);	
		/* free the neighbor list itself */
		free(NeighObject);	
    
    /* free paramlist */
    free_free_param (&FreeParamAllConfig[i]);
		/* free memory of KIM objects */
  	KIM_API_free(pkim, &status);
 		if (KIM_STATUS_OK > status) 
			KIM_API_report_error(__LINE__, __FILE__,"destroy", status);
	}
	free(FreeParamAllConfig);
  free(pkimObj);
	/* every thing is great */
}


void free_free_param (FreeParamType* FreeParam) 
{
	/*local variable*/
	int i;

	for(i = 0; i < FreeParam->Nparam; i++) {
		free(FreeParam->name[i]);	
		free(FreeParam->shape[i]);
	}
	free(FreeParam->name);	
	free(FreeParam->value);
	free(FreeParam->shape);	
	free(FreeParam->rank); 
	if (FreeParam->nestedvalue != NULL) {
		free(FreeParam->nestedvalue);
	}
}
