/*-----------------------------------------------------------------------
 * ! Copyright 2007, 2008, 2009 Sylvain Barbot
 * !
 * ! This file is part of RELAX
 * !
 * ! RELAX is free software: you can redistribute it and/or modify
 * ! it under the terms of the GNU General Public License as published by
 * ! the Free Software Foundation, either version 3 of the License, or
 * ! (at your option) any later version.
 * !
 * ! RELAX is distributed in the hope that it will be useful,
 * ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 * ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * ! GNU General Public License for more details.
 * !
 * ! You should have received a copy of the GNU General Public License
 * ! along with RELAX.  If not, see <http://www.gnu.org/licenses/>.
 * !
 * ! \author Sagar Masuti 
 * !----------------------------------------------------------------------*/
#include "config.h"

#ifdef PAPI_PROF
/*-------------------------------------- Includes ----------------------------------------------*/

#include <papi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------------------------*/


/*--------------------------------- General Declarations ---------------------------------------*/

#define MAX_NAME_LENGTH 16

/*----------------------------------------------------------------------------------------------*/


/*------------------------------- Internal structure -------------------------------------------*/
typedef struct _stProfData 
{
	long 			lStartTime ;   	/* The start time  */
	char 			*pcName ;	/* The name of the timer */
	struct _stProfData 	*pNext ;        /* Pointer to the next element in the list */
} ST_PROF_DATA ;

/* -------------------------------------------------------------------------------------------- */

/* ---------------------------------- Global variables ---------------------------------------- */

ST_PROF_DATA *pHead = NULL ; 		/* The head node of the list */ 

int iInitialized = 0 ; 			/* To check whether the PAPI is already initialized */

/*--------------------------------------------------------------------------------------------- */

/* ------------------- Forward declaration of static functions -------------------------------- */

static void 		profDeleteTimer (char	*pcName) ; 	

static ST_PROF_DATA* 	profGetTimer (char 	*pName) ;

static int	 	profAddTimer (char	*pName,
                  		      long      lTime) ;

/*----------------------------------------------------------------------------------------------*/


/*----------------------------------- Main Functions -------------------------------------------*/

/**
 * @brief : 	This function starts recording the time information. This function can be called 
		from other files. This function initializes PAPI library. 
 * @param       pcProfName[in]  Pointer to the name of the timer.
 * @retval 			Doesnt return anything. 
 */

void papistartprofiling_ (char	pcProfName[])
{
	int 	iRet ;
	long  	lTime ;

	pcProfName[MAX_NAME_LENGTH] = '\0' ;	
	iRet = 0 ; 

	/* Check if the PAPI library is already initialized, if not initialize it */
	if (0 == iInitialized)
	{
		if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
		{
			return ;
		}
	 	/* Indicate the library is initialized already */
		iInitialized = 1 ;	
	}

	
	/* Get the time before any further processing */
	lTime = PAPI_get_real_usec () ;
	
	if (NULL == profGetTimer (pcProfName))
	{
		/*Add it */
		profAddTimer (pcProfName, lTime) ;
	}
	return ;
}


void papiendprofiling_ (char	 pcName[])
{
	int 		iRetVal ;
	long 		lEndTime ; 
	ST_PROF_DATA 	*pTimer ;
 	float 		fElapsed ;

	iRetVal = 0 ;
	pTimer = NULL ;
	fElapsed = 0 ;	
	pcName[MAX_NAME_LENGTH] = '\0' ;	

	lEndTime = PAPI_get_real_usec () ;
	
	pTimer = profGetTimer (pcName) ;
	if (NULL == pTimer)
	{
		/*Timer doesnt exist*/
		return ;
	}
	
	/*Calculate the time*/
	fElapsed  = ((float)(lEndTime - pTimer->lStartTime)) / 1000000.0 ;
	printf ("Time taken to execute %s : %f\n", pcName, fElapsed) ;

	/*Remove the timer*/
	profDeleteTimer (pcName) ;	

	return ;
}

/*----------------------------------------------------------------------------------------------*/

/*---------------------------------- Utility Functions -----------------------------------------*/

/**
 * @brief :  	 This function finds the timer with the given name and return the pointed to the timer 
 *           	 in the linked list.
 * @param  	 pcName[in] 	 Pointer to the name of the timer.
 * @retval	 		 Returns the pointer to the timer or NULL in case not found.
 *
 **/

static ST_PROF_DATA* profGetTimer (char *pcName)
{
	int 		iRetValue ;
	ST_PROF_DATA 	*pTemp ;
	ST_PROF_DATA 	*pTimer ;
	
	iRetValue = 0 ;
	pTemp = pHead ;
	pTimer = NULL ;
	
	while (NULL != pTemp)
	{
		if (0 == strcmp (pcName, pTemp->pcName))
		{
			pTimer = pTemp ;
			break ;
		}
		pTemp = pTemp->pNext ;
	}
	
	return pTimer ;
}

/**
 * @brief : This function adds the timer information to the linked list.
 * @param 	pcName[in]	Pointer to the name of the timer.
 * @param 	lTime		Its the start time of the timer to be recorded.
 * @retval 			Returns the integer value 0 to indicate failure and 1 to
 *				indicate the timer is successfully added to the list.
 */

static int profAddTimer (char 		*pcName, 
     		         long 		lTime)
{
	int 		iRetVal ;
	ST_PROF_DATA 	*pTemp ;

	iRetVal = 0 ;
	pTemp = NULL ;
	
	if (NULL != profGetTimer (pcName))
	{
		/*Timer already exits*/
		return 0 ;		
	}
	
	pTemp = (ST_PROF_DATA *) calloc (1, sizeof (ST_PROF_DATA)) ;
	if (NULL == pTemp)
	{
		printf ("profAddTimer : Failed to allocate memory\n") ;
		return 0 ;
	}
	
	pTemp->pcName = (char *) calloc ((strlen (pcName) + 1), sizeof (char)) ;
	if (NULL == pTemp->pcName)
	{
		printf ("profAddTimer : Failed to allocate memory\n") ;
                return 0 ;
	}

	strcpy (pTemp->pcName, pcName) ; 
	pTemp->lStartTime = lTime ;

	if (NULL == pHead)
	{
		pHead = pTemp ;
		pHead->pNext = NULL ;
	}
	else
	{
		pTemp->pNext = pHead ;
		pHead = pTemp ;
	}

	return 1 ;
}

/**
 * @brief : 	This function deletes the timer node from the linked list.
 * @param       pcName[in]      Pointer to the name of the timer.
 * @retval                      Doesnt return any value.
 *
 */

static void profDeleteTimer (char	*pcName) 	
{	
	int iRetValue ;
        ST_PROF_DATA *pTemp ;
        ST_PROF_DATA *pPrev ;

        iRetValue = 0 ;
        pTemp = pHead ;
	pPrev = NULL ;

	if (0 == strcmp (pcName, pTemp->pcName))
	{
		pHead = pTemp->pNext ; 
		free (pTemp->pcName) ;
		free (pTemp) ;
	}
        while (NULL != pTemp)
        {
                if (0 == strcmp (pcName, pTemp->pcName))
                {
			pPrev->pNext = pTemp->pNext ;
			free (pTemp->pcName) ;
                	free (pTemp) ;  			             
                }
		pPrev = pTemp ;
                pTemp = pTemp->pNext ;
        }
}

/*----------------------------------------------------------------------------------------------*/

#endif
