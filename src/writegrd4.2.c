/************************************************************************
* writegrd routine to write a grd file in pixel registration            *
************************************************************************/
/************************************************************************
* Creator: David T. Sandwell    Scripps Institution of Oceanography    *
* Date   : 06/23/95             Copyright, David T. Sandwell           *
************************************************************************/
/************************************************************************
* Modification history:                                                 *
*   Revised for GMT3.4 December 28, 2002 - David Sandwell               *
*   Revised for GMT4.2 May 10, 2007 - David Sandwell                    *
*   Modified for pixel registration April 18, 2008 - Sylvain Barbot     *
************************************************************************/

# include <math.h>
# include <gmt.h>

void writegrd_(rdat,nx,ny,rlt0,rln0,dlt,dln,rland,rdum,title,fileout)

  float *rdat;            /* real array for output */
  int *nx;                /* number of x points */
  int *ny;                /* number of y points */
  double *rlt0;            /* starting latitude */
  double *rln0;            /* starting longitude */
  double *dlt;             /* latitude spacing */
  double *dln;             /* longitude spacing */
  double *rland;           /* land value */
  double *rdum;            /* dummy value */
  char  *title;           /* title */
  char  *fileout;         /* filename of output file */
  
  {
   int i;
   double xmin, xmax, xinc, ymin, ymax, yinc, zmin, zmax;
   struct GRD_HEADER grd;
   int argc2 = 1;
   char *argv2[2] = {"writegrd",0};

/* Initialize with default values */
 
   GMT_begin (argc2,argv2);
   GMT_grd_init(&grd, argc2, argv2, FALSE);

/* Calculate header parameters */

   xmax = *rln0 + ((*nx)-1) * *dln;
   xmin = *rln0;
   if(xmax < xmin) {
     xmin = xmax;
     xmax = *rln0;
     }
   xinc = fabs((double)*dln);
   ymax = *rlt0 + ((*ny)-1) * *dlt;
   ymin = *rlt0;
   if(ymax < ymin) {
     ymin = ymax;
     ymax = *rlt0;
     }
   yinc = fabs((double)*dlt);

/*  calculate zmin and zmax and zinc and set dummy values to NaN. */

   zmin = +fabs((double)*rdum);
   zmax = -fabs((double)*rdum);

   for (i = 0; i < *nx * *ny; i++) {
     if((rdat[i] == *rdum) || (rdat[i] == *rland)) rdat[i] = GMT_f_NaN;
     else {
        if(rdat[i] < zmin) zmin = rdat[i];
        if(rdat[i] > zmax) zmax = rdat[i];
     }
   }

/* update the header using values passed */

   strncpy(grd.title,title,GRD_TITLE_LEN); 
   grd.nx = *nx;
   grd.ny = *ny;
   grd.node_offset = FALSE;
   grd.x_min = xmin;
   grd.x_max = xmax;
   grd.x_inc = xinc;
   grd.y_min = ymin;
   grd.y_max = ymax;
   grd.y_inc = yinc;
   grd.z_min = zmin;
   grd.z_max = zmax;

/* grd.type = 10;
   grd.z_id = 15;
   grd.ncid = 15;*/

/*  write the file */

   GMT_write_grd(fileout, &grd, rdat, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE);

/*   GMT_end (argc2,argv2); */

  }

