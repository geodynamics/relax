# include <gmt.h>

/* Fortran callable routine to write a grd file in pixel registration */
/* June 23, 1995 - David Sandwell */
/* Revised for GMT3.4 December 28, 2002 - David Sandwell */
/* Modified for node registration - March 19, 2008 - Sylvain Barbot */

void writegrd(rdat,nx,ny,rlt0,rln0,dlt,dln,rland,rdum,title,fileout)

  float *rdat;            /* real array for output */
  int *nx;                /* number of x points */
  int *ny;                /* number of y points */
  double *rlt0;            /* starting latitude */
  double *rln0;            /* starting longitude */
  double *dlt;             /* latitude spacing */
  double *dln;             /* longitude spacing */
  double *rland;            /* land value */
  double *rdum;            /* dummy value */
  char  *title;           /* title */
  char  *fileout;         /* filename of output file */
  
  {
   int i;
   double xmin, xmax, xinc, ymin, ymax, yinc, zmin, zmax;
   int update = FALSE;
   struct GRD_HEADER grd;
   int argc = 0;
   char **argv = NULL;

/* Initialize with default values */
 
   GMT_grdio_init(); 
   GMT_make_dnan(GMT_d_NaN);
   GMT_make_fnan(GMT_f_NaN);
   
   GMT_grd_init(&grd, argc, argv, update);

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

   zmin = fabs((double)*rdum);
   zmax = -fabs((double)*rdum);

   for (i = 0; i < *nx * *ny; i++) {
     if((rdat[i] == *rdum) || (rdat[i] == *rland)) rdat[i] = GMT_f_NaN;
     else {
        if(rdat[i] < zmin) zmin = rdat[i];
        if(rdat[i] > zmax) zmax = rdat[i];
     }
   }

/* update the header using values passed */

   strncpy(grd.title,title,80); 
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

/*  write the file */
   GMT_write_grd(fileout, &grd, rdat, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE );
   
  }
