#include <stdio.h>
#include <stdlib.h>
#include <proj_api.h>
#include <string.h>

/*
 * proj routine to convert arrays of UTM coordinates
 * to longitude/latitude using the PROJ.4 library
 *
 * to do: check the output in the south hemisphere
 *
 * sylvain barbot (22/05/10) - original form
 */

void proj_(double *x, double *y, int * n, 
           double * lon0, double * lat0, int * zone) {

  projPJ pj_utm, pj_latlong;
  int i;
  char zonestr[3];
  char cmd_utm[100], cmd_latlong[100];
  char * to;

  // convert integer zone to string zone
  i=sprintf(zonestr, "%d", (*zone));

  // construct conversion command (+proj=utm +zone=11)
  to = stpcpy(cmd_utm,"+proj=utm +zone=");
  to = stpcpy(to,zonestr);
  //printf("%s\n",cmd_utm);

  // construct conversion command (+proj=latlong +zone=11)
  to = stpcpy(cmd_latlong,"+proj=latlong +zone=");
  to = stpcpy(to,zonestr);
  //printf("%s\n",cmd_latlong);

  if (!(pj_utm = pj_init_plus(cmd_utm)) ){
    printf("error initializing input projection driver. exiting.");
    exit(1);
  }
  if (!(pj_latlong = pj_init_plus(cmd_latlong)) ){
    printf("error initializing output projection driver. exiting.");
    exit(1);
  }

  // convert to radians
  (*lon0)*=DEG_TO_RAD;
  (*lat0)*=DEG_TO_RAD;

  pj_transform(pj_latlong, pj_utm, 1, 1, lon0, lat0, NULL);

  // add UTM coordinates of the origin
  for (i=0;i<(*n);i++){
    x[i]+=(*lon0);
    y[i]+=(*lat0);
  }
  pj_transform(pj_utm, pj_latlong, (*n), 1, x, y, NULL);

  // convert longitude and latitude to degrees
  for (i=0;i<(*n);i++){
    x[i]*=RAD_TO_DEG;
    y[i]*=RAD_TO_DEG;
  }
}
