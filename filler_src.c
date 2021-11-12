#include <sys/types.h>
#include <math.h>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define CLIP(x,x1,x2) (MIN(MAX((x),(x1)),(x2)))

void procima(
	     int *ima,
	     int npix,
	     float *xs,
	     float *ys,
	     float *rads,
	     int nstar,
	     int binning)
{
  for(int i=0;i<nstar;i++)
    {
      float currad = MAX(rads[i], 1);
      float xcen  = xs[i];
      float ycen = ys[i];
      int xl = CLIP(floor(xs[i] - currad), 0, npix-1); 
      int xr = CLIP(ceil(xs[i] + currad), 0, npix-1); 
      int yl = CLIP(floor(ys[i] - currad), 0, npix-1);
      int yr = CLIP(ceil(ys[i] + currad), 0, npix-1);
      for (int curx=xl; curx<=xr; curx++)
	{
	  for(int cury=yl; cury<=yr; cury++)
            {
	      for(int ii=0;ii<binning ; ii++)
                { 
		  float curx1 = curx -0.5 +  0.5 * (ii+ 1.0)/binning;
		  for(int jj=0;jj<binning; jj++)
		    {
		      float cury1 = cury -0.5 +  0.5 * (jj+1.0 )/binning;
		      if ((curx1-xcen)*(curx1-xcen)+(cury1-ycen)*(cury1-ycen)<currad*currad)
			{
			  ima[curx*npix+cury] = ima[curx*npix+cury]
			     | (1<<(ii+binning*jj));
			}
		    }
		}
	    }
	}
    }
}



void procima_pix(
	     int *ima,
	     int npix,
	     float *xs,
	     float *ys,
	     int nstar,
	     int binning)
{
  for(int i=0;i<nstar;i++)
    {
      float xcen  = xs[i];
      float ycen = ys[i];
      int xi = round(xs[i]); 
      int yi = round(ys[i]);
      if ((xi<0)|| (xi>=npix) ||(yi<0)||(yi>=npix))
	{continue;}
      for(int ii=0;ii<binning ; ii++)
	{ 
	  float curx1 = xi -0.5 +  0.5 * (ii+ 1.0)/binning;
	  for(int jj=0;jj<binning; jj++)
	    {
	      float cury1 = yi -0.5 +  0.5 * (jj+1.0 )/binning;
	      if (((curx1-xcen)<.5/binning)&&((cury1-ycen)<.5/binning))
		{
		  ima[xi*npix+yi] = ima[xi*npix+yi]
		    | (1<<(ii+binning*jj));
		}
	    }
	}
    }
}
