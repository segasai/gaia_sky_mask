from cffi import FFI

ffibuilder = FFI()

ffibuilder.set_source(
    "_filler",
    r""" // passed to the real C compiler
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
   int nstar)
   {
      for(int i=0;i<nstar;i++)
      {
         float currad = MAX(rads[i],1);
         float xcen  = xs[i];
         float ycen = ys[i];
         int xl = CLIP(floor(xs[i]-currad),0,npix-1); 
         int xr = CLIP(ceil(xs[i]+currad),0,npix-1); 
         int yl = CLIP(floor(ys[i]-currad),0,npix-1);
         int yr = CLIP(ceil(ys[i]+currad),0,npix-1);
         for (int curx=xl; curx<=xr; curx++)
         {
            for(int cury=yl; cury<=yr; cury++)
            {
                if ((curx-xcen)*(curx-xcen)+(cury-ycen)*(cury-ycen)<currad*currad)
    {
          ima[curx*npix+cury]=1;
            }
            }
         }
    }
   }

    """,
    libraries=['m'],  #
)  # or a list of libraries to link with

ffibuilder.cdef("""     // some declarations from the man page

void procima(
   int *ima,
   int npix,
   float *xs,
   float *ys,
   float *rads,
   int nstar);
""")

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
