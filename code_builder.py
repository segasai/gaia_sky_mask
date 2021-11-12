from cffi import FFI

ffibuilder = FFI()

ffibuilder.set_source(
    "_filler",
    open('filler_src.c').read(),
    libraries=['m'],
)

ffibuilder.cdef("""
void procima(
   int *ima,
   int npix,
   float *xs,
   float *ys,
   float *rads,
   int nstar,
int binning);
void procima_pix(
   int *ima,
   int npix,
   float *xs,
   float *ys,
   int nstar,
int binning);
""")

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
