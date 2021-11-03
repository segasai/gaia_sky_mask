from cffi import FFI

ffibuilder = FFI()

ffibuilder.set_source(
    "_filler",
    open('filler_src.c').read(),
    libraries=['m'],  #
)  # or a list of libraries to link with

ffibuilder.cdef("""     // some declarations from the man page

void procima(
   int *ima,
   int npix,
   float *xs,
   float *ys,
   float *rads,
   int nstar,
int binning);
""")

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
