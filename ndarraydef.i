
%include "numpy.i"

%init %{
    import_array();
%}

%apply (unsigned char* IN_ARRAY1, int DIM1) {(unsigned char *bmap, int blen)}

%apply (unsigned char* INPLACE_ARRAY2, int DIM1, int DIM2) {(unsigned char *mask, int hm, int wm)}
%apply (unsigned short* INPLACE_ARRAY2, int DIM1, int DIM2) {(unsigned short *data, int hd, int wd)}
%apply (unsigned char* IN_ARRAY2, int DIM1, int DIM2) {(unsigned char *mask2, int hm2, int wm2)}

%apply (unsigned char* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(unsigned char *mask3d, int zm3d, int hm3d, int wm3d)}
%apply (unsigned short* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(unsigned short *data3d, int zd3d, int hd3d, int wd3d)}
%apply (unsigned char* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(unsigned char *ptmask, int npts, int hptm, int wptm)}

