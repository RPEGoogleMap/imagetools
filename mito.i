%module mito

%include "stddef.i"

%{
    #define SWIG_FILE_WITH_INIT
    #include "cpp/mito.h"
%}

%include "ndarraydef.i"

%include "cpp/mito.h"
