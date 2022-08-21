%module imagetools

%include "stddef.i"

%{
    #define SWIG_FILE_WITH_INIT
    #include "cpp/imagetools.h"
%}

%include "ndarraydef.i"

%include "cpp/imagetools.h"
