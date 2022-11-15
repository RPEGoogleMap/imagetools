%module laminb

%include "stddef.i"

%{
    #define SWIG_FILE_WITH_INIT
    #include "cpp/laminb.h"
%}

%include "ndarraydef.i"

%include "cpp/laminb.h"
