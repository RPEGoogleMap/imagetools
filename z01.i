%module z01

%include "stddef.i"

%{
    #define SWIG_FILE_WITH_INIT
    #include "cpp/z01.h"
%}

%include "ndarraydef.i"

%include "cpp/z01.h"
