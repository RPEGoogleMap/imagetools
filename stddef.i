// This tells SWIG to treat char ** as a special case
%typemap(in) char ** {
  /* Check if is a list */
  if (PyList_Check($input)) {
    long long sz = PyList_Size($input);
    long long i = 0;
    $1 = (char **) malloc((sz+1)*sizeof(char *));
    // $1 = new char *[sz+1];
    for (i = 0; i < sz; i++) {
      PyObject *o = PyList_GetItem($input, i);
      if (PyString_Check(o)) {
      	$1[i] = PyString_AsString(o);
      } else
      if (PyUnicode_Check(o)) {
        $1[i] = PyBytes_AS_STRING(PyUnicode_AsEncodedString(o, "utf-8", "Error ~"));
      } else {
        PyErr_SetString(PyExc_TypeError, "list must contain strings");
        SWIG_fail;
      }
    }
    $1[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError, "not a list");
    SWIG_fail;
  }
}

// This cleans up the char ** array we malloc'd before the function call
%typemap(freearg) char ** {
  free($1);
  // delete [] $1;
}

%include "std_vector.i"
%include "std_string.i"

namespace std {
   %template(vectori) vector<int>;
   %template(vectorl) vector<long long>;
   %template(vectord) vector<double>;
   %template(vectori2d) vector<vector<int>>;
   %template(vectori3d) vector<vector<vector<int>>>;
};

