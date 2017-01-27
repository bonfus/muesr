#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include <numpy/arrayobject.h>

// support numpy 1.6 - this macro got renamed and deprecated at once in 1.7
#ifndef NPY_ARRAY_IN_ARRAY
#define NPY_ARRAY_IN_ARRAY NPY_IN_ARRAY
#endif

#if NPY_API_VERSION <= NPY_1_6_API_VERSION
#define PyArray_SHAPE PyArray_DIMS
#endif

static char module_docstring[] = "This module provides the Fields and the DipolarTensor functions.";
static char py_lfcext_fields_docstring[] = "The Fields calculation.";
static char py_lfcext_dt_docstring[] = "Dipolar tensor calculation.";

//Arbitrary size sum for incommensurate magnetic orders
void INCASS(const double *, const double *, const double *, const double *,
          const double *, const int * , const double *, const double , 
          const unsigned int , const double , unsigned int, unsigned int,
          double *, double *, double *);
          
//Arbitrary size sum with rotations
void RASS(const double *, const double *, const double *, const double *,
          const double *, const int * , const double *, const double , 
          const unsigned int , const double, unsigned int , const double *, 
          unsigned int , double *, double *, double *);
          
// Dipolar tensor with the positions given in in_positions
void DT(const double *in_positions, 
          const double *in_muonpos, const int * in_supercell, const double *in_cell, 
          const double radius, unsigned int size,
          double *out_field);
          

//Arbitrary size sum
void  ASS(const double *in_positions, 
          const double *in_fc, const double *in_K, const double *in_phi,
          const double *in_muonpos, const int * in_supercell, const double *in_cell, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius, 
          unsigned int size,
          double *out_field_cont, double *out_field_dip, double *out_field_lor);          


static PyObject * py_lfcext_fields(PyObject *self, PyObject *args) {
  /* input variables */
  char* calc_type = NULL;
  unsigned int nnn=0;
  unsigned int nangles=0;
  
  double r=0.0;
  double rcont=0.0;
    
  PyObject *opositions, *oFC, *oK, *oPhi;
  PyObject *omu, *osupercell, *ocell;
  PyObject *orot_axis = NULL;
  
  PyArrayObject *positions, *FC, *K, *Phi;
  PyArrayObject *mu, *supercell, *cell;
  PyArrayObject *rot_axis = NULL;

  /* local variables */
  int num_atoms=0;
  int icalc_type=0;
  
  
  double * in_positions=NULL;
  double * in_fc=NULL;
  double * in_K=NULL;
  double * in_phi=NULL;
  double * in_muonpos=NULL;
  int * in_supercell=NULL;
  double * in_cell=NULL;
  double * in_axis = NULL;

  /* return variables */
  double * dip=NULL, *cont=NULL, *lor=NULL;
  npy_intp *out_dim=NULL;
  PyArrayObject *odip, *ocont, *olor;



  /* put arguments into variables */
  if (!PyArg_ParseTuple(args, "sOOOOOOOdId|IO", &calc_type, 
                                        &opositions, &oFC, &oK, &oPhi,
                                        &omu, &osupercell, &ocell,
                                        &r,&nnn,&rcont,
                                        &nangles,&orot_axis)) {
    return NULL;
  }
  
  /* turn inputs into numpy array types */
  positions = (PyArrayObject *) PyArray_FROMANY(opositions, NPY_DOUBLE, 2, 2,
                                              NPY_ARRAY_IN_ARRAY);
  FC = (PyArrayObject *) PyArray_FROMANY(oFC, NPY_COMPLEX128, 2, 2,
                                              NPY_ARRAY_IN_ARRAY);
  K = (PyArrayObject *) PyArray_FROMANY(oK, NPY_DOUBLE, 1, 1,
                                              NPY_ARRAY_IN_ARRAY);
  Phi = (PyArrayObject *) PyArray_FROMANY(oPhi, NPY_DOUBLE, 1, 1,
                                              NPY_ARRAY_IN_ARRAY);
  mu = (PyArrayObject *) PyArray_FROMANY(omu, NPY_DOUBLE, 1, 1,
                                              NPY_ARRAY_IN_ARRAY);
                                              
  supercell = (PyArrayObject *) PyArray_FROMANY(osupercell, NPY_INT32,
                                                   1, 1, NPY_ARRAY_IN_ARRAY);
  cell = (PyArrayObject *) PyArray_FROMANY(ocell, NPY_DOUBLE, 2, 2,
                                             NPY_ARRAY_IN_ARRAY);

  /* Validate data */
  if (!calc_type || !positions || !FC || !K || !Phi || !mu || !supercell || !cell) {
    Py_XDECREF(positions);
    Py_XDECREF(FC);
    Py_XDECREF(K);
    Py_XDECREF(Phi);
    Py_XDECREF(mu);
    Py_XDECREF(supercell);
    Py_XDECREF(cell);
    PyErr_Format(PyExc_RuntimeError,
                    "Error parsing numpy arrays.");                 
    return NULL;
  }


  
  
  
  /* Select calculation type */
  if (strcmp(calc_type, "s")==0 || strcmp(calc_type, "sum")==0) {
    icalc_type = 1;
    nangles = 1;
  }  else if (strcmp(calc_type, "r")==0 || strcmp(calc_type, "rotate")==0) {
    icalc_type = 2;
  }  else if (strcmp(calc_type, "i")==0 || strcmp(calc_type, "incommensurate")==0) {
    icalc_type = 3;
  }  else  {
    Py_DECREF(positions);
    Py_DECREF(FC);
    Py_DECREF(K);
    Py_DECREF(Phi);
    Py_DECREF(mu);
    Py_DECREF(supercell);
    Py_DECREF(cell);    
    PyErr_Format(PyExc_ValueError,
                   "Valid calculations are 's', 'r', 'i'.  Unknown value %s", calc_type);
    return NULL;
  }

  /* Check optional arguments are present for calculation 'i' or 'r' */
  if (icalc_type >= 2) {
      if(nangles==0) {
        Py_DECREF(positions);
        Py_DECREF(FC);
        Py_DECREF(K);
        Py_DECREF(Phi);
        Py_DECREF(mu);
        Py_DECREF(supercell);
        Py_DECREF(cell);          
        PyErr_Format(PyExc_ValueError,
                        "Number of angles required!");
        return NULL;  
      }
  }

  
  /* Check optional arguments are present for calculation 'r' */
  if (icalc_type == 2) {
    if (!orot_axis) {
      Py_DECREF(positions);
      Py_DECREF(FC);
      Py_DECREF(K);
      Py_DECREF(Phi);
      Py_DECREF(mu);
      Py_DECREF(supercell);
      Py_DECREF(cell);          
      PyErr_Format(PyExc_ValueError,
                    "Axis for rotation required!");
      return NULL;        
    }
  }

  if (nnn > 200) {
    PyErr_Format(PyExc_RuntimeError,
                    "Error, number of nearest neighbours exceedingly large.");
    return NULL;
  }
  if (nangles < 1) {
    PyErr_Format(PyExc_RuntimeError,
                    "Error, number of angles must be greater than 0.");
    return NULL;
  }  
    
  
  /* Additional validation */
  if (r<=0.0) {
    PyErr_Warn( PyExc_UserWarning, "Radius for sums is <= 0!"); 
  }
  if (rcont<=0.0) {
    PyErr_Warn( PyExc_UserWarning, "Radius for contact hyperfine couplig is <= 0!"); 
  }


  /* Parse ptional argument for 'r' */
  if (icalc_type == 2) {
    rot_axis = (PyArrayObject *) PyArray_FROMANY(orot_axis, NPY_DOUBLE, 1, 1,
                                             NPY_ARRAY_IN_ARRAY);
    if (!rot_axis) {
      Py_DECREF(positions);
      Py_DECREF(FC);
      Py_DECREF(K);
      Py_DECREF(Phi);
      Py_DECREF(mu);
      Py_DECREF(supercell);
      Py_DECREF(cell);          
      PyErr_Format(PyExc_ValueError,
                    "Axis for rotation required but not parsed!");
      return NULL;        
    }
  }
  
  /* General validity checks */
  
  npy_intp pndims = PyArray_NDIM(positions);
  npy_intp fcndims = PyArray_NDIM(FC);

  npy_intp * pShape = PyArray_SHAPE(positions);
  npy_intp * fcShape = PyArray_SHAPE(FC);
  // The first condition is impossible since it is implied by the 
  //  PyArray_FROMANY function above. That function guaranties that the 
  //  comparison of the shape is fine.
  
  if (pndims != fcndims || 
        pShape[0]!=fcShape[0] ||
        pShape[1]!=fcShape[1]) {
    Py_DECREF(positions);
    Py_DECREF(FC);
    Py_DECREF(K);
    Py_DECREF(Phi);
    Py_DECREF(mu);
    Py_DECREF(supercell);
    Py_DECREF(cell);       
    Py_XDECREF(rot_axis);   // Null in case it is optional
    
    PyErr_SetString(PyExc_RuntimeError, "positions and FC arrays must have "
		    "same shape.");
    return NULL;
  }
  
  
  //npy_intp * phiSpahe = (npy_intp *) malloc(PyArray_NDIM(Phi) * sizeof(npy_intp));
  npy_intp * phiSpahe = PyArray_SHAPE(Phi);
  if (phiSpahe[0] != pShape[0]) {
    Py_DECREF(positions);
    Py_DECREF(FC);
    Py_DECREF(K);
    Py_DECREF(Phi);
    Py_DECREF(mu);
    Py_DECREF(supercell);
    Py_DECREF(cell);       
    Py_XDECREF(rot_axis);   // Null in case it is optional      
    
    PyErr_SetString(PyExc_RuntimeError, "positions and Phi arrays must have "
		    "same shape[0].");
    return NULL;
  }
  
  num_atoms = pShape[0]; 
  

  /* allocate variables needed for simulations */
  
  in_positions = (double *) malloc(3*num_atoms*sizeof(double));
  in_fc = (double *) malloc(6*num_atoms*sizeof(double));
  in_K = (double *) malloc(3*sizeof(double));
  in_phi = (double *) malloc(num_atoms*sizeof(double));
  in_muonpos = (double *) malloc(3*sizeof(double));
  in_supercell = malloc(3*sizeof(int));
  in_cell = (double *) malloc(9*sizeof(double));
  
  if (!in_positions || !in_fc || !in_K || !in_phi || !in_muonpos ||
      !in_supercell || !in_cell) {
    PyErr_SetString(PyExc_MemoryError, "Malloc failed.");
    Py_DECREF(positions);
    Py_DECREF(FC);
    Py_DECREF(K);
    Py_DECREF(Phi);
    Py_DECREF(mu);
    Py_DECREF(supercell);
    Py_DECREF(cell);       
    Py_XDECREF(rot_axis);   // Null in case it is optional      
    if (in_positions != NULL){
        free(in_positions);
    }    
    if (in_fc != NULL){
        free(in_fc);
    }    
    if (in_K != NULL){
        free(in_K);
    }    
    if (in_phi != NULL){
        free(in_phi);
    }    
    if (in_muonpos != NULL){
        free(in_muonpos);
    }    
    if (in_supercell != NULL){
        free(in_supercell);
    }    
    if (in_cell != NULL){
        free(in_cell);
    }
    return NULL;          
  }
  
  for (int i=0; i< num_atoms; i++){
    in_positions[3*i+0] = *(npy_float64 *) PyArray_GETPTR2(positions, i,0);
    in_positions[3*i+1] = *(npy_float64 *) PyArray_GETPTR2(positions, i,1);
    in_positions[3*i+2] = *(npy_float64 *) PyArray_GETPTR2(positions, i,2);
  }
  
  
  npy_cdouble v;
  for (int i=0; i< num_atoms; i++){
    v = *(npy_cdouble *)PyArray_GETPTR2(FC, i,0);
    in_fc[6*i+0] = v.real;
    in_fc[6*i+1] = v.imag;
    v = *(npy_cdouble *)PyArray_GETPTR2(FC, i,1);
    in_fc[6*i+2] = v.real;
    in_fc[6*i+3] = v.imag;    
    v = *(npy_cdouble *)PyArray_GETPTR2(FC, i,2);
    in_fc[6*i+4] = v.real;
    in_fc[6*i+5] = v.imag;
  }

  
  in_K[0] = *(npy_float64 *) PyArray_GETPTR1(K,0);
  in_K[1] = *(npy_float64 *) PyArray_GETPTR1(K,1);
  in_K[2] = *(npy_float64 *) PyArray_GETPTR1(K,2);
  
  
  for (int i=0; i< num_atoms; i++){
    in_phi[i] = *(npy_float64 *) PyArray_GETPTR1(Phi,i);
  }
  
  
  in_muonpos[0] = *(npy_float64 *) PyArray_GETPTR1(mu,0);
  in_muonpos[1] = *(npy_float64 *) PyArray_GETPTR1(mu,1);
  in_muonpos[2] = *(npy_float64 *) PyArray_GETPTR1(mu,2);
  
  
  in_supercell[0] = *(npy_int64 *)PyArray_GETPTR1(supercell, 0);
  in_supercell[1] = *(npy_int64 *)PyArray_GETPTR1(supercell, 1);
  in_supercell[2] = *(npy_int64 *)PyArray_GETPTR1(supercell, 2);
  
  
  in_cell[0] = *(npy_float64 *) PyArray_GETPTR2(cell,0,0);
  in_cell[1] = *(npy_float64 *) PyArray_GETPTR2(cell,0,1);
  in_cell[2] = *(npy_float64 *) PyArray_GETPTR2(cell,0,2);
  in_cell[3] = *(npy_float64 *) PyArray_GETPTR2(cell,1,0);
  in_cell[4] = *(npy_float64 *) PyArray_GETPTR2(cell,1,1);
  in_cell[5] = *(npy_float64 *) PyArray_GETPTR2(cell,1,2);
  in_cell[6] = *(npy_float64 *) PyArray_GETPTR2(cell,2,0);
  in_cell[7] = *(npy_float64 *) PyArray_GETPTR2(cell,2,1);
  in_cell[8] = *(npy_float64 *) PyArray_GETPTR2(cell,2,2);
  
  
  if (icalc_type == 2) {
    in_axis = (double *) malloc(3*sizeof(double));
    in_axis[0] = *(npy_float64 *)PyArray_GETPTR1(rot_axis, 0);
    in_axis[1] = *(npy_float64 *)PyArray_GETPTR1(rot_axis, 1);
    in_axis[2] = *(npy_float64 *)PyArray_GETPTR1(rot_axis, 2);
  }
  
  

  /* allocate output arrays */
  int nd = 1;
  if (icalc_type >= 2) {
    /* finish creating return variables */
    out_dim = (npy_intp *) malloc(2 * sizeof(npy_intp));
    out_dim[0] = (npy_intp) nangles;
    out_dim[1] = (npy_intp) 3;
    nd = 2;
  } else {
    out_dim = (npy_intp *) malloc(1 * sizeof(npy_intp));
    out_dim[0] = (npy_intp) 3;
    nd = 1;
  }
  
  odip = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);
  ocont = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);
  olor = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);
  if (out_dim != NULL) {
    free(out_dim); out_dim = NULL;
  }

  if (!odip || !ocont || !olor) {

    Py_XDECREF(odip);   
    Py_XDECREF(ocont);  
    Py_XDECREF(olor);
    if (in_positions != NULL)
      free(in_positions);
    if (in_fc != NULL)
      free(in_fc);
    if (in_K != NULL)
      free(in_K);
    if (in_phi != NULL)
      free(in_phi);
    if (in_muonpos != NULL)
      free(in_muonpos);
    if (in_supercell != NULL)
      free(in_supercell);
    if (in_cell != NULL)
      free(in_cell);          
    if (in_axis != NULL)
      free(in_axis);          
    PyErr_SetString(PyExc_MemoryError, "Cannot create output arrays.");
    return NULL;
  }


  dip = (double *) PyArray_DATA(odip);
  cont = (double *) PyArray_DATA(ocont);
  lor = (double *) PyArray_DATA(olor);
  
  // long computation starts here. No python object is touched so free thread execution
  Py_BEGIN_ALLOW_THREADS
  switch (icalc_type)
  {
    case 1:
      ASS(in_positions, in_fc, in_K, in_phi, in_muonpos, in_supercell, 
        in_cell,r, nnn,rcont,num_atoms,cont,dip,lor);
      break;
    case 2:
      RASS(in_positions, in_fc, in_K, in_phi, in_muonpos, in_supercell, 
        in_cell,r, nnn,rcont,num_atoms,in_axis,nangles,cont,dip,lor);
      break;
    case 3:
      INCASS(in_positions, in_fc, in_K, in_phi, in_muonpos, in_supercell, 
        in_cell,r, nnn,rcont,num_atoms,nangles,cont,dip,lor);
    
  }
  Py_END_ALLOW_THREADS

  if (in_positions != NULL){
      free(in_positions);
      in_positions = NULL;
  }
  if (in_fc != NULL)
    free(in_fc);
  if (in_K != NULL)
    free(in_K);
  if (in_phi != NULL)
    free(in_phi);
  if (in_muonpos != NULL)
    free(in_muonpos);
  if (in_supercell != NULL)
    free(in_supercell);
  if (in_cell != NULL)
    free(in_cell);          
  if (in_axis != NULL)            // Null in case it is optional     
    free(in_axis);  

  Py_DECREF(positions);
  Py_DECREF(FC);
  Py_DECREF(K);
  Py_DECREF(Phi);
  Py_DECREF(mu);
  Py_DECREF(supercell);
  Py_DECREF(cell);
  Py_XDECREF(rot_axis);   // Null in case it is optional      

  return Py_BuildValue("NNN", ocont, odip, olor);
}


static PyObject * py_lfcext_dt(PyObject *self, PyObject *args) {
    
  double r=0.0;
    
  PyObject *opositions, *omu, *osupercell, *ocell;
  
  PyArrayObject *positions,  *mu, *supercell, *cell;
  
  PyArrayObject *odt;
  int num_atoms=0;
  int * in_supercell;

  /* put arguments into variables */
  if (!PyArg_ParseTuple(args, "OOOOd", 
                            &opositions, &omu, &osupercell, &ocell,&r))
  {
    return NULL;
  }
  
  /* turn inputs into numpy array types */
  positions = (PyArrayObject *) PyArray_FROMANY(opositions, NPY_DOUBLE, 2, 2,
                                              NPY_ARRAY_IN_ARRAY);
  mu = (PyArrayObject *) PyArray_FROMANY(omu, NPY_DOUBLE, 1, 1,
                                              NPY_ARRAY_IN_ARRAY);
                                              
  supercell = (PyArrayObject *) PyArray_FROMANY(osupercell, NPY_INT32,
                                                   1, 1, NPY_ARRAY_IN_ARRAY);
  cell = (PyArrayObject *) PyArray_FROMANY(ocell, NPY_DOUBLE, 2, 2,
                                             NPY_ARRAY_IN_ARRAY);

  /* Validate data */
  if ( !positions || !mu || !supercell || !cell) {
    Py_XDECREF(positions);
    Py_XDECREF(mu);
    Py_XDECREF(supercell);
    Py_XDECREF(cell);
    PyErr_Format(PyExc_RuntimeError,
                    "Error parsing numpy arrays.");                 
    return NULL;
  }


  /* General validity checks */
  
  npy_intp * pShape = PyArray_SHAPE(positions);
  
  num_atoms = pShape[0];
  
  if (num_atoms == 0) {
    PyErr_Format(PyExc_RuntimeError,
                    "No magnetic atoms specified.");                 
    return NULL;      
  }


  /* to avoid problems with int types */
  in_supercell = malloc(3*sizeof(int));
  in_supercell[0] = *(npy_int64 *)PyArray_GETPTR1(supercell, 0);
  in_supercell[1] = *(npy_int64 *)PyArray_GETPTR1(supercell, 1);
  in_supercell[2] = *(npy_int64 *)PyArray_GETPTR1(supercell, 2);
  
 
  int nd = 2; 
  npy_intp * out_dim = (npy_intp *) malloc(2 * sizeof(npy_intp));
  out_dim[0] = (npy_intp) 3;
  out_dim[1] = (npy_intp) 3;
  
  odt = (PyArrayObject *) PyArray_ZEROS(nd, out_dim, NPY_DOUBLE,0);
  
  free(out_dim);

  if (!odt) {
    Py_XDECREF(positions);
    Py_XDECREF(mu);
    Py_XDECREF(supercell);
    Py_XDECREF(cell);
    if (in_supercell != NULL)
      free(in_supercell);
    PyErr_SetString(PyExc_MemoryError, "Cannot create output array.");
    return NULL;
  }
  
  
  DT( (double *) PyArray_DATA(positions),
      (double *) PyArray_DATA(mu),
      in_supercell, 
      (double *) PyArray_DATA(cell), 
      r, num_atoms, 
      (double *) PyArray_DATA(odt));

  Py_DECREF(positions);
  Py_DECREF(mu);
  Py_DECREF(supercell);
  Py_DECREF(cell);
  
  free(in_supercell);
  return Py_BuildValue("N", odt);

}

static PyMethodDef lfcext_methods[] =
{
  {"Fields", (PyCFunction)py_lfcext_fields, METH_VARARGS, py_lfcext_fields_docstring},
  {"DipolarTensor", (PyCFunction)py_lfcext_dt, METH_VARARGS, py_lfcext_dt_docstring},
  {NULL}  /* sentinel */
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "lfcext",            /* m_name */
  module_docstring,    /* m_doc */
  -1,                  /* m_size */
  lfcext_methods,      /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};
#endif


PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit_lfcext(void)
#else
initlfcext(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
  PyObject *module = PyModule_Create(&moduledef);
#else
  (void) Py_InitModule("lfcext", lfcext_methods);
#endif

  import_array(); /* Must be present for NumPy */

#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}
