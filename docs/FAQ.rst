Frequently Asked Questions
==========================

Symmetry is not correctly recognized by spglib!
  Try lowering the numerical precision threshold by running, 
  for example, ::

      >>> symsearch(yoursample, precision=1e-3)

  If that does not work I'm keen to think that there is an error in 
  your input structure since spglib is a well tested piece of code.
  If your input is correct file a bug at `spglib.sf.net <http://spglib.sf.net>`_.


How do I convert the hyperfine filed to reasonable units?
  See ContactTerm
