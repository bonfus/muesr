from muesr.core.nprint import print_cell


  
        
def print_cell(sample):
    """
    Print base cell
    """
    if sample._check_lattice():
        print_cell(sample._cell)
            
            
        
   
