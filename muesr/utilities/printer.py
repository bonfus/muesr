from muesr.core.nprint import print_cell as pc


        
def print_cell(sample):
    """
    Print base cell
    """
    if sample._check_lattice():
        pc(sample._cell)
        return True
            
            
        
   
