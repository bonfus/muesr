from muesr.io.cif.cif import parse_cif,parse_magn_operation_xyz_string
from muesr.io.cif.cell import cellpar_to_cell
from muesr.core.magmodel import MM

def mago_load_from_mcif(sample, fname):
    #http://www.physics.byu.edu/faculty/campbell/events/cmsworkshop2014/campbell_magCIF.pdf
    
    floateps = np.finfo(np.float).eps*10
    
    try:
        f = open(fname,'r')
    except:
        nprintmsg('efile')
        return
    
    data = parse_cif(f)
    
    f.close()
    
    nprint('WARNING! This function is VERY experimental! Do not use for production!','warn')
    
    nprint('Parsing data from: ' + data[0][0])
    
    values = data[0][1]
    
    # find magnetic cell
    
    a = values['_cell_length_a']
    b = values['_cell_length_b']
    c = values['_cell_length_c']
    alpha = values['_cell_angle_alpha']
    beta = values['_cell_angle_beta']
    gamma = values['_cell_angle_gamma']
    
    mag_cell = cellpar_to_cell([a, b, c, alpha, beta, gamma], (0,0,1), None)

    
    # Find magnetic atoms
    
    mag_atoms_labels = values['_atom_site_moment_label']

    # load mag moments
    mag_atoms_moments = np.zeros([len(mag_atoms_labels),3],dtype=np.float)
    
    for i in range(len(mag_atoms_labels)):
        
        mag_atoms_moments[i,0] = values['_atom_site_moment_crystalaxis_x'][i]
        mag_atoms_moments[i,1] = values['_atom_site_moment_crystalaxis_y'][i]
        mag_atoms_moments[i,2] = values['_atom_site_moment_crystalaxis_z'][i]
        
    # THESE ARE IN CRYSTAL AXIS COORDINATE SYSTEM!!!
    # bohr magneton units are used
    # the magnetic metric tensor is M = L.G.L^(-1), which is unitless. 
    # NOW GO TO REDUCED LATTICE COORDINATE SYSTEM TO DO THE SYMMETRY
    L = np.diag([1./a,1./b,1./c])    
    mag_atoms_moments = np.dot(mag_atoms_moments,L)
    
    mag_atoms_pos = np.zeros([len(mag_atoms_labels),3])
    mag_atoms_sym = []
    i = 0
    for label in mag_atoms_labels:
        j = values['_atom_site_label'].index(label)
        mag_atoms_pos[i,0] = values['_atom_site_fract_x'][j]
        mag_atoms_pos[i,1] = values['_atom_site_fract_y'][j]
        mag_atoms_pos[i,2] = values['_atom_site_fract_z'][j]
        mag_atoms_sym.append(values['_atom_site_type_symbol'][j])
        i += 1
    
    
    # Find transformation matrix between the two cells
    #  OCell = A * MCell
    
    unit_cell = sample.cell   
    fcs = np.zeros([unit_cell.get_number_of_atoms(),3],dtype=np.complex)
    
    atoms_types = [x.lower() for x in mag_atoms_sym] #lowercase
    set_this_fc = lambda x: x.lower() in atoms_types
    
    # count magnetic atoms and add their positions
    cart_pos = unit_cell.get_positions()
    ocell_cart_pos = []
    ocell_crys_pos = []
    for i, atom in enumerate(unit_cell):
        if set_this_fc(atom[0]): # check if we should set this fc
            ocell_cart_pos.append(cart_pos[i])
            ocell_crys_pos.append(atom[2])

    
    # get propagation vector
    k = np.fromstring(values['_magnetic_propagation_vector_kxkykz'][0],dtype=float, sep=',')
    
    # define coefficients for linear system
    ac = [[] for x in range(len(ocell_crys_pos))]
    bc = [[] for x in range(len(ocell_crys_pos))]
    
    
    for j, m_a_p in enumerate(mag_atoms_pos):
        for cent in values['_space_group_symop.magn_centering_xyz']:
            rc,tc,trc = parse_magn_operation_xyz_string(cent)
            cm_a_p = (np.dot(rc,m_a_p)+tc)
            
            for s in values['_space_group_symop.magn_operation_xyz']:
                
                r,t,tr = parse_magn_operation_xyz_string(s)
                
                symp = (np.dot(r,cm_a_p)+t)
                
                #print('Symp is: '+ str(symp))
                
                mag_atom_cart_pos = np.dot(symp, mag_cell)
                # clean noise
                mag_atom_cart_pos[np.abs(mag_atom_cart_pos)< floateps] = 0.
                
                mag_atom_crys_pos = np.dot(mag_atom_cart_pos,
                                np.linalg.inv(unit_cell.get_cell()))
                # clean noise
                mag_atom_crys_pos[np.abs(mag_atom_crys_pos)< floateps] = 0.                                        

                #print('trc,np.linalg.det(rc),tr,np.linalg.det(r),np.dot(r, np.dot(rc,mag_atoms_moments[j])), np.dot(rc,mag_atoms_moments[j])')
                
                #print(trc,np.linalg.det(rc),tr,np.linalg.det(r),np.dot(r, np.dot(rc,mag_atoms_moments[j])), np.dot(rc,mag_atoms_moments[j]))
                
                #print('r')
                #print(r)
                
                crysfc = trc*np.linalg.det(rc)*tr*np.linalg.det(r)*np.dot(r, np.dot(rc,mag_atoms_moments[j]))
                
                # Go to cartesian coordinates.
                cart_fcs = np.dot(crysfc,mag_cell)  ##### THIS IS WRONG!!
                
                for i, p in enumerate(ocell_crys_pos):
                    
                    spos = np.remainder(np.round(mag_atom_crys_pos,decimals=4),1.)


                    if np.allclose( spos , p, rtol=1e-03):
                        val = mag_atom_crys_pos - spos
                        cs = np.cos(2.*np.pi*np.dot(k,val))
                        sn = np.sin(2.*np.pi*np.dot(k,val))
                        
                        ac[i].append( [cs,sn])
                        bc[i].append( cart_fcs)
                        break
    
    
    j = 0
    for i, atom in enumerate(unit_cell):
        if set_this_fc(atom[0]): # check if we should set this fc
            minsqx = np.linalg.lstsq(np.array(ac[j]),np.array(bc[j])[:,0],rcond=1e-7)
            minsqy = np.linalg.lstsq(np.array(ac[j]),np.array(bc[j])[:,1],rcond=1e-7)
            minsqz = np.linalg.lstsq(np.array(ac[j]),np.array(bc[j])[:,2],rcond=1e-7)
            
            
            xre,xim = np.round(minsqx[0],decimals=6)
            yre,yim = np.round(minsqy[0],decimals=6)
            zre,zim = np.round(minsqz[0],decimals=6)
            fcs[i]=[np.complex(xre,xim),np.complex(yre,yim),np.complex(zre,zim)]
            #print ('fcs[%d]'%i)
            #print (fcs[i])
            j += 1
    

    nmm = MM(unit_cell.get_number_of_atoms(),unit_cell.get_cell())
    nmm.k = k
    nmm.fc_set(fcs)
    
    
    sample.mm = nmm


#    def _find_coefficients(self,ac,bc):
#        # try without phase
#        minsqx = np.linalg.lstsq(np.array(ac[j]),np.array(bc[j])[:,0],rcond=1e-7)
#        minsqy = np.linalg.lstsq(np.array(ac[j]),np.array(bc[j])[:,1],rcond=1e-7)
#        minsqz = np.linalg.lstsq(np.array(ac[j]),np.array(bc[j])[:,2],rcond=1e-7)
#        
#        if (minsqx[1] < 1e-10 and minsqy[1] < 1e-10 and minsqz[1] < 1e-10):
#        
#            xre,xim = np.round(minsqx[0],decimals=6)
#            yre,yim = np.round(minsqy[0],decimals=6)
#            zre,zim = np.round(minsqz[0],decimals=6)
#            
#        else:
#            raise RuntimeError
#            
#            
#        return [np.complex(xre,xim),np.complex(yre,yim),np.complex(zre,zim)]        
    
