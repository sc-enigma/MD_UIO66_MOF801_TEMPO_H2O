def get_mof801_params():
    # [atom_type] = m
    mass = {}
    mass['itcFF_Zr'] = 91.224
    mass['itcFF_O1'] = 16.000
    mass['itcFF_O2'] = 16.000
    mass['itcFF_C1'] = 12.011
    mass['itcFF_C2'] = 12.011
    mass['itcFF_H1'] =  1.008
    mass['itcFF_H7'] =  1.008
 
    # [atom_type] = q
    charge = {}
    charge['itcFF_Zr'] =  1.968
    charge['itcFF_O1'] = -0.533
    charge['itcFF_O2'] = -0.902
    charge['itcFF_C1'] =  0.630
    charge['itcFF_C2'] = -0.366
    charge['itcFF_H1'] =  0.131
    charge['itcFF_H7'] =  0.131
    
    # [atom_type-atom_type] = [funct, r0, k]
    bond_params = {}

    bond_params['itcFF_Zr-itcFF_O1'] = [1, 0.2232, 287290.200]
    bond_params['itcFF_Zr-itcFF_O2'] = [1, 0.2255, 107733.800]   
    bond_params['itcFF_C1-itcFF_O1'] = [1, 0.1273, 451872.000]
    bond_params['itcFF_O2-itcFF_H7'] = [1, 0.1080, 304106.800]
    bond_params['itcFF_C2-itcFF_H1'] = [1, 0.1080, 304106.800]
    bond_params['itcFF_C1-itcFF_C2'] = [1, 0.1269, 293928.300]
    bond_params['itcFF_C2-itcFF_C2'] = [1, 0.1523, 293928.300]
    
    # [atom_type-atom_type-atom_type] = [funct, angle, k]
    angle_params = {}
    angle_params['itcFF_O1-itcFF_Zr-itcFF_O1'] = [[1, 74.400, 115.776], [1, 120.000, 115.776]]
    angle_params['itcFF_O2-itcFF_Zr-itcFF_O2'] = [[1, 71.100, 115.776], [1, 100.400, 115.776]]
    angle_params['itcFF_O1-itcFF_Zr-itcFF_O2'] = [[1, 85.700, 115.776], [1, 124.300, 115.776], [1, 127.600, 115.776], [1, 143.500, 115.776]]
    angle_params['itcFF_C1-itcFF_O1-itcFF_Zr'] = [[1, 132.000, 231.637]]
    angle_params['itcFF_O1-itcFF_C1-itcFF_O1'] = [[1, 129.000,1213.360]]   
    angle_params['itcFF_C2-itcFF_C1-itcFF_O1'] = [[1, 117.300, 456.013]]   
    angle_params['itcFF_C2-itcFF_C2-itcFF_H1'] = [[1, 120.000, 309.616]]
    angle_params['itcFF_C1-itcFF_C2-itcFF_H1'] = [[1, 120.000, 309.616]]
    angle_params['itcFF_C1-itcFF_C2-itcFF_C2'] = [[1, 120.000, 231.637]]
    
    # [atom_type-atom_type-atom_type-atom_type] = [funct, angle, k, n]        - periodic
    # [atom_type-atom_type-atom_type-atom_type] = [funct, c1, c2, c3, c4, c5] - fourier
    dihedral_params = {}
    dihedral_params['itcFF_Zr-itcFF_O1-itcFF_C1-itcFF_C2'] = [9, 180.000, 86.837, 2]
    dihedral_params['itcFF_O1-itcFF_C1-itcFF_C2-itcFF_C2'] = [9, 180.000, 10.460, 2]
    dihedral_params['itcFF_O1-itcFF_C1-itcFF_C2-itcFF_H1'] = [9, 180.000, 12.552, 2]
    dihedral_params['itcFF_C1-itcFF_C2-itcFF_C2-itcFF_C1'] = [9, 180.000, 12.552, 2]
    dihedral_params['itcFF_C1-itcFF_C2-itcFF_C2-itcFF_H1'] = [9, 180.000, 12.552, 2]
    dihedral_params['itcFF_H1-itcFF_C2-itcFF_C2-itcFF_H1'] = [9, 180.000, 12.552, 2]

















    return mass, charge, bond_params, angle_params, dihedral_params