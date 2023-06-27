#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Abhinav Roy

Lab: MTD + McCue Research Group 
Department: Materials Science and Engineering
University: Northwestern University
"""
#--------------------------------------------------------------------------------------------------
def threshold(atom_index, coord, num_atoms, Nx, Ny, Nz):
    def common_elements(a,b):
        common=0
        for value in a:
            if value in b:
                common=1
        return common
    face_1 = []
    face_2 = []
    face_3 = []
    face_4 = []
    face_5 = []
    face_6 = []
    for m in range(num_atoms): 
        I,J,K = coord[m,0],coord[m,1],coord[m,2]
        if (atom_index[m] != 0):
            # Face 1  
            if (J == 0.0 and I >= 0.0  and K >= 0.0):
                face_1.append(atom_index[m])

            # Face 2
            if (J == float(Ny-0.5) and I >= 0.0  and K >= 0.0):
                face_2.append(atom_index[m])

            # Face 3
            if (K == 0.0 and I >= 0.0 and J >= 0):
                face_3.append(atom_index[m])

            # Face 4
            if (K == float(Nz-0.5) and I >= 0.0 and J >= 0):
                face_4.append(atom_index[m])

            # Face 5
            if (I == 0.0 and K >= 0.0 and J >= 0):
                face_5.append(atom_index[m])

            # Face 6
            if (I == float(Nx-0.5) and K >= 0.0 and J >= 0):
                face_6.append(atom_index[m])
    Ry = common_elements(face_1, face_2)
    Rx = common_elements(face_5, face_6)
    Rz = common_elements(face_3, face_4)

    return Rx, Ry, Rz


    
