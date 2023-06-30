#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Abhinav Roy

Lab: MTD + McCue Research Group 
Department: Materials Science and Engineering
University: Northwestern University
"""
import numpy as np
import os
def labeling(Nx, Ny, Nz, num_atoms, composition, seed, generate_init_struct=False):
    NodeNext = []
    Node = []
    NodeS = []
    NodeNext = []
    coord = np.zeros((num_atoms,3),dtype=float)
    particles=0
    ctr1 = ctr2 = 0
    rng = np.random.default_rng(seed=seed)
    count_A = 0
    final_comp = 0
    for k in range(2*Nz):
        '''
        Code segment for randomly populating an FCC lattice and storing atom positions 
        for tracking the nearest neighbors.
        '''
        ctr1 += 1
        ctr2 = 0
        for i in range(2*Nx):
            ctr2 +=1
            for j in range(Ny):
                I = 0.5*float(i)
                if ((ctr1 % 2 == 0 and ctr2 % 2 == 0) or (ctr1 % 2 != 0 and ctr2 % 2 != 0)):
                    J = float(j) 
                else:
                    J = float(j) + 0.5
                K = 0.5*float(k) 
                coord[particles,0] = I
                coord[particles,1] = J
                coord[particles,2] = K
                particles += 1
                Node.append(particles)
    for i in range(num_atoms):
        if (rng.random() < composition):
            NodeS.append(1)
            count_A += 1

        else:
            NodeS.append(0)
    final_comp = count_A/num_atoms
    # The segment is executed if initial unlabeled structure needs to be visualized.
    if generate_init_struct == True:
        infilename = f"initial3D_{Nx}x{Ny}x{Nz}_{composition:0.2f}"
        with open(os.path.join("./output_xyz",infilename + ".xyz"), "w") as f:
            f.write("%d\n" % (num_atoms))
            f.write("%s\n" % "Type X Y Z Radius")                
            for i in range(num_atoms):
                f.write("%s %f %f %f %f\n" % (NodeS[i],coord[i,0],coord[i,1],
                                            coord[i,2],0.1))
                
    # Loop for generating the nearest neighbor list.   
    for m in range(num_atoms):
        I,J,K = coord[m,0],coord[m,1],coord[m,2]
        nn1 = nn2 = nn3 = nn4 =nn5 =  nn6 = nn7 = nn8 = nn9 = nn10 = nn11 = nn12 = 0
        # Corner 1
        if (I == 0.0 and J == 0.0 and K == 0.0):
            nn2 = nn3 = nn4 = nn6 = nn7 = nn8 = nn10 = nn11 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
        # Corner 2
        elif (I == float(Nx-0.5) and J == 0.0 and K == 0.0):
            nn2 = nn3 = nn4 = nn5 = nn7 = nn8 = nn9 = nn11 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]
        # Corner 3
        elif (I == 0.0 and J == float(Ny-0.5)  and K == 0.0):
            nn1 = nn2 = nn4 = nn6 = nn7 = nn8 = nn9 = nn10 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]

        # Corner 4
        elif (I == float(Nx-0.5) and J == float(Ny-0.5) and K == 0.0):
            nn1 = nn2 = nn4 = nn5 = nn7 = nn8 = nn9 = nn10 = nn11 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]        
            
        # Corner 5
        elif (I == 0.0 and J == 0.0 and K == float(Nz-0.5)):
            nn1 = nn3 = nn4 = nn5 = nn6 = nn8 = nn10 = nn11 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
                                
        # Corner 6
        elif (I == float(Nx-0.5) and J == 0.0 and K == float(Nz-0.5)):
            nn1 = nn3 = nn4 = nn5 = nn6 = nn7 = nn9 = nn11 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]
        # Corner 7
        elif (I == 0.0 and J == float(Nz-0.5) and K == float(Nz-0.5)):
            nn1 = nn2 = nn3 = nn5 = nn6 = nn8 = nn9 = nn10 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]
        # Corner 8
        elif (I == float(Nx-0.5) and J == float(Ny-0.5) and K == float(Nz-0.5)):
            nn1 = nn2 = nn3 = nn5 = nn6 = nn7 = nn9 = nn10 = nn11 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]
        # Edge 1
        elif ((I == 0.0 and K == 0.0)  and (J != 0.0 and J != float(Nz-0.5))):
            nn2 = nn4 = nn6 = nn7 = nn8 = nn10 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]
        # Edge 2
        elif ((I == float(Nx-0.5) and K == 0.0)  and (J != 0.0 and J != float(Nz-0.5))):
            nn2 = nn4 = nn5 = nn7 = nn8 = nn9 = nn11 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]

        # Edge 3
        elif ((J == 0.0 and K == 0.0)  and (I != 0.0 and I != float(Nx-0.5))):
            nn2 = nn3 = nn4 = nn7 = nn8 = nn11 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]

        # Edge 4
        elif ((J == float(Ny-0.5) and K == 0.0)  and (I != 0.0 and I != float(Nx-0.5))):
            nn1 = nn2 = nn4 = nn5 = nn6 = nn9 = nn10 = 0 
            for p in range(num_atoms):
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p] 

        # Edge 5
        elif ((I == 0.0 and K == float(Nz-0.5))  and (J != 0.0 and J != float(Nz-0.5))):
            nn1 = nn3 = nn5 = nn6 = nn8 = nn10 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]        

        # Edge 6
        elif ((I == float(Nx-0.5) and K == float(Nz-0.5))  and (J != 0.0 and J != float(Nz-0.5))):
            nn1 = nn3 = nn5 = nn6 = nn7 = nn9 = nn11 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]
        
        # Edge 7
        elif ((J == 0.0 and K == float(Nz-0.5))  and (I != 0.0 and I != float(Nx-0.5))):
            nn1 = nn3 = nn4 = nn5 = nn6 = nn11 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]

        # Edge 8
        elif ((J == float(Ny-0.5) and K == float(Nz-0.5))  and (I != 0.0 and I != float(Nx-0.5))):
            nn1 = nn2 = nn3 = nn5 = nn6 = nn9 = nn10 = 0 
            for p in range(num_atoms):
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]

        
        # Edge 9
        elif ((I == 0.0 and J == 0.0)  and (K != 0.0 and K != float(Nz-0.5))):
            nn3 = nn4 = nn6 = nn8 = nn10 = nn11 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]


        # Edge 10
        elif ((I == float(Nx-0.5) and J == 0.0)  and (K != 0.0 and K != float(Nz-0.5))):
            nn3 = nn4 = nn5 = nn7 = nn9 = nn11 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]
        
        # Edge 11
        elif ((I == 0.0 and J == float(Ny-0.5))  and (K != 0.0 and K != float(Nz-0.5))):
            nn1 = nn2 = nn6 = nn8 = nn9 = nn10 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]
        # Edge 12
        elif ((I == float(Nx-0.5) and J == float(Ny-0.5))  and (K != 0.0 and K != float(Nz-0.5))):
            nn1 = nn2 = nn5 = nn7 = nn9 = nn10 = nn11 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]

        # Face 1
        elif (J == 0.0 and (I != 0.0 and I != float(Nx-0.5))  and (K != 0.0 and K != float(Nz-0.5))):
            nn3 = nn4 = nn11 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]

        # Face 2
        elif (J == float(Ny-0.5) and (I != 0.0 and I != float(Nx-0.5))  and (K != 0.0 and K != float(Nz-0.5))):
            nn1 = nn2 = nn9 = nn10 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]


        # Face 3
        elif (K == 0.0 and (I != 0.0 and I != float(Nx-0.5))  and (J != 0.0 and J != float(Ny-0.5))):
            nn2 = nn4 = nn7 = nn8 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]

        # Face 4
        elif (K == float(Nz-0.5) and (I != 0.0 and I != float(Nx-0.5))  and (J != 0.0 and J != float(Ny-0.5))):
            nn1 = nn3 = nn5 = nn6 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]

        # Face 5
        elif (I == 0.0 and (J != 0.0 and J != float(Ny-0.5))  and (K != 0.0 and K != float(Nz-0.5))):
            nn6 = nn8 = nn10 = nn12 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]

        # Face 6
        elif (I == float(Nx-0.5) and (J != 0.0 and J != float(Ny-0.5))  and (K != 0.0 and K != float(Nz-0.5))):
            nn5 = nn7 = nn9 = nn11 = 0
            for p in range(num_atoms):
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]

        #if ((I != 0.0 and I != float(Nx-0.5)) and (J != 0.0 and J != float(Ny-0.5)) and (K != 0.0 and K != float(Nz-0.5))):
        else:
            for p in range(num_atoms):
                # Keeping all the I index constant
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn1 = Node[p]
                if (((coord[p,1] == J+0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn2 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K+0.5)) and coord[p,0] == I):
                    nn3 = Node[p]
                if (((coord[p,1] == J-0.5) and (coord[p,2] == K-0.5)) and coord[p,0] == I):
                    nn4 = Node[p]
                # Keeping all the J index constant
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn5 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K+0.5)) and coord[p,1] == J):
                    nn6 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn7 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,2] == K-0.5)) and coord[p,1] == J):
                    nn8 = Node[p]
                # Keeping all the K index constant
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn9 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J+0.5)) and coord[p,2] == K):
                    nn10 = Node[p]
                if (((coord[p,0] == I+0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn11 = Node[p]
                if (((coord[p,0] == I-0.5) and (coord[p,1] == J-0.5)) and coord[p,2] == K):
                    nn12 = Node[p]
        nnlist = [nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8,nn9,nn10,nn11,nn12]
        NodeNext.append(nnlist)

    # Cluster labelling segment of the code            
    NodeL = [0]*len(Node)
    NodeLP = []
    clusters = 0
    for i in range(num_atoms):
        if (NodeS[i] != 0):
            idx = [j for j, e in enumerate(NodeNext[i]) if e != 0]
            nlst = []
            for j in idx:
                if (NodeS[NodeNext[i][j]-1] != 0):
                    nlst.append(j)
            if (len(nlst) != 0):
                Nodes = [NodeNext[i][j] for j in nlst]
                NodeNextL = [NodeL[j-1] for j in Nodes]        
                if (sum(NodeNextL) == 0):
                    clusters += 1
                    NodeL[i] = clusters
                    NodeLP.append(clusters)
                else:
                    M = [e for e in NodeNextL if e != 0]       
                    N = []
                    for j in M:
                        N.append(NodeLP[j-1])
                    lst = [NodeLP[j-1] for j in N]   
                    NodeLP_min = min(lst)
                    NodeL[i] = NodeLP_min
                    for j in N:
                        NodeLP[j-1] = NodeLP_min
        
            else:
                clusters += 1
                NodeL[i] = clusters
                NodeLP.append(clusters)

    for i in range(1,len(NodeLP)+1):
        N=i
        while (NodeLP[N-1]<N):
            N = NodeLP[N-1]
        NodeLP[i-1] = N
    NodeLP1 = []
    NodeLP1 = NodeLP.copy()
    NodeLP1.sort()
    lst = [a > b for a,b in zip(NodeLP1[1:],NodeLP1[0:-1])]
    relabl = [a*b for a,b in zip(NodeLP1[1:],lst)]
    lst2 = [e for e in relabl if e != 0]
    relabl = []
    relabl.append(NodeLP1[0])
    for j in lst2:
        relabl.append(j)
    for i in range(1,len(relabl)+1):
        for j,e in enumerate(NodeLP):
            if (e == relabl[i-1]):
                NodeLP[j] = i   
    for i in range(1,len(NodeLP)+1):
        for j,e in enumerate(NodeL):
            if (i == e):
                NodeL[j] = NodeLP[i-1]
    # Return the coordinates of atoms and the assigned cluster labels.
    return NodeL, coord, final_comp

# To be executed for generating labelled FCC structures for visualization purpose.       
if __name__ == '__main__':
    # To record the code execution time
    from time import perf_counter
    import sys
    t_start = perf_counter()
    # Define the variables for which the labelled structure needs to be generated. 
    comp = 0.15
    sys_size = int(sys.argv[1])
    seed = int(sys.argv[2])
    Nx = Ny = Nz = sys_size
    num_atoms = 4*Nx*Ny*Nz
    # Variables for storing the coordinates of atoms and the assigned labels.
    coord = np.zeros((num_atoms,3),dtype=float)
    NodeL = []
    NodeL, coord, final_comp = labeling(Nx, Ny, Nz, num_atoms, comp, seed, generate_init_struct=True)
    # Storing the final structure in a XYZ file for visualization.
    outfilename = f"final3D_{Nx}x{Ny}x{Nz}_{comp:0.2f}" 
    with open(os.path.join("./output_xyz",outfilename + ".xyz"), "w") as f:
        f.write("%d\n" % (num_atoms))
        f.write("%s\n" % "Type X Y Z Radius")                
        for i in range(num_atoms):
            f.write("%s %f %f %f %f\n" % (NodeL[i],coord[i,0],coord[i,1],
                                        coord[i,2],0.1))
    # Report the total time for code execution.
    t_stop = perf_counter()
    print("Time elapsed: %.4f sec" % (t_stop-t_start))
        
