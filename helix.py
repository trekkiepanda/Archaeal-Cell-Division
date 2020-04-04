from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import math
import time
# from scipy.spatial import KDTree
import sys, getopt
# from numpy import cross, eye, dot
# from scipy.linalg import expm3, norm
# from scipy.optimize import minimize
# from mpl_toolkits.mplot3d import Axes3D

'''
This program generates proteins made out of 4 beads that attach to each other via 2 patches
and places them in a spiral position. It generates the Lammps read-in file.
'''

def proteingenerator(radius_o,x_start,x_step,n,cell_radius,tau):
    # Generate unrotated protein positions
    P1 = [radius_o*np.cos(tau)-x_start+n*x_step, 0 ,-cell_radius-np.sin(tau)*radius_o]
    P2 = [-radius_o*np.cos(tau)-x_start+n*x_step, 0 ,-cell_radius+np.sin(tau)*radius_o]
    P3 = [-x_start+n*x_step+np.sin(tau)*radius_o*np.sqrt(3), 0 ,-cell_radius+np.sqrt(3)*radius_o*np.cos(tau)]
    return P1,P2,P3


def plane_rotation(theta,coordinates):
    x = coordinates[0]
    y = np.cos(theta)*coordinates[1]-np.sin(theta)*coordinates[2]
    z = np.sin(theta)*coordinates[1]+np.cos(theta)*coordinates[2]
    return x,y,z

def Proteinrotator(P1,P2,P3,angle):
    P1_r = plane_rotation(angle,P1)
    P2_r = plane_rotation(angle,P2)
    P3_r = plane_rotation(angle,P3)
    protein = [P1_r,P2_r,P3_r]
    return protein

def bondcalculator(protein,protein_relaxed):
    res = np.zeros([3,3])
    for i, p in enumerate(protein):
        for j, pr in enumerate(protein_relaxed):
            res[i,j] = np.sqrt(np.sum((np.array(p)-pr)**2))
#     print('bond1 : '+str(res[0,0])+'\n')
#     print('bond2 : '+str(res[0,1])+'\n')
#     print('bond3 : '+str(res[0,2])+'\n')
#     print('bond4 : '+str(res[1,0])+'\n')
#     print('bond5 : '+str(res[1,1])+'\n')
#     print('bond6 : '+str(res[1,2])+'\n')
#     print('bond7 : '+str(res[2,0])+'\n')
#     print('bond8 : '+str(res[2,1])+'\n')
#     print('bond9 : '+str(res[2,2])+'\n')
    return res

def packingcalc(sigma,packing): # calculates packing of the cytoplasm
    Sigma = int(math.ceil(sigma))
    if Sigma < 4:
        print("Please change Sigma Value!")
    else:
        sig = np.arange(4,20,1)
        max = [0.56638, 0.5273, 0.4912, 0.435, 0.4296, 0.453, 0.3611, 0.3334, 0.3047, 0.2698, 0.2536,0.2339,0.2092,0.1971,0.1773,0.1584,0.1362]
        for i, item in enumerate(sig):
            if item == Sigma:
                index = i
        max_packing = max[index]
    if (packing > max_packing):
        packing = max_packing
        print("Packing is greater than max! Packing has been set = {0:1.1f}% \n".format(max_packing*100))
        return sigma, packing
    else:
        return sigma, packing


def cytoplasm(packing, cell_radius, actual_cyto_radius):
    x_c = []
    y_c = []
    z_c = []
    x_cyto = []
    y_cyto = []
    z_cyto = []
    s = 600                                                  # length of cytoplasm 'box'
    num_per_row = int(s/(2*actual_cyto_radius))              # number of cyto particles per row
    # initialising a box of cytoplasm particles in centre of cell
    x_adjust = 0                                             # adjustments, making sure the cyto 'box' is in the middle of the cell
    y_adjust = 11 * actual_cyto_radius
    z_adjust = 20 * actual_cyto_radius
    end_of_box = -0.5*s + 0.5 *actual_cyto_radius
    for k in range(num_per_row):
        z_val = end_of_box + z_adjust + actual_cyto_radius * (2 * np.sqrt(6) * k)/3
        for i in range (num_per_row):
            if (i%2 == 0):
                for j in range (num_per_row):
                    x_c.append(end_of_box + x_adjust+ 2 * actual_cyto_radius * j)
                    y_c.append(end_of_box + y_adjust+ np.sqrt(3) * actual_cyto_radius * i)
                    z_c.append(z_val)
            if (i%2 == 1):
                for j in range (num_per_row):
                    x_c.append(end_of_box  + x_adjust+ 2 * actual_cyto_radius * j + actual_cyto_radius)
                    y_c.append(end_of_box  + y_adjust+ np.sqrt(3) * actual_cyto_radius * i)
                    z_c.append(z_val)

    min_radius = cell_radius -  2 * actual_cyto_radius               # making sure cyto doesn't touch membrane
    V_cell = (4/3)*np.pi*cell_radius**3                                # volume of whole cell
    V_cyto = (4/3)*np.pi*actual_cyto_radius**3                                # volume per cyto particle
    num_cyto = packing*V_cell/V_cyto
    m = 0
    n = 0
    while n < num_cyto and m < len(x_c):
        r = np.sqrt (x_c[m] **2 + y_c[m] **2 + z_c[m] **2)
        if r < min_radius:
            x_cyto.append(x_c[m])
            y_cyto.append(y_c[m])
            z_cyto.append(z_c[m])
            n += 1
        m += 1
    actual_packing = len(x_cyto) * V_cyto /V_cell                      # packing of the cytoplasm in the cell
    print("cytoplasm packing =  {0:1.1f}% \n".format(actual_packing*100))
    return x_cyto, y_cyto, z_cyto

def main():

    ################ Parameters ################
    seed = int(sys.argv[6]) # seed
    # filament is made of (particle tripletts = proteins)
    radius = 0.5 						# particle radius
    dr = 0.2     						# overlap length of particles
    radius_o = radius-(dr/2)  					# distance between center of one particle to the middle between two particles
    # each triplett/protein is connected by 9 bonds into filaments
    bondstrength = 600						# bond strength of filament								*******VARY*******
    bondstrength2 = 0						# bond strength is set to zero to account for disassembly events
    # the preferred curvature of a filament is defined by alpha0, a0, and R
    alpha0_input = float(sys.argv[1])
    alpha0 = alpha0_input*np.pi/180				# preferred inner curvature (1-alpha0 is actual curvature) of filament			*******VARY*******
    a0 = 0.37  							# distance between centers of particles is dconst=a0+2radius				*******VARY*******

    R = (radius+a0/2)/np.sin(alpha0/2)				# target ring radius
    print((1-alpha0)/R)
    xstep = 0.01 						# shift distance of helix in cell (pitch)						*******VARY*******
    tau = 0 						# target ring twist
    # simulation parameters
    np.random.seed(421) 					# random seed for permutations
    lbox=math.fabs(100)
    thermo_relax = 1000						# thermo for "run: relax"
    run_relax = 10000						# time duration of "run: relax" -> leads to onset of "run: initial"
    thermo_initial = 1000					# thermo for "run: initial"
    run_initial = 200000					# time duration of "run: initial" -> leads to onset of "run: bonds transform"
    thermo_transform = 1000					# thermo for "run: bonds transform"
    run_transform = int(sys.argv[2])				# time duration of "run: bonds transform" -> leads to onset of "run: bonds disassemble"
    thermo_disassemble = 1000					# thermo for "run: bonds disassemble"
    run_disassemble = int(sys.argv[3])				# time duration of a single "run: bonds disassemble" = time it takes to disassemble a single protein at each end of the filament
    thermo_final = 1000						# thermo for "run: final"
    run_final = 200000						# time duration of "run: final"
    ################ Spherical Membrane ################
    fname_in = 'spherecoords_48002.xyz' 			# read in positions of membrane particles
    f_in = open(fname_in,'rU')
    tmp = f_in.readline()
    tmp = f_in.readline()					#? two times necessary?
    x_membrane = []
    y_membrane = []
    z_membrane = []
    r_membrane = []
    while True:
        tmp = f_in.readline()
        if tmp:
            x_membrane.append(float(tmp.split()[1]))
            y_membrane.append(float(tmp.split()[2]))
            z_membrane.append(float(tmp.split()[3]))
            r_membrane.append(np.sqrt(float(tmp.split()[1])**2+float(tmp.split()[2])**2+float(tmp.split()[2])**2))
        else:							#? shouldn't the last "tmp.split()[2]" be "tmp.split()[3]"?
            break

    cell_radius = np.mean(r_membrane)-1 				# radius of the cell (-1 such that the generated filament definitely is inside the cell)
    print('cell_radius = ' + str(cell_radius) + '\n')
    Nmembrane = len(x_membrane)					# number of membrane particles

    ################    removing the membrane particles within so cytoplasm doesn't overlap  #################
    x_mem = []
    y_mem = []
    z_mem = []
    r_mem = []
    for k in range(Nmembrane):
        r = np.sqrt(x_membrane[k]**2 + y_membrane[k]**2 + z_membrane[k]**2)
        if (r > cell_radius and r < cell_radius + 10):                    # only membrane particles that are on the cell
            x_mem.append(x_membrane[k])
            y_mem.append(y_membrane[k])
            z_mem.append(z_membrane[k])
            r_mem.append(r_membrane[k])

    ################ Generate helical positions ################
    anglestep = 2*np.arcsin((2*radius+a0)/(2*cell_radius))
    bead_per_circ = int(np.floor(2*np.pi/anglestep))		# number of proteins per circulation of cell (a0 is always fixed)
    Ncircles = 2						# number of circulations
    Nprotein = Ncircles*bead_per_circ				# number of proteins per filament

    xstart = Nprotein/2*xstep					#? is it: xstart = Nprotein/(Ncircles*xstep)?
    spiral = []

    for n in range(Nprotein):
        P1,P2,P3 = proteingenerator(radius_o,xstart,xstep,n,cell_radius,tau)
        protein = [P1,P2,P3]
        protein_rot = Proteinrotator(P1,P2,P3,n*anglestep)
        spiral.extend(protein_rot[:3])

    ################ Bond distances ################
    # stiff filament
    [P1,P2,P3] = proteingenerator(radius_o,0,0,0,cell_radius,tau)
    protein = [P1,P2,P3]
    protein_rot = Proteinrotator(P1,P2,P3,anglestep)
    res_s = bondcalculator(protein,protein_rot)

    # elastic filament
    [P1,P2,P3] = proteingenerator(radius_o,0,0,0,R,tau)
    protein = [P1,P2,P3]
    protein_rot = Proteinrotator(P1,P2,P3,alpha0)
    res_e =  bondcalculator(protein,protein_rot)

    ############### cytoplasm particles #################
    sigma_cyto =   int(sys.argv[4])                                              # sigma value of cytoplasm, *VARY* IDEALLY BET 4-15
    desired_packing = int(sys.argv[5])  /100                                         # packing fraction od cyto, *VARY* ideally max 30%
    sigma_cyto, packing = packingcalc(sigma_cyto, desired_packing)               # calculations for packing
    cyto_radius = sigma_cyto  * 1.12/2                                  # radius of cytoplasm particles

    x_cyto, y_cyto, z_cyto = cytoplasm (packing,cell_radius, cyto_radius)      # generates cytoplasm

    ################ Number of particles ################
    Nmem = len(x_mem)					# number of membrane particles
    Ncyto = len(x_cyto)                 # number of cytoplasm particles
    Nspiral = Nprotein*3					# number of particles in filament
    Nbonds = (Nprotein-1)*9					# number of bonds in filament
    Natoms = Nmem+Nspiral +Ncyto					# total number of particles="atoms"
    Nparticletypes = Nspiral+1		+1			# total number of particle types = number of particles in filament + 1 for membrane and cytoplasm
    print('number of membrane particles = ' + str(Nmem) + '\n')
    print('number of filament proteins = ' + str(Nprotein) + '\n')
    print('number of filament particles = ' + str(Nspiral) + '\n')
    print('number of cytoplasm particles = ' + str(Ncyto) + '\n')
    #x_coords = [s[0] for s in spiral]
    #y_coords = [s[1] for s in spiral]
    #plt.scatter(x_coords,y_coords)
    #plt.show()

    ################  Write polyspiral.in ################
    print('writing polyspiral.in\n')
    outfile='polyspiral.in'
    f=open(outfile,'w')
    f.write('LAMMPS data file generated by tc387 with a python script for hex membrane\n')
    f.write('\n')
    f.write(str(Natoms) + '   atoms\n')
    f.write(str(Nbonds) + '   bonds\n')
    f.write('\n')
    f.write(str(Nparticletypes)+' atom types\n')
    f.write(str(Nbonds)+' bond types\n')
    f.write('\n')
    f.write(str(-lbox)+' '+str(lbox)+' xlo xhi \n')
    f.write(str(-lbox)+' '+str(lbox)+' ylo yhi \n')
    f.write(str(-lbox)+' '+str(lbox)+' zlo zhi \n')
    f.write('\n')
    f.write('Masses\n')
    f.write('\n')
    for ii in range(Nparticletypes):
        f.write(str(ii+1) + ' 1\n')
    f.write('\n')
    f.write('Atoms # hybrid\n')
    f.write('\n')

    for ii in range(Nmem):
        # membrane particles
        r = np.sqrt(x_mem[ii]**2+y_mem[ii]**2+z_mem[ii]**2)
        f.write(str(ii+1) + '\t' + '1' + '\t' + str(x_mem[ii]) + '\t' + str(y_mem[ii]) + '\t' + str(z_mem[ii]) + '\t1\t1\t0\t' +   str(x_mem[ii]/r) + '\t' + str(y_mem[ii]/r) + '\t' + str(z_mem[ii]/r) + '\t' + str(ii+1) + '\n')
    for ii in range(Ncyto):
        # cytoplasm particles
        f.write(str(Nmem+ii+1) + '\t' + '2' + '\t' + str(x_cyto[ii]) + '\t' + str(y_cyto[ii]) + '\t' + str(z_cyto[ii]) + '\t1\t1\t0\t0\t0\t0\t' + str(ii+1) + '\n')
    for ii in range(Nprotein):
        # inner particles
        f.write(str(1 +Nmem +Ncyto +ii*3) + '\t' + str(3+3*ii) + '\t' + str(spiral[3*ii][0]) + '\t' + str(spiral[3*ii][1]) + '\t' + str(spiral[3*ii][2]) + '\t1\t1\t0\t0\t0\t0\t' + str(ii+1) + '\n')
        # outer particles
        f.write(str(2 +Nmem +Ncyto +ii*3) + '\t' + str(4+3*ii) + '\t' + str(spiral[3*ii+1][0]) + '\t' + str(spiral[3*ii+1][1]) + '\t' + str(spiral[3*ii+1][2]) + '\t1\t1\t0\t0\t0\t0\t' + str(ii+1) + '\n')
        # upper particles
        f.write(str(3 +Nmem +Ncyto +ii*3) + '\t' + str(5+3*ii) + '\t' + str(spiral[3*ii+2][0]) + '\t' + str(spiral[3*ii+2][1]) + '\t' + str(spiral[3*ii+2][2]) + '\t1\t1\t0\t0\t0\t0\t' + str(ii+1) + '\n')
    f.write('\n')
    f.write('Bonds \n')
    f.write('\n')
    for ii in range(Nprotein-1):
        # bond inner inner
        f.write(str(9*ii+1) + '\t' + str(9*ii+1) + '\t' +  str(1+Nmem +Ncyto +ii*3) + '\t' + str(1+Nmem +Ncyto +(ii+1)*3) + '\n')
        # bond upper outer
        f.write(str(9*ii+2) + '\t' + str(9*ii+2) + '\t' +   str(1+Nmem +Ncyto +ii*3) + '\t' + str(2+Nmem +Ncyto +(ii+1)*3)  + '\n')
         # bond upper inner
        f.write(str(9*ii+3) + '\t' + str(9*ii+3) + '\t' +   str(1+Nmem +Ncyto +ii*3)  + '\t' + str(3+Nmem +Ncyto +(ii+1)*3)  + '\n')
        # bond outer upper
        f.write(str(9*ii+4) + '\t' + str(9*ii+4) + '\t' +  str(2+Nmem +Ncyto +ii*3) + '\t' + str(1+Nmem +Ncyto +(ii+1)*3) + '\n')
        # bond outer outer
        f.write(str(9*ii+5) + '\t' + str(9*ii+5) + '\t' +  str(2+Nmem +Ncyto +ii*3) + '\t' + str(2+Nmem +Ncyto +(ii+1)*3) + '\n')
        # bond outer inner
        f.write(str(9*ii+6) + '\t' + str(9*ii+6) + '\t' +  str(2+Nmem +Ncyto +ii*3) + '\t' + str(3+Nmem +Ncyto +(ii+1)*3) + '\n')
        # bond inner upper
        f.write(str(9*ii+7) + '\t' + str(9*ii+7) + '\t' +  str(3+Nmem +Ncyto +ii*3) + '\t' + str(1+Nmem +Ncyto +(ii+1)*3) + '\n')
        # bond inner outer
        f.write(str(9*ii+8) + '\t' + str(9*ii+8) + '\t' + str(3+Nmem +Ncyto +ii*3) + '\t' + str(2+Nmem +Ncyto +(ii+1)*3) + '\n')
        # bond inner inner
        f.write(str(9*ii+9) + '\t' + str(9*ii+9) + '\t' +  str(3+Nmem +Ncyto +ii*3) + '\t' + str(3+Nmem +Ncyto +(ii+1)*3) + '\n')

    ################ write in.local ################
    print('writing in.local\n')
    infile = 'in.local'
    with open(infile,'rU') as fi:
        data = fi.readlines()
    outfile2='in.local'
    fo=open(outfile2,'w')
    for line in data:

        if len(line.strip().split())>0 and line.strip().split()[0] == 'group' and line.strip().split()[1] == 'spiral':
            fo.write('group\t\tspiral\t\ttype 3:' + str(Nparticletypes) + '\t\t# rigid polymer \n')

        elif len(line.strip().split())>0 and line.strip().split()[0] == '#pw359:':
            fo.write('#pw359: Initialise LJ/expand to zero for all possible combinations \n')
            fo.write('pair_coeff\t*\t*\tlj/cut\t0\t0\n')
            for ii in range(Nprotein):
                fo.write('pair_coeff\t1\t' + str(3*ii+3) + '\tlj/cut\t${eps_ms}\t${sigma_ms}\t${rc_ms}\n')		# attraction membrane and filament bead 3
                fo.write('pair_coeff\t1\t' + str(3*ii+4) + '\tlj/cut\t${eps_ms}\t${sigma_ms}\t${rc_ms}\n')		# attraction membrane and filament bead 4
                fo.write('pair_coeff\t1\t' + str(3*ii+5) + '\tlj/cut\t${eps_mu}\t${sigma_mu}\t${rc_mu}\n')		# volume exclusion membrane and filament bead 5
            for ii in range(Nprotein):
                fo.write('pair_coeff\t2\t' + str(3*ii+3) + '\tlj/cut\t${eps_cyto}\t${sigma_cyto}\t${rc_cyto}\n')		# volume exclusion and filament bead 3
                fo.write('pair_coeff\t2\t' + str(3*ii+4) + '\tlj/cut\t${eps_cyto}\t${sigma_cyto}\t${rc_cyto}\n')		# volume exclusion and filament bead 4
                fo.write('pair_coeff\t2\t' + str(3*ii+5) + '\tlj/cut\t${eps_cyto}\t${sigma_cyto}\t${rc_cyto}\n')		# volume exclusion cytoplasm  and filament bead 5
            fo.write('pair_coeff\t3*\t3*\tlj/cut\t${eps_ss}\t${sigma_ss}\t${rc_ss}\n')					# volume exclusion among filament beads
            fo.write('pair_coeff\t1\t1\tmembrane\t${eps}\t${sigma}\t${rmin}\t${rc}\t${zeta}\t${mu}\t${theta0_11}\n') 	# volume exclusion among membrane beads
            fo.write('pair_coeff\t1\t2\tlj/cut\t${eps_cyto}\t${sigma_cyto}\t${rc_cyto}\n') 	# volume exclusion among membrane and cytoplasm
            fo.write('pair_coeff\t2\t2\tlj/cut\t${eps_cc}\t${sigma_cc}\t${rc_cc}\n') 	# volume exclusion among cytoplasm

        elif len(line.strip().split())>0 and line.strip().split()[0] == '#bonds_initial':
            fo.write('#bonds_initial\n')
            for ii in range(Nprotein-1):
                fo.write('bond_coeff\t' + str(9*ii+1) + '\t' + str(bondstrength) + '\t' + str(res_s[0,0]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+2) + '\t' + str(bondstrength) + '\t' + str(res_s[0,1]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+3) + '\t' + str(bondstrength) + '\t' + str(res_s[0,2]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+4) + '\t' + str(bondstrength) + '\t' + str(res_s[1,0]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+5) + '\t' + str(bondstrength) + '\t' + str(res_s[1,1]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+6) + '\t' + str(bondstrength) + '\t' + str(res_s[1,2]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+7) + '\t' + str(bondstrength) + '\t' + str(res_s[2,0]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+8) + '\t' + str(bondstrength) + '\t' + str(res_s[2,1]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+9) + '\t' + str(bondstrength) + '\t' + str(res_s[2,2]) + '\n' )

        elif len(line.strip().split())>0 and line.strip().split()[0] == '#run_relax':
            fo.write('#run_relax\n')
            fo.write('timestep\t0.01\n')
            fo.write('thermo\t\t' + str(thermo_relax) + '\n')
            fo.write('run\t\t' + str(run_relax) + '\n')

        elif len(line.strip().split())>0 and line.strip().split()[0] == '#run_initial':
            fo.write('#run_initial\n')
            fo.write('timestep\t0.01\n')
            fo.write('thermo\t\t' + str(thermo_initial) + '\n')
            fo.write('run\t\t' + str(run_initial) + '\n')

        elif len(line.strip().split())>0 and line.strip().split()[0] == '#run_bonds_transform':
            fo.write('#run_bonds_transform\n')
            for ii in range(Nprotein-1):
                fo.write('bond_coeff\t' + str(9*ii+1) + '\t' + str(bondstrength) + '\t' + str(res_e[0,0]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+2) + '\t' + str(bondstrength) + '\t' + str(res_e[0,1]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+3) + '\t' + str(bondstrength) + '\t' + str(res_e[0,2]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+4) + '\t' + str(bondstrength) + '\t' + str(res_e[1,0]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+5) + '\t' + str(bondstrength) + '\t' + str(res_e[1,1]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+6) + '\t' + str(bondstrength) + '\t' + str(res_e[1,2]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+7) + '\t' + str(bondstrength) + '\t' + str(res_e[2,0]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+8) + '\t' + str(bondstrength) + '\t' + str(res_e[2,1]) + '\n' )
                fo.write('bond_coeff\t' + str(9*ii+9) + '\t' + str(bondstrength) + '\t' + str(res_e[2,2]) + '\n' )
            fo.write('\n')
            fo.write('timestep\t0.01\n')
            fo.write('thermo\t\t' + str(thermo_transform) + '\n')
            fo.write('run\t\t' + str(run_transform) + '\n')

        elif len(line.strip().split())>0 and line.strip().split()[0] == '#run_bonds_disassemble':
            fo.write('#run_bonds_disassemble\n')
            fil_length = (Nprotein-1)/2
            if fil_length%1==0:
                print('fil_length = ' + str(fil_length) + '\n')
                for cnt in range(int(fil_length)):
                    fo.write('bond_coeff\t' + str(9*cnt+1) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,0]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+2) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,1]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+3) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,2]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+4) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,0]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+5) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,1]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+6) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,2]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+7) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,0]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+8) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,1]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+9) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,2]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+0)) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,0]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+1)) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,1]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+2)) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,2]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+3)) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,0]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+4)) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,1]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+5)) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,2]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+6)) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,0]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+7)) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,1]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+8)) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,2]) + '\n' )
                    fo.write('pair_coeff\t1\t' + str(3*cnt+3) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('pair_coeff\t1\t' + str(3*cnt+4) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('pair_coeff\t1\t' + str(3*cnt+5) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('pair_coeff\t1\t' + str(3*(Nprotein-1)+3-3*cnt) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('pair_coeff\t1\t' + str(3*(Nprotein-1)+4-3*cnt) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('pair_coeff\t1\t' + str(3*(Nprotein-1)+5-3*cnt) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('\n')
                    fo.write('timestep\t0.01\n')
                    fo.write('thermo\t\t' + str(thermo_disassemble) + '\n')
                    fo.write('run\t\t' + str(run_disassemble) + '\n')
                    fo.write('\n')
            else:
                print('fil_length = ' + str(fil_length) + '\n')
                fil_length = fil_length -0.5
                print('new fil_length = ' + str(fil_length) + '\n')
                for cnt in range(int(fil_length)):
                    fo.write('bond_coeff\t' + str(9*cnt+1) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,0]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+2) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,1]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+3) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,2]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+4) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,0]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+5) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,1]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+6) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,2]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+7) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,0]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+8) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,1]) + '\n' )
                    fo.write('bond_coeff\t' + str(9*cnt+9) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,2]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+0)) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,0]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+1)) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,1]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+2)) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,2]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+3)) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,0]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+4)) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,1]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+5)) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,2]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+6)) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,0]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+7)) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,1]) + '\n' )
                    fo.write('bond_coeff\t' + str((Nprotein-1)*9 - (9*cnt+8)) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,2]) + '\n' )
                    fo.write('pair_coeff\t1\t' + str(3*cnt+3) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('pair_coeff\t1\t' + str(3*cnt+4) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('pair_coeff\t1\t' + str(3*cnt+5) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('pair_coeff\t1\t' + str(3*(Nprotein-1)+3-3*cnt) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('pair_coeff\t1\t' + str(3*(Nprotein-1)+4-3*cnt) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('pair_coeff\t1\t' + str(3*(Nprotein-1)+5-3*cnt) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                    fo.write('\n')
                    fo.write('timestep\t0.01\n')
                    fo.write('thermo\t\t' + str(thermo_disassemble) + '\n')
                    fo.write('run\t\t' + str(run_disassemble) + '\n')
                    fo.write('\n')
                cnt = cnt+1
                fo.write('bond_coeff\t' + str(9*cnt+1) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,0]) + ' \n' )
                fo.write('bond_coeff\t' + str(9*cnt+2) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,1]) + ' \n' )
                fo.write('bond_coeff\t' + str(9*cnt+3) + '\t' + str(bondstrength2) + '\t' + str(res_e[0,2]) + ' \n' )
                fo.write('bond_coeff\t' + str(9*cnt+4) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,0]) + ' \n' )
                fo.write('bond_coeff\t' + str(9*cnt+5) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,1]) + ' \n' )
                fo.write('bond_coeff\t' + str(9*cnt+6) + '\t' + str(bondstrength2) + '\t' + str(res_e[1,2]) + ' \n' )
                fo.write('bond_coeff\t' + str(9*cnt+7) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,0]) + ' \n' )
                fo.write('bond_coeff\t' + str(9*cnt+8) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,1]) + ' \n' )
                fo.write('bond_coeff\t' + str(9*cnt+9) + '\t' + str(bondstrength2) + '\t' + str(res_e[2,2]) + ' \n' )
                fo.write('pair_coeff\t1\t' + str(3*cnt+3) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                fo.write('pair_coeff\t1\t' + str(3*cnt+4) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                fo.write('pair_coeff\t1\t' + str(3*cnt+5) + '\t lj/cut \t ${eps_mu} \t ${sigma_mu} \t ${rc_mu} \n' )
                fo.write('\n')

        elif len(line.strip().split())>0 and line.strip().split()[0] == '#run_final':
            fo.write('#run_final\n')
            fo.write('timestep\t0.01\n')
            fo.write('thermo\t\t' + str(thermo_final) + '\n')
            fo.write('run\t\t' + str(run_final) + '\n')

        elif len(line.strip().split())>0 and line.strip().split()[0] == '#integrators':
            fo.write('fix\t\tfLANG all langevin 1.0 1.0 1 ' + str(seed) + ' zero yes omega yes\n')
            fo.write('fix\t\tfNVE\tmem\tnve/sphere update dipole\n')
            fo.write('fix\t\tfNVE2\tcyto\tnve/sphere update dipole\n')

        elif len(line.strip().split())>0 and line.strip().split()[0] == 'bond_coeff':
            print('skip')

        elif len(line.strip().split())>0 and line.strip().split()[0] == 'pair_coeff':
            print('skip')

        elif len(line.strip().split())>0 and line.strip().split()[0] == 'timestep':
            print('skip')

        elif len(line.strip().split())>0 and line.strip().split()[0] == 'thermo':
            print('skip')

        elif len(line.strip().split())>0 and line.strip().split()[0] == 'run':
            print('skip')

        else:
            fo.write(line)

    fo.close()

if __name__ == "__main__":
    main()
