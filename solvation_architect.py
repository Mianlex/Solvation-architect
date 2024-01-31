import numpy as np
from ase.build import fcc111, add_adsorbate, sort
from ase.io import read, write
import random
import copy

def read_slab(file_path):
    """
    Read the slab from a file and return the ASE Atoms object.
    """
    return read(file_path)

def add_water_molecules(slab, water, suflayerindex, watern, space):
    """
    Add water molecules to the slab randomly.
    The first parameter is the atomic objects of pure slab surface.
    The second parameter is the atomic objects of single water molecule.
    The thrid parameter is the index list of first slab surface layer.
    The 4th parameer is number of water molecules you want to add at surface.
    The 5th parameer is space parameter that the z length region to fill water molecules.
    """
    i = 0
    while i < watern:
        #generate the structure as random as possible
        index = random.choice(suflayerindex)
        dev = random.choice([-0.1*np.random.rand(), 0.1*np.random.rand()])
        water.rotate((np.random.rand(), np.random.rand(), np.random.rand()),
                     (np.random.rand(), np.random.rand(), np.random.rand()))
        #this add_adsorbate parameters has tuned and worked well for 4x4x4 metal surface with 40 H2O
        #Tuning it to match your requirment
        add_adsorbate(slab, water, np.random.rand()*space + 2.3,
                      (slab.positions[index][0]+dev, slab.positions[index][1]+dev),
                      offset=(0.13*np.random.rand(), 0.13*np.random.rand()))
        if i > 0:
            #Check the distance between the added water molecule
            flag_go=False
            num=slab.get_atomic_numbers()
            O=[x for x, y in list(enumerate(num)) if y ==8] #index for Oxygen atoms
            H=[x for x, y in list(enumerate(num)) if y ==1] #index for H atoms
            OH=[x for x, y in list(enumerate(num)) if y==8 or y==1] #index for O and H
            dOO=1.21*1.1 #1.1 noise of bond length
            dOH=0.96 #cannot control noise because  H---O-H  hydrogen bond
            dHH=0.8*1.2 # 1.2 noise of bond length
            for u in enumerate(O):
                disa=slab.get_distances(u[1],O)#OOdistance
                disa=np.sort(disa)
                mindis=disa[1]
                if mindis<dOO:
                    flag_go=True
                    break
                disb=slab.get_distances(u[1],H)#OHdistance
                disb=np.sort(disb)
                mindisb=disb[0]
                if mindisb<dOH:
                    flag_go=True
                    break
            for h in enumerate(H):
                dish=slab.get_distances(h[1],H)#HHdistance
                dish=np.sort(dish)
                mindish=dish[1]
                if mindish<dHH:
                    flag_go=True
                    break

            if flag_go:
                del slab[-3:]
            else:
                i = i + 1
        else:
            i = i + 1
#    return slab

def add_ads(slab,mol,index,height):
    """
    Compulsory add adsorbate.
    The first parameter is the atomic objects of slab surface(+random water).
    The second parameter is the atomic objects of adsorbate.
    The thrid parameter is the index of specific atom you want to put adsorbate.
    The 4th parameer is height parameter of adsorbate.
    Tuning height parameter can put the adsorbate at water phase, not always at metal surface.
    """
    a=np.random.rand()
    add_adsorbate(slab,mol,height+a,(slab.positions[index][0],slab.positions[index][1]))

def remove_closest_water(slab, MCindex, num_water_to_remove):
    """
    Remove the closest water molecules to a specific carbon atom.
    The first parameter is the atomic objects of surface+water+adsorbate.
    The second parameter is the index of the centeral atoms in adsorbate.In glycerol, it is the middle carbon.
    The thrid parameter is the number of water molecules that the closest to ads to remove.
    """
    num = slab.get_atomic_numbers()
    O = [x for x, y in enumerate(num) if y == 8][:-3]  # index for Oxygen atoms in H2O # 3 is the num O in Glycerol
    disco = slab.get_distances(MCindex, O)  # C--O distance
    sorted_indices = np.argsort(disco)
    related_numbers = [O[i] for i in sorted_indices][0:num_water_to_remove]
    delindex = []
    for i in related_numbers:
        delindex.append(i)#O atoms
        delindex.append(i+1)#H atoms
        delindex.append(i+2)#H atoms

    del slab[delindex]


    
def separate_and_merge(slab,slab2,num_ads_atoms):
    """
    Separate and merge the slab, water trajectories, and adsorbate trajectories.
    The first parameter is the atomic objects of surface+water+adsorbate.
    The sceond para is the pure slab surface.
    The thrid parameter is the number of atoms of adsorbate.
    """
    #slab = read(slab_traj)
    con = slab
    con2 = copy.deepcopy(slab)
    del con[-1*num_ads_atoms:]#delete adsorbate
    del con[0:64]#delete metal
    
    
    #sort H & O in O , H order
    custom_order = {
        "O": 0,
        "H": 1
    }
    tagslist = sorted(con.get_chemical_symbols(), key=lambda x: custom_order[x])
    con = sort(con)
    con = sort(con, tags=tagslist)#first iteration sort
    con = sort(con, tags=tagslist)#second iteration sort
    con.write('watertraj.traj') #write out the water traj

    del con2[:-1*num_ads_atoms]#catch adsorbate

    slab2.extend(con)
    slab2.extend(con2)
    slab2.write('POSCAR', format="vasp")

# Example usage
slab = read('slab.traj')
slab2 = read('slab.traj') #copy.deepcopy(slab)
water = read('H2O.traj')
mol = read('gly.traj')
num_ads_atoms=len(mol)
print(num_ads_atoms)
heights = np.unique([sl.z for sl in slab])
#fix_heights = np.sort(heights)[0:2]
#c = FixAtoms([atom.index for atom in atoms if atom.z in fix_heights])
#atoms.set_constraint(c)
suflayerindex=[atom.index for atom in slab if atom.z in heights[-16:]]#first layer
watern = 45
space = 11
add_water_molecules(slab, water, suflayerindex, watern,space)

index=58 # central Pt index
height=2.3 # adsorbate height
add_ads(slab, mol, index,height)
MCindex = 200 # middle carbon index of glycerol
num_water_to_remove = 5
remove_closest_water(slab, MCindex, num_water_to_remove)
separate_and_merge(slab,slab2,num_ads_atoms)
