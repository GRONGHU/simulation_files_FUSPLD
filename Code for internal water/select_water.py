import numpy as np
import MDAnalysis as mda
import pyvista as pv
import pyKVFinder
from numpy.linalg import norm
import matplotlib
import matplotlib.pyplot as plt
import time
from MDAnalysis.analysis import distances
import pickle
import pymeshfix
import mdtraj as md

start =time.time()
# Input
PDB="mdtest_nowion.pdb"
GRO="mdtest_all.gro"
PLY="SURFACE_OUTPUT.ply"
u = mda.Universe(GRO)
reader =pv.get_reader(PLY)
mesh = reader.read()

#water molecules in clefts 
NEBOUR_WATER_ow=u.select_atoms("name OW and around 5 protein")
NEBOUR_WATER_ow_resid=NEBOUR_WATER_ow.resids
NEBOUR_WATER_mw=u.select_atoms("name MW and around 5 protein")
NEBOUR_WATER_mw_resid=NEBOUR_WATER_mw.resids
NEBOUR_WATER_hw1=u.select_atoms("name HW1 and around 5 protein")
NEBOUR_WATER_hw1_resid=NEBOUR_WATER_hw1.resids
NEBOUR_WATER_hw2=u.select_atoms("name HW2 and around 5 protein")
NEBOUR_WATER_hw2_resid=NEBOUR_WATER_hw2.resids
neibourwater_resid=list(set(NEBOUR_WATER_hw2_resid.tolist()+NEBOUR_WATER_hw1_resid.tolist()+NEBOUR_WATER_mw_resid.tolist()+NEBOUR_WATER_ow_resid.tolist()))
neibourwater_resid.sort()

points_ow = NEBOUR_WATER_ow.positions.tolist()
points_mw = NEBOUR_WATER_mw.positions.tolist()
points_hw1 = NEBOUR_WATER_hw1.positions.tolist()
points_hw2 = NEBOUR_WATER_hw2.positions.tolist()

points_poly_ow = pv.PolyData(points_ow) 
points_poly_mw = pv.PolyData(points_mw) 
points_poly_hw1 = pv.PolyData(points_hw1) 
points_poly_hw2 = pv.PolyData(points_hw2) 

try:
   select_ow=points_poly_ow.select_enclosed_points(mesh)
   select_mw=points_poly_mw.select_enclosed_points(mesh)
   select_hw1=points_poly_hw1.select_enclosed_points(mesh)
   select_hw2=points_poly_hw2.select_enclosed_points(mesh)
except:
   print("mesh_fixed")
   meshfix=pymeshfix.MeshFix(mesh)
   meshfix.repair()
   mesh_fix=meshfix.mesh
   select_ow=points_poly_ow.select_enclosed_points(mesh_fix)
   select_mw=points_poly_mw.select_enclosed_points(mesh_fix)
   select_hw1=points_poly_hw1.select_enclosed_points(mesh_fix)
   select_hw2=points_poly_hw2.select_enclosed_points(mesh_fix)


internalwater_ow=select_ow['SelectedPoints'].tolist()
internalwater_mw=select_mw['SelectedPoints'].tolist()
internalwater_hw1=select_hw1['SelectedPoints'].tolist()
internalwater_hw2=select_hw2['SelectedPoints'].tolist()

internalwater_resid_ow=[]
internalwater_resid_mw=[]
internalwater_resid_hw1=[]
internalwater_resid_hw2=[]

for i in range(len(internalwater_ow)):
    if internalwater_ow[i]==1:
        internalwater_resid_ow.append(NEBOUR_WATER_ow_resid[i])
for i in range(len(internalwater_mw)):
    if internalwater_mw[i]==1:
        internalwater_resid_mw.append(NEBOUR_WATER_mw_resid[i])
for i in range(len(internalwater_hw1)):
    if internalwater_hw1[i]==1:
        internalwater_resid_hw1.append(NEBOUR_WATER_hw1_resid[i])
for i in range(len(internalwater_hw2)):
    if internalwater_hw2[i]==1:
        internalwater_resid_hw2.append(NEBOUR_WATER_hw2_resid[i])

Select_water_resid=list(set(internalwater_resid_mw) & set(internalwater_resid_ow) & set(internalwater_resid_hw1) & set(internalwater_resid_hw2))

#remaining water molecules in water molecules
OTHERresid=[]
OTHERresid_mass_center=[]
for i in neibourwater_resid:
    if i not in Select_water_resid:
        OTHERresid.append(i)
        OTHERresid_mass_center.append(u.select_atoms("resid "+str(i)).center_of_mass())

results = pyKVFinder.run_workflow(PDB)
Carray=results.cavities
Carraylist=[]
countl2=0

for i in range(len(Carray)):
    LISTTEM=Carray[i].tolist()
    Carraylist.append(LISTTEM)

atomic = pyKVFinder.read_pdb(PDB)
probe_out = 4.0
step = 0.6
vertices = pyKVFinder.get_vertices(atomic, probe_out=probe_out, step=step)

waterinCav=[]
Xlimit=len(Carraylist)
Ylimit=len(Carraylist[0])
Zlimit=len(Carraylist[0][1])
print(Xlimit,Ylimit,Zlimit)
for item in OTHERresid_mass_center:
    countl2=0
    Xindex=int(abs(item[0]-vertices[0,0])/step)
    Yindex=int(abs(item[1]-vertices[0,1])/step)
    Zindex=int(abs(item[2]-vertices[0,2])/step)
    #print(Xindex,Yindex,Zindex)
    if Xindex+1>=Xlimit-1 or Yindex+1>=Ylimit-1 or Zindex+1>=Zlimit-1:
        waterinCav.append(0)
        continue
    for i in range(Xindex-1,Xindex+2):
        for j in range(Yindex-1,Yindex+2):
            for k in range(Zindex-1,Zindex+2):
                if Carraylist[i][j][k]>=2:
                    countl2=countl2+1

    if countl2>=1:
        waterinCav.append(1)
    else:
        waterinCav.append(0)

Select_water_resid_CAV=[]
for i in range(len(OTHERresid)):
    if waterinCav[i]==1:
        Select_water_resid_CAV.append(OTHERresid[i])

end=time.time()
print('Running time: %s Seconds'%(end-start))    

       
#save
output_cavity = 'cavity.pdb'
pyKVFinder.export(output_cavity, results.cavities, results.surface, vertices, step=step)
file_save=open("Select_water.npy","ab")
pickle.dump([Select_water_resid,Select_water_resid_CAV],file_save)
file_save.close()