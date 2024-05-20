import numpy as np
import MDAnalysis as mda
import time
import pickle
import sys
import os
import logging
from pathos.multiprocessing import ProcessingPool as Pool
from pathlib import Path

#Multiprocess version
def SASAfilter_multip(Internal_Water,u,start,end,SASA_limit):
    Atomnumber=len(u.select_atoms("protein"))
    Internal_water_filter=[]
    IW_resid=""
    
    if len(Internal_Water)==0:
        return Internal_water_filter
    else:
        for i in Internal_Water:
            IW_resid=IW_resid+" "+str(i)
    protein_IW=u.select_atoms("protein or (resid"+IW_resid +" and not type DUMMY )")

    temp_gro ="temp_"+str(start)+"_"+str(end)+".gro"
    try:
        protein_IW.write(temp_gro)
        traje=md.load(temp_gro,top=temp_gro)
        sasa = md.shrake_rupley(traje)
        IW_sasa=[] 
        for i in range(len(Internal_Water)):
            IW_sasa.append(sum(sasa[0][Atomnumber+i*3:Atomnumber+3+i*3]))
        for j in range(len(IW_sasa)):
                if IW_sasa[j]<=SASA_limit:
                    Internal_water_filter.append(Internal_Water[j])
    finally:
        if os.path.exists(temp_gro):
            os.remove(temp_gro)

    return Internal_water_filter

def getPDBwithWater_multip(Internal_water,u,start,end,SASA_limit,output=False):
    Internal_water_filter=SASAfilter_multip(Internal_water,u,start,end,SASA_limit)
    temp=Internal_water_filter.copy()
    item=1
    while True:
        Internal_water_filter=SASAfilter_multip(Internal_water_filter,u,start,end,SASA_limit)
        if temp==Internal_water_filter:
            break
        else:
            temp=Internal_water_filter.copy()
            item=item+1
    if output==True:
        IW_resid_filter=""
        for j in Internal_water_filter:
            IW_resid_filter=IW_resid_filter+" "+str(j)
            protein_IW_filter=u.select_atoms("protein or (resid"+IW_resid_filter+" and not type DUMMY )")
            protein_IW_filter.write("WaterinAsa_filter"+str(start)+"_"+str(end)+".gro")
    return Internal_water_filter


def process_range(args):
    import MDAnalysis as mda
    import pickle
    import sys
    import logging
    import os
    start, end, GRO,XTC,SURF,SASA_limit = args
    u = mda.Universe(GRO,XTC)
    logging.basicConfig(filename=f"process_range_{start}_{end}.log", filemode='w', level=logging.INFO)
    SURF_water=[]
    file=open(SURF,"rb")
    for i in range(len(u.trajectory)):
        waterlist=pickle.load(file)
        water=waterlist[0]+waterlist[1]
        SURF_water.append(list(set(water)))
    file.close()
    results = []
    for ts in u.trajectory[start:end]:
        frame=u.trajectory.frame
        neibourwater_resid = SURF_water[frame].copy()
        neibourwater_resid.sort()
        Internal_water_filter =getPDBwithWater_multip(neibourwater_resid, u,start, end, SASA_limit, output=False)
        logging.info(f"Frame: {frame}, IW: {len(Internal_water_filter)}")
        results.append(Internal_water_filter)
    return results

GRO="md.gro"
XTC="md1.xtc"
SURF="Select_water.npy"
SASA_limit=1.2*0.20

u = mda.Universe(GRO,XTC)
n_frames = len(u.trajectory)
n_processes = 4
frames_per_process = n_frames // n_processes
frame_ranges = [(i * frames_per_process, (i + 1) * frames_per_process) for i in range(n_processes)]
frame_ranges[-1] = (frame_ranges[-1][0], n_frames)
with Pool(n_processes) as p:
    args = [(start, end, GRO, XTC,SURF, SASA_limit) for start, end in frame_ranges]
    results = p.map(process_range, args)

file=open("InternalWater.npy","wb")
pickle.dump(results,file)
file.close()