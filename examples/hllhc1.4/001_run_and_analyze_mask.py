# %% Import packages
from cpymad.madx import Madx


# %% Check the executable
import sys
print(sys.executable)

# %% Changing folder
import os
# please change it conveniently

# os.chdir('/afs/cern.ch/work/s/sterbini/beambeam_macros/examples')
print(os.getcwd())
# %% Unmask the mask
beam = 1
i_octupoles = 100.
emittance_in_um = 2.3
n_particles = 2.25e11
chromaticity = 15
xing_angle_urad = 245.
seedran = 1

fname_mask = 'hl14_template.mask'

with open(fname_mask) as fid:
    mask_content = fid.read()

mask_content = mask_content.replace(r'%BEAM%', str(1))
mask_content = mask_content.replace(r'%OCT%', f'{i_octupoles:e}')
mask_content = mask_content.replace(r'%EMIT_BEAM', f'{emittance_in_um:e}')
mask_content = mask_content.replace(r'%NPART', f'{n_particles:e}')
mask_content = mask_content.replace(r'%CHROM%', f'{chromaticity:e}')
mask_content = mask_content.replace(r'%XING', f'{xing_angle_urad:e}')
mask_content = mask_content.replace(r'%SEEDRAN', f'{seedran:d}')
# %% Dump the unmasked mask on file
with open(fname_mask.split('.mask')[0]+'_unmask.mask', 'w') as fid:
    fid.write(mask_content)

# %% split the mask
# I am assuming that the start of the file is '! %%'
assert(mask_content[0:4]=="! %%")
aux=mask_content.split("! %%")
title=[]
body=[]
for i in aux:
    title.append(i.split('\n')[0])
    body.append("".join([a+'\n' for a in i.split('\n')[1:]]))
# I remove the first (empty line, see assertion above)
title=title[1:]
body=body[1:]
# %% built a pandas df
import pandas as pd
myDF=pd.DataFrame(body,index=title, columns=['Code string'])

# %% Run MADX
madx = Madx()

# %%
import time
myGlobals=[]
mySmallDF=myDF
myString=''
for block in mySmallDF.iterrows():
    print(block[0])
    start_time = time.time()
    myString=myString+ '! %%' +block[0]+ '\n' + block[1]['Code string'][0:-1]
    with madx.batch():
        madx.input('! %%' +block[0]+ '\n' + block[1]['Code string'][0:-1])
    execution_time_s=time.time()-start_time
    myDict={}
    myDict=dict(madx.globals)
    myDict['execution time [s]']=execution_time_s
    myGlobals.append(myDict)
profileDF=pd.DataFrame(myGlobals, index=mySmallDF.index)
# %%
with open('inverse.mask', 'w') as fid:
    fid.write(myString)
# %%
# %%
madx.input(aux[2])


# %%
