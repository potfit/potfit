import numpy as np
import shutil
import os

natom = 512 

#
#num_param = 11
#num_config = 500
#num_block = 10
#num_config_per_block = num_config/num_block
#


# construct training set 
with open('lammps.dump','r') as fin:
  raw_data = fin.readlines()
  numc = 0
  ll = 0

  while ll < range(len(raw_data)):

    if(ll == len(raw_data)):
      break
    
    row = raw_data[ll]
    if row[0:9] == 'ITEM: BOX':
      # parsing parameters 
      box_size = []
      for ii in range(3):
        ll += 1
        tmp_line = raw_data[ll]
        tmp_line = tmp_line.split()
        box_size.append(float(tmp_line[1]))

      # parsing positions and forces 
      fname = 'training_set'+str(numc)
      with open (fname, 'w') as fout:
        fout.write('#N %d 1\n'%natom)
        fout.write('#C Si\n')
        fout.write('#X %18.10e 0.0 0.0\n'%box_size[0])
        fout.write('#Y  0.0 %18.10e 0.0\n'%box_size[1])
        fout.write('#Z  0.0 0.0 %18.10e\n'%box_size[2])
        fout.write('#E 0.0\n')
        fout.write('#F\n')
       
        ll += 1
        for ii in range(natom): 
          ll += 1
          tmp_line = raw_data[ll]
          tmp_line = tmp_line.split()
          fout.write('0 ')
          for jj in range(2,8):
            fout.write(tmp_line[jj]+'  ')
          fout.write('\n')  

      # update numc info
      numc += 1
    # update line info
    ll += 1


