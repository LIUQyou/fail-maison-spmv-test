import ssgetpy
import os
r_max=23947347
r_min=534388
result=ssgetpy.search(rowbounds=(r_min,r_max), limit=1000)
for r in result:
    if r.id <= 1856:
        continue
    if r.dtype=='complex':
        continue
    if r.dtype=='matrix':
        continue
    if r.dtype=='coordinate':
        continue
    if r.dtype=='hermitian':
        continue
    cmd='wget https://suitesparse-collection-website.herokuapp.com/MM/'+r.group+'/'+r.name+'.tar.gz'
    os.system(cmd)
    cmd='mv '+r.name+'.tar.gz ../spm'
    os.system(cmd)
    cmd='tar -zxvf ../spm/'+r.name+'.tar.gz -C ../spm'
    os.system(cmd)
    cmd='bash profile.sh csr ../spm/' + r.name + '/' + r.name + '.mtx'
    os.system(cmd)