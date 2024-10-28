from signature_gb import *

# cyclic5
F = FreeAlgebra(QQ,['v','w','x','y','z'])
v,w,x,y,z = F.gens()
gens = [ 
    v+w+x+y+z,
    v*w+w*x+x*y+v*z+y*z,
    v*w*x+w*x*y+v*w*z+v*y*z+x*y*z,
    v*w*x*y+v*w*x*z+v*w*y*z+v*x*y*z+w*x*y*z,
    v*w*x*y*z-1,
    ]
cyclic5 = LabelledModule(gens)