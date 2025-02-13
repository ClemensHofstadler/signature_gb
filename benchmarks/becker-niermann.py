from signature_gb import *

# becker-niermann 
F = FreeAlgebra(QQ,['x','y','z'])
x,y,z = F.gens()
gens = [ 
    x**2+x*y**2*z-2*x*y+y**4+y**2+z**2,
    -x**3*y**2+x*y**2*z+x*y*z**3-2*x*y+y**4,
    -2*x**2*y+x*y**4+y*z**4-3,
    ]
becker = LabelledModule(gens)