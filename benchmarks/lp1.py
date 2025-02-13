from signature_gb import *

# lp1 
F = FreeAlgebra(QQ,['x','y','z'])
x,y,z = F.gens()
gens = [
    z**4 + y*x*y*x - x*y**2*x - 3*z*y*x*z, 
    x**3 + y*x*y - x*y*x, 
    z*y*x - x*y*z + z*x*z
    ]
lp1 = LabelledModule(gens)
