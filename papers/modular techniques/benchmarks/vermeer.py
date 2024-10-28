from signature_gb import *

# vermeer 
F = FreeAlgebra(QQ,['w','p','z','t','s'])
w,v,u,y,x = F.gens()
gens = [
    v**2+u**2-2*v*y+y**2-2*u*x+x**2-1,
    -u**3+v**2,
    -3*v*u**2+3*u**2*y-2*v*u+2*v*x,
    6*w**2*v*u**2-3*w*u**2-2*w*v+1,
    ]
vermeer = LabelledModule(gens)