from signature_gb import *

# cyclic6
F = FreeAlgebra(QQ,['x1','x2','x3','x4','x5','x6'])
x1,x2,x3,x4,x5,x6 = F.gens()
gens = [ 
    x1+x2+x3+x4+x5+x6,
    x1*x2+x2*x3+x3*x4+x4*x5+x1*x6+x5*x6,
    x1*x2*x3+x2*x3*x4+x3*x4*x5+x1*x2*x6+x1*x5*x6+x4*x5*x6,
    x1*x2*x3*x4+x2*x3*x4*x5+x1*x2*x3*x6+x1*x2*x5*x6+x1*x4*x5*x6+x3*x4*x5*x6,
    x1*x2*x3*x4*x5+x1*x2*x3*x4*x6+x1*x2*x3*x5*x6+x1*x2*x4*x5*x6+x1*x3*x4*x5*x6+x2*x3*x4*x5*x6,
    x1*x2*x3*x4*x5*x6-1,
    ]
cyclic6 = LabelledModule(gens)