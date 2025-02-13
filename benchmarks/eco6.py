from signature_gb import *

# eco6 
F = FreeAlgebra(QQ,['x1','x2','x3','x4','x5','x6'])
x1,x2,x3,x4,x5,x6 = F.gens()
gens = [ 
    x1+x2+x3+x4+x5+1
    ,x5*x6-5
    ,x1*x5*x6+x4*x6-4
    ,x1*x4*x6+x2*x5*x6+x3*x6-3
    ,x1*x3*x6+x2*x4*x6+x3*x5*x6+x2*x6-2
    ,x1*x2*x6+x2*x3*x6+x3*x4*x6+x4*x5*x6+x1*x6-1
    ]
eco6 = LabelledModule(gens)