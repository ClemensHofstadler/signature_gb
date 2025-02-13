from signature_gb import *

# metab3
F = FreeAlgebra(QQ,['x1','x2','x3'])
x1,x2,x3 = F.gens()
gens = [
    x1*x3*x2*x3 - x1*x3**2*x2 - x2*x3*x1*x3 + x2*x3**2*x1 - x3*x1*x2*x3 + x3*x1*x3*x2 + x3*x2*x1*x3 - x3*x2*x3*x1,
    x1*x2*x1*x3 - x1*x2*x3*x1 - x1*x3*x1*x2 + x1*x3*x2*x1 - x2*x1**2*x3 + x2*x1*x3*x1 + x3*x1**2*x2 - x3*x1*x2*x1,
    x1*x2**2*x3 - x1*x2*x3*x2 - x2*x1*x2*x3 + x2*x1*x3*x2 - x2*x3*x1*x2 + x2*x3*x2*x1 + x3*x2*x1*x2 - x3*x2**2*x1,
    -x1*x3*x2*x3 + x1*x3**2*x2 + x2*x3*x1*x3 - x2*x3**2*x1 + x3*x1*x2*x3 - x3*x1*x3*x2 - x3*x2*x1*x3 + x3*x2*x3*x1,
    -x1*x2*x1*x3 + x1*x2*x3*x1 + x1*x3*x1*x2 - x1*x3*x2*x1 + x2*x1**2*x3 - x2*x1*x3*x1 - x3*x1**2*x2 + x3*x1*x2*x1,
    -x1*x2**2*x3 + x1*x2*x3*x2 + x2*x1*x2*x3 - x2*x1*x3*x2 + x2*x3*x1*x2 - x2*x3*x2*x1 - x3*x2*x1*x2 + x3*x2**2*x1
    ]
metab3 = LabelledModule(gens)