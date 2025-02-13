from signature_gb import *

# P6 
F = FreeAlgebra(QQ,['a','b','c'])
a,b,c = F.gens()
gens = [
    c*c*c + 2*c*c*b + 3*c*c*a + 5*b*c*c + 7*a*c*a, 
    b*c*c + 11*b*a*b + 13*a*a*a
    ]
P6 = LabelledModule(gens)