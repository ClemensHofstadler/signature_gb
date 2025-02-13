from signature_gb import *

# amrhein
F = FreeAlgebra(QQ,['a','b','c','d','e','f'])
a,b,c,d,e,f = F.gens()
gens = [ 
    2*f*b+2*e*c+d**2+a**2+a,
    2*f*c+2*e*d+2*b*a+b, 
    2*f*d+e**2+2*c*a+c+b**2,
    2*f*e+2*d*a+d+2*c*b, 
    f**2+2*e*a+e+2*d*b+c**2,
    2*f*a+f+2*e*b+2*d*c,
    ]
armhein = LabelledModule(gens)