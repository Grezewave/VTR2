import OSfunct
import os

def colorscale(VMD, cutoff, out):
    #return a rgb code color between blue and red, based on cutoff and VMD, bigger the VMD, bluest the color
    Redest = [255,0,0]
    Bluest = [0,0,255]
    R = (((-255)/cutoff)*VMD)+255
    G = 0
    B = (((255)/cutoff)*VMD)
    if out == 'l':
        color = [int(R),int(G),int(B)]
    elif out == 't':
        color = (float(B)/255,float(G),float(R)/255)
    return(color)
