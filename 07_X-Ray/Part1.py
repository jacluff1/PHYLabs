import Operations1 as op
import Operations2 as op2


#########################################################################

# PARAMETERS

###########################################################################

CallAll = False # prints and plots all lab results
plotA = True
printA = True
saveA = True

##################################################################################

if CallAll == True:
    #op.PartIII(case,plotA,saveA,printA)
    op.PartIII(1,True,True,True) # Fig 3.1 Find A,B
    op.PartIII(2,True,True,True) # Fig 3.2 Full Spectrum
    op.PartIII(3,True,True,True) # Fig 3.3 Mosley's Law
    op.PartIV(1,True,True,True) # Fig 4.1 Find Position and Intensity
    op.PartIV(2,True,True,True) # Fig 4.2 Power Law Fit
    op.PartV(1,True,True,True) # Fig 5.1 
    op.PartV(2,True,True,True) # Fig 5.1 Attenuation length
    op.PartV(3,True,True,True) # Fig 5.2 mu(E)
    op.PartV(4,True,True,True) # Fig 5.3 Absorption Edge
    op2.Part(1,True,True,True) # Fig 6.1 and 6.2

else:
    #op.PartIII(case,plotA,saveA,printA)
    #op.PartIII(1,plotA,saveA,printA) # Fig 3.1 Find A,B
    #op.PartIII(2,plotA,saveA,printA) # Fig 3.2 Full Spectrum
    #op.PartIII(3,plotA,saveA,printA) # Fig 3.3 Mosley's Law
    #op.PartIV(1,plotA,saveA,printA) # Fig 4.1 Find Position and Intensity
    #op.PartIV(2,plotA,saveA,printA) # Fig 4.2 Power Law Fit
    op.PartV(1,plotA,saveA,printA) # Fig 5.1 Spectrum
    #op.PartV(2,plotA,saveA,printA) # Fig 5.2 Attenuation length
    #op.PartV(3,plotA,saveA,printA) # Fig 5.3 Attenuation
    #op.PartV(4,plotA,True,printA) # Fig 5.3 Absorption Edge
    #op2.Part(1,plotA,saveA,printA) # Fig 6.1 and 6.2
    
