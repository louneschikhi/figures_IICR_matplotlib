
from pylab import plot, figure, savefig, text, axis, xlabel, ylabel, axes, arrow, grid, title
import math
import scipy
import matplotlib.pyplot as plt


################################################################################
##### Function representing the change in population size that explains 
##### the distribution of coalescence times for a sample of size 2
##### when the two genes are sampled in the SAME deme
##### t: the time expressed in units of N
##### t=1 thus corresponds to N generatiions in the past
################################################################################


def lambdaS(t) :
    num = (1- beta)*math.exp(-alpha*t) + (alpha -1)* math.exp(-beta*t)
    denom = (alpha- gamma)*math.exp(-alpha*t) + (gamma -beta)* math.exp(-beta*t)
    res = num/denom
    return res

################################################################################

################################################################################
##### Function representing the change in population size that explains 
##### the distribution of coalescence times for a sample of size 2
##### when the two genes are sampled in the SAME deme
##### t: the time expressed in units of N
##### t=1 thus corresponds to N generations in the past
##### n: the number of islands
##### M: the number of immigrant genes per generation per island
##### The model assumes that each island had N haploid genomes
##### If m is the migration rate then M= Nm
################################################################################

def lambdaMnS (t,M,n) :
    gamma = M/float(n-1)
    delta = pow((1.0 + n*gamma), 2.0) - 4*gamma
    alpha = (1.0 + n*gamma + math.sqrt(delta))/2.0
    beta =  (1.0 + n*gamma - math.sqrt(delta))/2.0
    a = (gamma - alpha)/(beta - alpha)
    c = gamma/(beta - alpha)
    num = (1.0- beta)*math.exp(-alpha*t) + (alpha -1.0)* math.exp(-beta*t)
    denom = (alpha- gamma)*math.exp(-alpha*t) + (gamma -beta)* math.exp(-beta*t)
    res = num/denom
    return res


################################################################################


################################################################################
##### Function representing the change in population size that explains 
##### the distribution of coalescence times for a sample of size 2
#####
##### when the two genes are sampled in the DIFFERENT demes
#####
##### t: the time expressed in units of N
##### t=1 thus corresponds to N generations in the past
##### n: the number of islands
##### M: the number of immigrant genes per generation per island
##### The model assumes that each island had N haploid genomes
##### If m is the migration rate then M= Nm
################################################################################

def lambdaMnD (t,M,n) : 
    gamma = M/float(n-1)
    delta = pow((1.0 + n*gamma), 2.0) - 4*gamma
    alpha = (1.0 + n*gamma + math.sqrt(delta))/2.0
    beta =  (1.0 + n*gamma - math.sqrt(delta))/2.0
    a = (gamma - alpha)/(beta - alpha)
    c = gamma/(beta - alpha)
    num = beta*math.exp(-alpha*t) - alpha*math.exp(-beta*t)
    denom = gamma*math.exp(-alpha*t) - gamma*math.exp(-beta*t)
    res = num/denom
    return res

################################################################################

####################################
# Plot with variable M  
# and n= 50
# of lambdaMnS
####################################

x = [0.01*i for i in range(601)] # vector for the x-axis
n=50

figure()

color = "black"

plt.axhline(y=50, xmin=0, xmax=5, linestyle="--", color =color)

for M in [50, 20, 10, 5, 2, 1, 0.5] :
    y = []
    for i in range(len(x)) :
        y.append(lambdaMnS(x[i], M, n))
    plot(x,y, color=color)


#text(scipy.mean(x), scipy.mean(y), "$\lambda(t)$" ) 
text(4.7, 130, "M=0.5")
text(3.5, 95, "M=1.0")
text(2.5, 75, "M=2")
text(1.4, 64, "M=5")
text(0.6, 57, "M=10")
text(0.08, 53, "M=20")
text(-0.50, 42, "M=50")

title("(a) IICR for genes sampled in the same deme")
axis([min(x)-max(x)/10, max(x)*1.1, min(y)-max(y)/10, max(y)*1.1])
xlabel("scaled time ($\t{t}$)")
ylabel("IICR $\lambda(t)$")

savefig("Fig_IICRs_n50varM_50_05.png")


####################################
# Plot with variable M  
# and n= 50
# of lambdaMnD
####################################

####
# note that here we create x by multiplying by (i+1) 
# to avoid computing the IICR at x=0
####

x = [0.001*(i+1) for i in range(10001)] # vector for the x-axis
n=50

figure()
color = "black"

plt.axhline(y=50, xmin=0, xmax=5, linestyle="--", color=color)

for M in [50, 20, 10, 5, 2, 1, 0.5] :
    y = []
    for i in range(len(x)) :
        y.append(lambdaMnD(x[i], M, n))
    plot(x,y, color=color)


text(0.5, 350, "M=0.5")
text(0.35, 250, "M=1.0")
text(0.22, 180, "M=2")
text(0.10, 145, "M=5")
text(0.095, 300, "M=10")
text(0.05, 500, "M=20")
text(-0.13, -70, "M=50")

#ax.arrow(0, 0, 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')

axes()
#arrow(0.08, 300, 0.05, lambdaMnD(0.05, 10, 50), width=0.001) #head_width=0.05, head_length=0.1)
arrow(0.1, 290, -0.05, (lambdaMnD(0.05, 10, 50) - 290)+30, width=0.001, head_width=0.01, head_length=24, color=color, linestyle="dotted")

arrow(0.05, 470, -0.02, lambdaMnD(0.03, 20, 50) - 470+30, width=0.001, head_width=0.01, head_length=25, color=color, linestyle="dotted")

arrow(-0.03, -10, 0.04, lambdaMnD(0.01, 50, 50)+10-60, width=0.001, head_width=0.01, head_length=15, color=color, linestyle="dotted")


#grid(True)
title("IICR for genes sampled in different demes")
axis([-0.15, 1, -100, 1500])
xlabel("scaled time ($\t{t}$)")
ylabel("IICR $\lambda(t)$")

savefig("Fig_IICRd_n50varM_50_05.png")

################################################
##########   figure with subplots  #############
# Plot with variable M  
# and n= 50
# of lambdaMnS
# and lambdaMnD

################################################




x = [0.01*i for i in range(601)] # vector for the x-axis
n=50

fig = figure()

### First subplot

figure.add_subplot(211)

color = "black"

plt.axhline(y=50, xmin=0, xmax=5, linestyle="--", color =color)

for M in [50, 20, 10, 5, 2, 1, 0.5] :
    y = []
    for i in range(len(x)) :
        y.append(lambdaMnS(x[i], M, n))
    plot(x,y, color=color)


#text(scipy.mean(x), scipy.mean(y), "$\lambda(t)$" ) 
text(4.7, 130, "M=0.5")
text(3.5, 95, "M=1.0")
text(2.5, 75, "M=2")
text(1.4, 64, "M=5")
text(0.6, 57, "M=10")
text(0.08, 53, "M=20")
text(-0.50, 42, "M=50")

title("(a) IICR for genes sampled in the same deme")
axis([min(x)-max(x)/10, max(x)*1.1, min(y)-max(y)/10, max(y)*1.1])
xlabel("scaled time ($\t{t}$)")
ylabel("IICR $\lambda(t)$")



### Second subplot

figure.add_subplot(212)


x = [0.001*(i+1) for i in range(10001)] # vector for the x-axis
n=50

figure()
color = "black"

plt.axhline(y=50, xmin=0, xmax=5, linestyle="--", color=color)

for M in [50, 20, 10, 5, 2, 1, 0.5] :
    y = []
    for i in range(len(x)) :
        y.append(lambdaMnD(x[i], M, n))
    plot(x,y, color=color)


text(0.5, 350, "M=0.5")
text(0.35, 250, "M=1.0")
text(0.22, 180, "M=2")
text(0.10, 145, "M=5")
text(0.095, 300, "M=10")
text(0.05, 500, "M=20")
text(-0.13, -70, "M=50")

#ax.arrow(0, 0, 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')

axes()
#arrow(0.08, 300, 0.05, lambdaMnD(0.05, 10, 50), width=0.001) #head_width=0.05, head_length=0.1)
arrow(0.1, 290, -0.05, (lambdaMnD(0.05, 10, 50) - 290)+30, width=0.001, head_width=0.01, head_length=24, color=color, linestyle="dotted")

arrow(0.05, 470, -0.02, lambdaMnD(0.03, 20, 50) - 470+30, width=0.001, head_width=0.01, head_length=25, color=color, linestyle="dotted")

arrow(-0.03, -10, 0.04, lambdaMnD(0.01, 50, 50)+10-60, width=0.001, head_width=0.01, head_length=15, color=color, linestyle="dotted")


#grid(True)
title("(b) IICR for genes sampled in different demes")
axis([-0.15, 1, -100, 1500])
xlabel("scaled time ($\t{t}$)")
ylabel("IICR $\lambda(t)$")


savefig("Fig_IICRsd_n50varM_50_05.png")




