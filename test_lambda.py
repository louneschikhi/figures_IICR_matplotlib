
from pylab import plot, figure, savefig, text, axis, xlabel, ylabel, axes, arrow, grid, title, show
import math
import scipy
import matplotlib.pyplot as plt


################################################################################
##### Function representing the change in population size that explains 
##### the distribution of coalescence times for a sample of size 2
##### when the two genes are sampled in the SAME deme
##### t: the time expressed in units of N
##### t=1 thus corresponds to N generatiions in the past
##### NOTE: THIS FUNCTION WAS A TEST AND ASSUMED THAT 
##### THE PARAMETERS alpha, betwa, etc, had been computed before
##### YOU CAN SKIP IT.
################################################################################


def lambdaS(t) :
    num = (1- beta)*math.exp(-alpha*t) + (alpha -1)* math.exp(-beta*t)
    denom = (alpha- gamma)*math.exp(-alpha*t) + (gamma -beta)* math.exp(-beta*t)
    res = num/denom
    return res

################################################################################

################################################################################
##### Function representing the change in population size that explains 
##### the distribution of coalescence times for a sample of size 2 (genes)
##### when the two genes are sampled in the SAME deme
#####
##### t: the time expressed in units of N
##### t=1 thus corresponds to N generations in the past
##### n: the number of islands/demes
##### M: the number of immigrant genes per generation per island
##### The model assumes that each island had N haploid genomes
##### If m is the migration rate then M= Nm
##### See Mazet, Rodriguez and Chikhi (2015) TPB for details.
#####
##### This function correspoinds to equation (X) in Mazet et al. (submitted)
################################################################################

def lambdaMnS(t,M,n) :
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

##### See Mazet, Rodriguez and Chikhi (2015) TPB for details.
#####
##### This function corresponds to equation (X) in Mazet et al. (submitted)
################################################################################

def lambdaMnD(t,M,n) : 
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
n=50  # nb of demes

figure()  # prepare the plot

color = "black" # just because we do not want color for the publication

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

title("$(a) IICR: \lambda_s(t)$ (n=50, variable migration rates)")
axis([min(x)-max(x)/10, max(x)*1.1, min(y)-max(y)/10, max(y)*1.1])
xlabel("scaled time ($\t{t=T/N}$)")
ylabel("$IICR: \lambda_s(t)$")

savefig("Fig_IICRs_n50varM_50_05.png")

####################################
# Plot variable n with 
# the STRONG  MIGRATION limit
# M=100 and M=500
####################################

x = [0.001 * i for i in range(1001)]

# prepare figure plot:
figure()
color = "black"
fontsize=12

for M in [100, 500] :
    for n in [50, 100, 200, 400, 500] :
        y = []
        for i in range(len(x)) :
            y.append(lambdaMnS(x[i], M, n))
        if M==100 :
	  plot(x,y, color=color, linestyle=":")
	else :
	  plot(x,y, color=color)
        plt.axhline(y=n, xmin=0, xmax=5, linestyle="--", color =color)

## Add text to figure
text(0.1, 450, "n=500, M=100", fontsize=fontsize)
text(0.1, 350, "n=400, M=100", fontsize=fontsize)
text(0.1, 150, "n=200, M=100", fontsize=fontsize)
text(0.1, 80, "n=100, M=100", fontsize=fontsize)
text(0.1, 30, "n=50, M=100", fontsize=fontsize)

text(-0.17, 450, "n=500, M=500", fontsize=fontsize)
text(-0.17, 350, "n=400, M=500", fontsize=fontsize)
text(-0.17, 150, "n=200, M=500", fontsize=fontsize)
text(-0.17, 80, "n=100, M=500", fontsize=fontsize)
text(-0.17, 30, "n=50, M=500", fontsize=fontsize)

## Prepare title, axes, labels etc.
title("$(b) IICR: \lambda_s(t)$ (strong migration)")
axis([min(x)-max(x)/5, max(x)*0.6, min(y)-max(y)/10, max(y)*1.15])
xlabel("scaled time ($\t{t=T/N}$)")
ylabel("$IICR: \lambda_s(t)$")

savefig("Fig_IICRs_nvarM_100_500.png")

####################################
# Plot with variable M  
# and n= 50
# of lambdaMnD
####################################

####
# note that here we create x by multiplying by (i+1) 
# to avoid computing the IICR_d at x=0 (infinite...)
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

axes()

arrow(0.1, 290, -0.05, (lambdaMnD(0.05, 10, 50) - 290)+30, width=0.001, head_width=0.01, head_length=24, color=color, linestyle="dotted")
arrow(0.05, 470, -0.02, lambdaMnD(0.03, 20, 50) - 470+30, width=0.001, head_width=0.01, head_length=25, color=color, linestyle="dotted")
arrow(-0.03, -10, 0.04, lambdaMnD(0.01, 50, 50)+10-60, width=0.001, head_width=0.01, head_length=15, color=color, linestyle="dotted")

title("$IICR: \lambda_d(t)$ (genes sampled in different demes)")
axis([-0.15, 1, -100, 1500])
xlabel("scaled time ($\t{t=T/N}$)")
ylabel("$IICR: \lambda_d(t)$")

savefig("Fig_IICRd_n50varM_50_05.png")


################################################
##########   figure with subplots  #############
# Plot with variable M  
# and n= 50
# of lambdaMnS
# with low and high migration rates
# This is just to put the two plots above into 
# one single figure
################################################


fig = plt.figure(figsize=(6,9))
fontsize = 8
#fig.tick_params(axis='both', labelsize=10)

### First subplot

ax1 = fig.add_subplot(211)
#ax1 = fig.add_subplot(311)

x = [0.01*i for i in range(601)] # vector for the x-axis
n=50

color = "black"


ax1.axhline(y=50, xmin=0, xmax=5, linestyle="--", color =color)

for M in [50, 20, 10, 5, 2, 1, 0.5] :
    y = []
    for i in range(len(x)) :
        y.append(lambdaMnS(x[i], M, n))
    ax1.plot(x,y, color=color)

text(4.7, 130, "M=0.5", fontsize=fontsize)
text(3.5, 95, "M=1.0", fontsize=fontsize)
text(2.5, 75, "M=2", fontsize=fontsize)
text(1.4, 64, "M=5", fontsize=fontsize)
text(0.6, 57, "M=10", fontsize=fontsize)
text(0.08, 53, "M=20", fontsize=fontsize)
text(-0.50, 42, "M=50", fontsize=fontsize)

title("$(a) IICR: \lambda_s(t)$ (n=50, variable migration rates)", fontsize=1.2*fontsize)
axis([min(x)-max(x)/10, max(x)*1.1, min(y)-max(y)/10, max(y)*1.1])
xlabel("scaled time ($\t{t=T/N}$)", fontsize=fontsize*1.1)
ylabel("$IICR: \lambda_s(t)$", fontsize=fontsize*1.1)

### Second subplot

ax2 = fig.add_subplot(212)
#ax2 = fig.add_subplot(313)

x = [0.001 * i for i in range(1001)]
color = "black"
fontsize=8

for M in [100, 500] :
    for n in [50, 100, 200, 400, 500] :
        y = []
        for i in range(len(x)) :
            y.append(lambdaMnS(x[i], M, n))
        if M==100 :
	  ax2.plot(x,y, color=color, linestyle=":")
	else :
	  ax2.plot(x,y, color=color)
        ax2.axhline(y=n, xmin=0, xmax=5, linestyle="--", color =color)

## Add text to figure
text(0.1, 450, "n=500, M=100", fontsize=fontsize)
text(0.1, 350, "n=400, M=100", fontsize=fontsize)
text(0.1, 150, "n=200, M=100", fontsize=fontsize)
text(0.1, 70, "n=100, M=100", fontsize=fontsize)
text(0.1, 20, "n=50, M=100", fontsize=fontsize)

text(-0.17, 450, "n=500, M=500", fontsize=fontsize)
text(-0.17, 350, "n=400, M=500", fontsize=fontsize)
text(-0.17, 150, "n=200, M=500", fontsize=fontsize)
text(-0.17, 70, "n=100, M=500", fontsize=fontsize)
text(-0.17, 20, "n=50, M=500", fontsize=fontsize)

## Prepare title, axes, labels etc.
title("$(b) IICR: \lambda_s(t)$ (strong migration)", fontsize=fontsize*1.2)
axis([min(x)-max(x)/5, max(x)*0.6, min(y)-max(y)/10, max(y)*1.15])
xlabel("scaled time ($\t{t=T/N}$)", fontsize=fontsize*1.1)
ylabel("$IICR: \lambda_s(t)$", fontsize=fontsize*1.1)

#fig.subplots_adjust(bottom=0.5)
plt.tight_layout()

savefig("Fig_IICRs_a_b.png", dpi=600)

