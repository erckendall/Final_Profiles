# Note: Because we are working with a cubic grid, performing a spherical average is subject to some considerable limitations.
# Firstly, the total size of the grid in terms of resolution determines the total number of different directions which can be probed. i.e. a grid with 3 grid points
# in each of the three dimensions can only be probed in 13 directions through the centre. As the grid size increases, more directions can be probed,
# corresponding to more pairs of points which draw a line through the centre.
# (Though in a cubic grid there will be separate pairs of points lying along the same line in some directions in general)
# Because of the spacing, some directions will not be resolved as well as others. The three primary directions will be resolved best,
# and then the 45 degree angled lines will also be appreciably resolved in general. For other directions the resolution will get progressively worse.
# One could perform a spherical interpolation, to try to resolve the other directions better, but would this yield better results than an average from
# the 13 main directions? Maybe, but also much more difficult computationally.
# And the diagonals are not all equal in 3d - 3 different resolutions
# Also we require a centre, which is unnatural on a cubic grid with even-valued resol.


import numpy as np
import matplotlib.pyplot as plt

fc = np.load('rho_#400.npy')
resol = np.shape(fc)[0]
gridlength = 10.
trunc = 90
gridspace = float(gridlength/resol)
est_cen = np.unravel_index(np.argmax(fc), fc.shape)
max_den = fc[est_cen]
rge = min(est_cen)
sol = []
ma = 10**(-22)
print est_cen

sm = min(resol-est_cen[0],est_cen[0])
sm = min(sm, resol-est_cen[1], est_cen[1], est_cen[2], resol-est_cen[2])

print (est_cen)


f = np.load('initial_f.npy')
dr = .00001
alpha = np.sqrt(max_den)
rhm = (1/np.sqrt(alpha))*0.69
mcore = np.sqrt(alpha)*0.918959850533
for i in np.arange(trunc):
    if (int(np.sqrt(alpha) * (i*gridspace / dr + 1))) < 900000:
        sol.append((alpha * f[int(np.sqrt(alpha) * (i*gridspace / dr + 1))])**2)
    else:
        sol.append(0)

north = []
east = []
south = []
west = []
up = []
down = []

for i in np.arange(sm):
    north.append(fc[est_cen[0],est_cen[1]+i, est_cen[2]])
    south.append(fc[est_cen[0],est_cen[1]-i, est_cen[2]])
    east.append(fc[est_cen[0]+i,est_cen[1], est_cen[2]])
    west.append(fc[est_cen[0]-i,est_cen[1], est_cen[2]])
    up.append(fc[est_cen[0], est_cen[1], est_cen[2]+i])
    down.append(fc[est_cen[0], est_cen[1], est_cen[2]-i])



n_trunc = []
s_trunc = []
e_trunc = []
w_trunc = []
u_trunc = []
d_trunc = []

avg = []
data = []
nfw = []
fitting = []

def fit(r):
    return (4.807*10**7*(ma/10**(-23))**(-2)*(rhm/0.02613)**(-4))/(1+0.091*(r/rhm)**2)**8


for i in np.arange(trunc):
    n_trunc.append(north[i])
    data.append(north[i])
    s_trunc.append(south[i])
    data.append(south[i])
    e_trunc.append(east[i])
    data.append(east[i])
    w_trunc.append(west[i])
    data.append(west[i])
    u_trunc.append(up[i])
    data.append(up[i])
    d_trunc.append(down[i])
    data.append(down[i])
    avg.append(np.average(data))
    data = []
    fitting.append(fit(i*gridspace))
    if i == 0:
        nfw.append(0)
    else:
        nfw.append(1/(i*gridspace)**3)

# init = avg[0]
#
# for i in np.arange(trunc):
#     rat = avg[i]/avg[0]
#     if rat <= 0.5:
#         print rat
#         print i
#         break


gradf = (avg[trunc-1]-avg[trunc-2])/gridspace
rf = trunc*gridspace
rhosrs3 = rf**4*gradf/(-3)

const = rf**3*avg[trunc-1]




# test = []
# for i in np.arange(trunc):
#     if i == 0:
#         test.append(0)
#     else:
#         test.append(rhosrs3/(i*gridspace)**3)

# plt.loglog(n_trunc, label='N')
# plt.loglog(s_trunc, label = 'S')
# plt.loglog(e_trunc, label = 'E')
# plt.loglog(w_trunc, label = 'W')
plt.loglog(avg, label='Avg.')
plt.loglog(sol, label='soliton')
# plt.loglog(nfw, label='1/r^3')
# plt.loglog(test, label='test')
plt.loglog(fitting, label='fit', linestyle=':')
plt.axes().axvline(x=(rhm/gridspace), color='r', label='r_c')
plt.axes().axvline(x=3*(rhm/gridspace), color='g', label='3*r_c')
plt.legend()
plt.xlim(1,trunc-1)
plt.ylim(1,max_den)
plt.show()



vol = 0
su = 0
cont = 0

vir = False
for i in np.arange(trunc):
    radius = gridspace*i
    cont = cont + avg[i]*((4./3.)*np.pi*((radius+gridspace)**3.-radius**3.))
    su = su + ((4./3.)*np.pi*((radius+gridspace)**3.-radius**3.))
    avdense = cont / su
    if avdense <= 374:
        print ('{}{}'.format('half-density radius = ', rhm))
        print('{}{}'.format('Virial radius is ', radius))
        rf = radius
        vm = 4/3*np.pi*radius**3*avdense
        print('{}{}'.format('Virial mass is ', vm))
        vir = True
        break

if vir == False:
    print 'Virial radius outside grid boundary'

while vir == False:
    rf = rf + gridspace
    avg.append(const/rf**3)
    cont = cont + (const/rf**3)*((4./3.)*np.pi*((rf+gridspace)**3.-rf**3.))
    su = su + ((4./3.)*np.pi*((rf+gridspace)**3.-rf**3.))
    avdense = cont / su
    # if rf >= rhm:
    #     print ('{}{}'.format('Core mass est is ', 4/3*np.pi*rf**3*rat))
    if avdense <= 374:
        vir = True
        print ('{}{}'.format('half-density radius = ', rhm))
        print('{}{}'.format('Virial radius is ', rf))
        vm = 4/3*np.pi*rf**3*avdense
        print('{}{}'.format('Virial mass is ', vm))




mass = np.load('masslist.npy')
print('{}{}'.format('Total initial mass = ', mass[0]))
print('{}{}'.format('Total final mass = ', mass[400]))
egy = np.load('egylist.npy')

mcoretest = (3./20.)**(1./2.)*(3./(4.*np.pi*374))**(-1./6.)*vm**(1./3.)
mcoretest2 = (3*vm/(20*rf))**0.5
print ('{}{}'.format('Calculated core mass is ',mcoretest))
print ('{}{}'.format('Actual core mass is ',mcore))
print ('{}{}'.format('Mc from mh and rv ', mcoretest2))
