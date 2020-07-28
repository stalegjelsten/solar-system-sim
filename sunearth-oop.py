# object oriented sun earth simulation
#%%
import numpy as np
# import pandas as pd
import spiceypy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astroquery.jplhorizons import Horizons
from IPython.display import HTML

# consts
G = 6.67428e-11                             # m3 / kg / s^2
msun = 1.98892e30                           # solar mass [kg]
year = 365*86400
au = 1.49598e11                             # AU in m
start_date = "2020-01-01"
duration = 1*year



def f(y,t,mstar):
    gm = G*mstar
    x = y[0:3]
    v = y[3:6]
    r = np.linalg.norm(x)
    dxdt = v
    dvdt = -gm*x/r**3
    dy = np.hstack((dxdt, dvdt))
    return dy

t0 = 0.
nt = int(np.round(duration/365,0))
t = np.linspace(t0,duration,nt)
r0 = au                             # initial dist sun earth
vp0 = 2*np.pi*au/year               # initial earth velocity
x0 = np.array([r0, 0., 0.])
v0 = np.array([0., vp0, 0.])
y0 = np.hstack((x0,v0))


y = odeint(f,y0,t,args=(msun,))




# fig, ax = plt.subplots(figsize=(8,6))
# ax.plot(y[:,0]/au,y[:,1]/au)
# ax.plot([0.],[0.],'o')
# ax.set_aspect("equal")

# fig = plt.figure(figsize=(5, 5))
# ax = fig.add_subplot(1, 1, 1)

# lns = []
# for i in range(len(y)):
#     ln = ax.plot([y[i,:2]/au], color="tab:orange", lw=2)
#     # ln, = ax.plot([0, np.sin(sol[i, 0])], [0, -np.cos(sol[i, 0])],
#     #               color='k', lw=2)
#     lns.append(ln)
# print(type(lns[2]))
# print(lns[-3])
# ax.set_aspect('equal', 'datalim')
# ax.grid()
# ani = animation.ArtistAnimation(fig, lns, interval=200)
# ani.save('ani.mp4',writer='ffmpeg',fps=1000/50)

# First set up the figure, the axis, and the plot element we want to animate
# fig = plt.figure()
# ax = plt.axes(xlim=(-2, 2), ylim=(-2, 2))


fig, ax = plt.subplots(figsize=(8,6))
line = ax.plot([], [], color="blue")
data_skip = 600
print(line)
# initialization function: plot the background of each frame
def init():
    ax.clear()
    ax.set_aspect("equal")
    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.plot(y[:,0]/au, y[:,1]/au, color="blue")

# animation function.  This is called sequentially
def animate(i):
    # x = np.linspace(0, 2, 1000)
    # y = np.sin(2 * np.pi * (x - 0.01 * i))
    # line.set_data(y[i:i+data_skip,0]/au, y[i+data_skip,1]/au)
    ax.plot(y[i:i+data_skip,0]/au, y[i:i+data_skip,1]/au, "o", color="blue")
    # ax.plot(y[i+data_skip,0]/au, y[i+data_skip,1]/au, 'o')
    # ax.plot(0., 0., 'o')
    # line.set_data(x, y)
    #return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,  \
                               frames=np.arange(0, len(y)-data_skip-1, data_skip), interval=20)

anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

# plt.show()

class Planet:
    def __init__(self, naifid, name, mass, state_vector, velocity_vector, \
        acc_vector):
        self.id = naifid
        self.name = name
        self.mass = mass
        self.state = state_vector
        self.velocity = velocity_vector
        self.accel = acc_vector
    


# %%

HTML(anim.to_html5_video())
# %%
