# object oriented sun earth simulation
#%%
import numpy as np
# import pandas as pd
import spiceypy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import datetime
import matplotlib.animation as animation
# from astroquery.jplhorizons import Horizons
# from IPython.display import HTML
from IPython.display import display
import ipywidgets as widgets

spiceypy.furnsh("../_kernels/pck/gm_de431.tpc")
spiceypy.furnsh("../_kernels/lsk/naif0012.tls")
spiceypy.furnsh("../_kernels/spk/de432s.bsp")
# consts
start_date_utc = datetime.datetime(year=2000, month=1, day=1)
end_date_utc = datetime.datetime(year=2000, month=9, day=1)
start_date_et = spiceypy.utc2et(start_date_utc.strftime("%Y-%m-%dT00:00:00"))
end_date_et = spiceypy.utc2et(end_date_utc.strftime("%Y-%m-%dT00:00:00"))

year = 365*86400                            # sec in year [s]
au = 1.49598e11                             # AU in m     [m]
duration = (end_date_utc-start_date_utc).days*86400
t0 = 0.
nt = int(np.round(duration/86400,0))
t = np.linspace(t0,duration,nt)
legend_objects = []
legend_titles = []

class Planet:

    data_skip = 1

    def __init__(self, naifid, name, orbiting):
        self.name = name
        self.barrycenter_id = int(str(naifid)[0])
        self.state, self.r_sun = spiceypy.spkgeo(targ=self.barrycenter_id, et=start_date_et, \
                                    ref="ECLIPJ2000", obs=10)
        self.state = self.state*1000                        # km -> m
        self.parent = orbiting
        self.y = []
        self.line = []
        orbiting.add_planet(self)
    
    def f(self, y, t, mstar):
        self.gm = solar_system.G*mstar
        self.x = y[0:3]
        self.v = y[3:]
        self.r = np.linalg.norm(self.x)
        self.dxdt = self.v
        self.dvdt = -self.gm*self.x/self.r**3
        self.dy = np.hstack((self.dxdt, self.dvdt))
        return self.dy
    
    def plot_orbit_2d(self):
        y = self.calculate_orbit()                      
        self.y = y
        self.line, = space_map_ax.plot(y[:,0],y[:,1])
        self.line.set_label(self.name)
        legend_objects.append(self.line)
        legend_titles.append(self.name)
        return self.line, self.y

    def calculate_orbit(self):
        y = odeint(self.f, self.state,t, args=(self.parent.sun_mass,))
        y = y/au                                        # divide by AU 
        return y
    
    def animate_orbit(self):
        for i in range(len(self.y)):
            space_map_anim[0].plot(self.y[i:i+Planet.data_skip,0]/au, \
                self.y[i:i+Planet.data_skip,1]/au, "o", color="blue")


    @classmethod
    def from_string(cls, planet_str):
        id, name, orbiting = planet_str.split('-')
        return cls(id, name, orbiting)

    
class solar_system:

    G = 6.67428e-11                             # m3 / kg / s^2
    def __init__(self, name, sun_mass):
        self.name = name
        self.sun_mass = sun_mass
        self.planets = []
        self.ax = []
        self.fig = []

    
    def add_planet(self, planet):
        self.planets.append(planet)

    def create_space_map(self):
        plt.style.use("dark_background")
        fig, ax = plt.subplots(figsize=(8,6))
        ax.grid(linestyle="dashed", alpha=0.4)
        ax.set_title(self.name + " Solar System")
        ax.plot([0.],[0.],'o', color="tab:orange")
        ax.set_aspect("equal")
        ax.set_ylabel("Distance in AU")
        ax.set_xlabel("Distance in AU")
        return ax, fig
    
    def populate_ss(self):
        for i in self.planets:
            i.plot_orbit_2d()


# Let's add our solar system as an instance

Sun = solar_system("Sun", 1.98892e30)

# create our solar system and save Axes object as space_map
space_map_ax, space_map_fig = Sun.create_space_map()


mercury = Planet(199, "Mercury", Sun)
venus = Planet(299, "Venus", Sun)
earth = Planet(399, "Earth", Sun)
mars = Planet(499, "Mars", Sun)
# jupiter = Planet(599, "Jupiter", Sun)
# saturn = Planet(699, "Saturn", Sun)
# uranus = Planet(799, "Uranus", Sun)
# neptune = Planet(899, "Neptune", Sun)

Sun.populate_ss()

space_map_ax.legend(legend_objects, legend_titles)
# plt.savefig("hei.png")
plt.show()
plt.close('all')

# space_map_anim = Sun.create_space_map()
# anim = animation.FuncAnimation(space_map_anim[1], earth.animate_orbit(),  \
#             frames=duration, interval=200)
# anim.save('ani.mp4',writer='ffmpeg',fps=1000/50)


#%%

# kind of functioning. Here we have no animation, but no errors either...

fig = plt.figure(2)
ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
scat = ax.scatter([], [], s=60)

def init():
    scat.set_offsets([])
    return scat,

def animate(i):
    data = [earth.y[i,0], earth.y[i,1]]
    scat.set_offsets(data)
    return scat,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(earth.y), 
                               interval=20, blit=False, repeat=False)

anim.save('animation.mp4')
# lns = []
# plots = []

# def update(i):
#     ln.set_xdata(y[i,0]/au)
#     ln.set_ydata(y[i,1]/au)
#     plots.set_offsets
    
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


# fig, ax = plt.subplots(figsize=(8,6))
# line = ax.plot([], [], color="blue")
# data_skip = 600
# print(line)
# # initialization function: plot the background of each frame
# def init():
#     ax.clear()
#     ax.set_aspect("equal")
#     ax.set_xlim([-1.2, 1.2])
#     ax.set_ylim([-1.2, 1.2])
#     ax.plot(y[:,0]/au, y[:,1]/au, color="blue")

# # animation function.  This is called sequentially
# def animate(i):
    # x = np.linspace(0, 2, 1000)
    # y = np.sin(2 * np.pi * (x - 0.01 * i))
    # line.set_data(y[i:i+data_skip,0]/au, y[i+data_skip,1]/au)
    # ax.plot(y[i:i+data_skip,0]/au, y[i:i+data_skip,1]/au, "o", color="blue")
    # ax.plot(y[i+data_skip,0]/au, y[i+data_skip,1]/au, 'o')
    # ax.plot(0., 0., 'o')
    # line.set_data(x, y)
    #return line,

# anim = animation.FuncAnimation(fig, animate, init_func=init,  \
#                                frames=np.arange(0, len(y)-data_skip-1, data_skip), interval=20)

# anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

# plt.show()

# %%

# HTML(anim.to_html5_video())
# %%
