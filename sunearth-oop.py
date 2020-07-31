#%% [markdown]
# # Solar system simulator
# 
# By Ståle Gjelsten.
#
# This script was written in the summer of 2020 to improve on a program I had written for
# [pico-8](https://www.lexaloffle.com/pico-8.php). Because of the pico-8's limited float and
# integer sizes (and the fact that I had to use [Euler method](https://en.wikipedia.org/wiki/Euler_method) 
# for integration), the orbit calculations were way off. This script is by no means an accurate representation
# of a planetary system -- it only calculates forces between each satelite and the star -- but at least
# it seems to be fairly close to other simulations and measurements.
#
# ## Theory on satelite orbits
# 
# The gravitational force $\vec{F}$ on a satelite with mass $m$ is given by:
# $$ {F}_G = - \frac{GMm}{r^2} $$
# Where $G$ is the gravitational constant, $M$ is the mass of the central object (the Sun in our case), 
# $m$ is the mass of the satelite and $r^2$ is the square of the distance between the objects.
# 
# We can also express the gravitational force as a vector:
# $$ \vec{F}_G = - \frac{GMm}{r^3}\vec{r} $$
# where $\vec{r}$ is the position vector and $r= \lvert x \rvert$.
#
# From Newtons second law we have $\sum{\vec{F}} = m\vec{a}$ and substituting for $\vec{F}_G$ gives:
# $$ m\vec{a} = - \frac{GMm}{r^3}\vec{r} \Leftrightarrow \vec{a} = - \frac{GM}{r^3}\vec{r} $$
#
# We can find the velocity and position of the satelite by integration as the acceleration of the satelite
# is the derivative of the velocity and the double derivative of the position.
# $$  \vec{a} = \vec{\dot v} = \vec{\ddot r} $$
#
# In this script we find the positions by numerical integration. To accomplish that we construct a vector $\vec{x}(t)$ such that:
# $$ \vec{x}(t) = \begin{bmatrix} x(t)\\ y(t) \\ z(t) \\ v_x(t) \\ v_y(t) \\ v_z(t) \end{bmatrix} $$
# This state vector contains both the position and velocity of the satelite.
# We also define the function $f$ such that $\frac{\mathrm{d}\vec{x}}{\mathrm{d}t} = f(\vec{x},t)$ and use numerical integration 
# to find the solution to our differential equation.
# $$ f(\vec{x},t) = \begin{bmatrix} v_x(t)\\ v_y(t)\\ v_z(t)\\ F_x/m\\ F_y/m\\ F_z/m \end{bmatrix} $$
#
# The initial values for our initial value problem are sourced from [NASAs SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html).
#
# Suggested reading on python simulation of the two-body problem:
#
# - [Scientific computing toolbox's introduction](https://faculty1.coloradocollege.edu/~sburns/toolbox/ODE_II.html) to 
# `odeint()` in python
# - [C.P. Dullemond of the University of Heidelberg's introduction to the N-body model](http://www.ita.uni-heidelberg.de/~dullemond/lectures/studtage_compastro_2018/Chapter_1.pdf)
# - [University of Oslo AST1100 lecture notes on the general solution to the two-body problem](https://www.uio.no/studier/emner/matnat/astro/AST1100/h13/undervisningsmateriale/ast1100-fullstendig.pdf)

#%%
import numpy as np
import spiceypy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import datetime
import matplotlib.animation as animation
# from IPython.display import HTML

# loading with kernels with planetary data for NASAs SPICE module
kernel_dir = "../_kernels/"        # directory with SPICE kernels

spiceypy.furnsh(kernel_dir + "pck/gm_de431.tpc")
spiceypy.furnsh(kernel_dir + "lsk/naif0012.tls")
spiceypy.furnsh(kernel_dir + "spk/de432s.bsp")
spiceypy.furnsh(kernel_dir + "pck/pck00010.tpc")

# figure variables
figure_size = (6, 4)                        # a tuple with size of animation
                                            # figure in inches

# time variables
start_date_utc = datetime.datetime(year=2000, month=1, day=1)
end_date_utc = datetime.datetime(year=2002, month=3, day=1)
sec_per_iteration = 86400                   # seconds per iteration of the
                                            # calculations and animation

# "universal" constants
year = 365*86400                            # sec in year [s]
au = 1.49598e11                             # AU in m     [m]


start_date_et = spiceypy.utc2et(start_date_utc.strftime("%Y-%m-%dT00:00:00"))
end_date_et = spiceypy.utc2et(end_date_utc.strftime("%Y-%m-%dT00:00:00"))
duration = (end_date_utc - start_date_utc).days * 86400  # sim. duration in [s]
nt = int(np.round(duration / sec_per_iteration, 0))      # number of steps
t = np.linspace(0, duration, nt)                         # evenly spaced time
                                                         #   array for calcs 

# TODO: consider moving these to the solar system class

legend_objects = []                                     
legend_titles = []

class Planet:
    # all planets are stored in the Planet class. 

    def __init__(self, naifid, name, orbiting):
        # planets are constructed with their NASA NAIF ID's (int), name (str)
        # and the body they are orbiting (planetary_system instance) 

        self.name = name
        self.naifid = naifid

        # not all planets' states are available in spiceypy.spkgeo, so we
        # use their barrycenters' to get state vectors

        self.barrycenter_id = int(str(naifid)[0])
        self.state, self.r_sun = spiceypy.spkgeo(targ=self.barrycenter_id, \
                                                 et=start_date_et, \
                                                 ref="ECLIPJ2000", obs=10)
        
        _, radii = spiceypy.bodvcd(naifid, "RADII",3)
        self.radii = np.average(radii)*1000                 # km -> m conv.
        self.vis_area = 2 * self.radii ** 2 * np.pi
        self.size = int(np.ceil(self.vis_area/orbiting.scaling_factor*50000))
        self.state = self.state*1000                        # km -> m conv.
        self.parent = orbiting
        self.gm = orbiting.G*orbiting.sun_mass
        self.y = []
        self.line = []
        self.anim_data = []
        self.scat = []
        self.color = []
        orbiting.add_planet(self)
    
    def f(self, y, t):
        x = y[0:3]
        v = y[3:]
        r = np.linalg.norm(x)
        dxdt = v
        dvdt = -self.gm*x/r**3
        self.dy = np.hstack((dxdt, dvdt))
        return self.dy
    
    def plot_orbit_2d(self):
        self.y = self.calculate_orbit()                      
        self.line, = self.parent.ax.plot(self.y[:,0],self.y[:,1], \
                                         color=self.color)
        self.line.set_label(self.name)
        legend_objects.append(self.line)
        legend_titles.append(self.name)
        return self.line, self.y

    def calculate_orbit(self):
        y = odeint(self.f, self.state,t)
        y = y/au                                        # divide by AU 
        return y
    
    def make_scatter_for_animation(self):

        x = self.y[:,0]
        y = self.y[:,1]
        self.data = np.hstack((x[:,np.newaxis], y[:, np.newaxis]))

        return self.data
    
    def animate_orbit(self):
        data = self.make_scatter_for_animation()
        return data


    @classmethod
    def from_string(cls, planet_str):
        id, name, orbiting = planet_str.split('-')
        return cls(id, name, orbiting)

    
class planetary_system:

    G = 6.67428e-11                             # m3 / kg / s^2
    def __init__(self, naifid, name, sun_mass):
        self.naifid = naifid
        self.name = name
        self.sun_mass = sun_mass
        self.planets = []
        self.num_planets = []
        _, radii = spiceypy.bodvcd(naifid, "RADII",3)
        self.radii = np.average(radii)*1000                 # km -> m conv.
        self.vis_area = 2 * self.radii ** 2 * np.pi
        self.size = 20
        self.scaling_factor = self.vis_area/self.size
        self.ax = []
        self.fig = []
        self.data = []
        self.scatters = []

    
    def add_planet(self, planet):
        self.planets.append(planet)

    def create_space_map(self):
        plt.style.use("dark_background")
        fig, ax = plt.subplots(figsize=figure_size)
        ax.grid(linestyle="dashed", alpha=0.4)
        ax.set_title(self.name + " Solar System")
        ax.plot([0.],[0.], marker='o', markersize=self.size, color="tab:orange")
        ax.set_aspect("equal")
        ax.set_ylabel("Distance in AU")
        ax.set_xlabel("Distance in AU")
        return ax, fig
    
    def populate_ss(self):

        cmap = plt.get_cmap('plasma')
        colors = cmap(np.linspace(0, 1, len(self.planets)))
        j = 0
        self.num_planets = len(self.planets)

        for i in self.planets:
            i.color = colors[j]
            i.plot_orbit_2d()
            j+=1

    def animate_ss(self):

        self.num_planets = len(self.planets)
        self.data = np.zeros((nt, self.num_planets, 2))
        self.color = np.zeros((nt, self.num_planets, 4))
        
        planet_no = 0
        
        for i in self.planets:
            data   = i.animate_orbit()
            for j in range(len(data)):
                self.color[j,planet_no,:] = i.color
                for k in range(2):
                    self.data[j,planet_no,k] = data[j,k]
            planet_no+=1

        def update_plot(i, data, scat):
            for j in range(self.num_planets):
                scat[j].set_offsets(data[i,j])
            time = start_date_utc + datetime.timedelta(days=i)
            time_str = time.strftime("%Y-%m-%d")
            timer.set_text("Date: " + time_str)
            return scat,timer

        timer = self.ax.text(0,0, "", ha="left", va="bottom", \
                             transform=self.ax.transAxes)
                            
        for i in range(self.num_planets):
            self.scatters.append(self.ax.scatter([],[], s=self.planets[i].size,\
                                 color=self.planets[i].color))
        
        ani = animation.FuncAnimation(self.fig, update_plot, frames=nt-1,\
            interval=33, fargs=(self.data, self.scatters))
        plt.show()
        return ani
        


# Let's add our solar system as an instance
Sun = planetary_system(10, "Sun", 1.98892e30)

# create our solar system and save Axes object as space_map
Sun.ax, Sun.fig = Sun.create_space_map()


# Let's add some planets
mercury = Planet(199, "Mercury", Sun)
venus = Planet(299, "Venus", Sun)
earth = Planet(399, "Earth", Sun)
mars = Planet(499, "Mars", Sun)
# jupiter = Planet(599, "Jupiter", Sun)
# saturn = Planet(699, "Saturn", Sun)
# uranus = Planet(799, "Uranus", Sun)
# neptune = Planet(899, "Neptune", Sun)

# Let's populate the solar solar system with our planets...
Sun.populate_ss()
Sun.ax.legend(legend_objects, legend_titles)

# ... and start the animation. Comment the next line if you just
# want a plot with the planetary orbits
animation = Sun.animate_ss()

# animation.save("solarsystem.gif", writer="imagemagick", fps=30)
# plt.savefig("solar_system.png")

# plt.close('all')
