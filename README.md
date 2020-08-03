# Solar system simulator

This repo contains a python script to simulate solar systems. The script has only been tested with our own [Solar System](https://en.wikipedia.org/wiki/Solar_System).

![Solar system animation example](solarsystem-sim-example.gif)

This script was written in the summer of 2020 to improve on a program I had written for [pico-8](https://www.lexaloffle.com/pico-8.php). Because of the pico-8's limited float and integer sizes (and the fact that I had to use [Euler method](https://en.wikipedia.org/wiki/Euler_method) for integration), the orbit calculations in that program were way off. *This script* is by no means an accurate representation of a planetary system -- it only calculates forces between each planet and the star -- but at least it seems to be fairly close to other simulations and measurements.

## Requirements

The script requires the following modules:

- [datetime](https://docs.python.org/3/library/datetime.html#module-datetime) for date calculations.
- [spiceypy](https://pypi.org/project/spiceypy/) for initial planetary orbit positions.
- [scipy](https://pypi.org/project/scipy/) for `odeint()` function to solve differential equations.
- [numpy](https://pypi.org/project/numpy/) because math. And ndarrays.
- [matplotlib](https://pypi.org/project/matplotlib/) to plot and animate the orbits.

## How the script works

Enter the simulation start and end dates in the start_date_utc variables. The sec_per_iteration variable controls the precision of the calculations and the speed of the animation. The default value of 86 400 is equal to one day (60 sec/min \* 60 min/hour \* 24 hours/day), which means we perform one calculation of the orbits per day and that every frame of the animation represents one day.

Initialize a planetary system by doing something like `Sun = planetary_system(10, "Sun", sun_mass_in_kg)`.

Initialize planets as `Planet` instances by doing something like `earth = Planet(399, "Earth", Sun)`, where Sun is the planetary_system instance you constructed earlier. 399 is the [NASA NAIF ID](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html) for the planet.

## How accurate is the simulation?

I have not tested the accuracy of the simulation yet, but the orbits seem to be quite accurate (and Earth takes 1 year to orbit the Sun). *The size of the planets and the Sun is not to scale*. The sizes of the planets are correct relative to each other, but they are magnified relative to the Sun. The size of the planets and the Sun is also magnified relative to the distances between them.

## Theory on planetary orbits

[Planetary orbits](https://en.wikipedia.org/wiki/Orbit#Planetary_orbits) are elliptical and orbit the system's [barycenter](https://en.wikipedia.org/wiki/Barycenter).

The [gravitational force](https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation) <img alt="$\vec{F}$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/1e542f84ff6a79def15ac7917bec74ec.svg?137b59b5ef" align="middle" width="13.17075374999999pt" height="31.799054100000024pt"/> on a planet with mass <img alt="$m$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/0e51a2dede42189d77627c4d742822c3.svg?ca97e77fbd" align="middle" width="14.433101099999991pt" height="14.15524440000002pt"/> is given by:

<p align="center"><img alt="$$ {F}_G = - \frac{GMm}{r^2} $$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/f755c890b14bdcd14a32bf24b2d8785b.svg?9e703b6cf3" align="middle" width="103.40032889999999pt" height="33.62942055pt"/></p>

Where <img alt="$G$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/5201385589993766eea584cd3aa6fa13.svg?c238b60897" align="middle" width="12.92464304999999pt" height="22.465723500000017pt"/> is the [gravitational constant](https://en.wikipedia.org/wiki/Gravitational_constant), <img alt="$M$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/fb97d38bcc19230b0acd442e17db879c.svg?18121ab4d8" align="middle" width="17.73973739999999pt" height="22.465723500000017pt"/> is the mass of the central object (the Sun in our case), <img alt="$m$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/0e51a2dede42189d77627c4d742822c3.svg?788ce23815" align="middle" width="14.433101099999991pt" height="14.15524440000002pt"/> is the mass of the planet and <img alt="$r^2$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/dd9ad1899e5c8220c8b4bbc13483d097.svg?910ee47a69" align="middle" width="14.42550119999999pt" height="26.76175259999998pt"/> is the square of the distance between the objects.  

We can also express the [gravitational force as a vector](https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation#Vector_form):

<p align="center"><img alt="$$ \vec{F}_G =  - \frac{GMm}{r^2} \cdot \frac{\vec{r}}{\lvert \vec{r} \rvert} = - \frac{GMm}{r^3}\vec{r} $$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/c564886d8e04a648668f3b49e599e217.svg?9035ecf466" align="middle" width="232.1647944pt" height="38.2504353pt"/></p>

where <img alt="$\vec{r}$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/b32c2a51de7b6df016e08d3c668bdf29.svg?e334d2c176" align="middle" width="10.747741949999993pt" height="23.488575000000026pt"/> is the position vector and <img alt="$r$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/89f2e0d2d24bcf44db73aab8fc03252c.svg?c5de3dfa52" align="middle" width="7.87295519999999pt" height="14.15524440000002pt"/> is the norm of <img alt="$\vec{r}$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/b32c2a51de7b6df016e08d3c668bdf29.svg?be5d10b922" align="middle" width="10.747741949999993pt" height="23.488575000000026pt"/> such that <img alt="$r= \lvert \vec{r} \rvert$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/326ae36b164b47913b45e59a4d2d081a.svg?ab85c33761" align="middle" width="46.79596679999999pt" height="24.65753399999998pt"/>.

From Newtons second law we have <img alt="$\sum{\vec{F}} = m\vec{a}$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/fa31f112c857f3891d84b6a758ea1682.svg?2a2fbb2528" align="middle" width="79.01440799999999pt" height="31.799054100000024pt"/> and substituting for <img alt="$\vec{F}_G$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/da14a7649c5f43262297c60caaf7a32b.svg?94d37e51c1" align="middle" width="20.805292199999986pt" height="31.799054100000024pt"/> gives:

<p align="center"><img alt="$$ m\vec{a} = - \frac{GMm}{r^3}\vec{r} \Leftrightarrow \vec{a} = - \frac{GM}{r^3}\vec{r} $$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/988d6b0515674d08a9cf09b469c9ceb7.svg?13ad331fee" align="middle" width="228.53684535pt" height="33.62942055pt"/></p>

We can find the velocity and position of the planet by integration as the acceleration of the planet is the derivative of the velocity and the double derivative of the position.

<p align="center"><img alt="$$  \vec{a} = \vec{\dot v}  = \vec{\ddot r} $$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/bc03515af3e320c46a8903985c3f7f93.svg?c6a23a3e5e" align="middle" width="72.33514695pt" height="15.645186149999999pt"/></p>

In this script we find the positions by numerical integration. To accomplish that we construct a vector <img alt="$\vec{x}(t)$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/ec2a7d3c63d41be3886d013bcc8a3056.svg?414df0910e" align="middle" width="28.33622219999999pt" height="24.65753399999998pt"/> such that:

<p align="center"><img alt="$$ \vec{x}(t) = \begin{bmatrix} x(t)\\ y(t) \\ z(t) \\ v_x(t) \\ v_y(t) \\ v_z(t) \end{bmatrix} $$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/579471fcd94ee308cee6dafb55846f11.svg?4c61809fbc" align="middle" width="107.13755249999998pt" height="118.35734295pt"/></p>

This [state vector](https://en.wikipedia.org/wiki/Orbital_state_vectors) contains both the position and velocity of the planet. We also define the function <img alt="$f$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/190083ef7a1625fbc75f243cffb9c96d.svg?7f1432ce91" align="middle" width="9.81741584999999pt" height="22.831056599999986pt"/> such that <img alt="$\frac{\mathrm{d}\vec{x}}{\mathrm{d}t} = f(\vec{x},t)$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/fcbb73f2b5144107a08cb236000135d6.svg?749e34c1e0" align="middle" width="83.8446774pt" height="29.485239299999986pt"/> and use numerical integration to find the solution to our differential equation.

<p align="center"><img alt="$$ f(\vec{x},t) = \begin{bmatrix} v_x(t)\\ v_y(t)\\ v_z(t)\\ F_x/m\\ F_y/m\\ F_z/m \end{bmatrix} $$" src="https://cdn.jsdelivr.net/gh/stalegjelsten/solar-system-sim@master/svgs/694b80136eb5c63c902f9b61d133b20a.svg?cea62a8d35" align="middle" width="130.57465079999997pt" height="118.35734295pt"/></p>

The initial values for our initial value problem are sourced from [NASAs SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html).

Suggested reading on python simulation of the two-body problem:

- [Scientific computing toolbox's introduction](https://faculty1.coloradocollege.edu/~sburns/toolbox/ODE_II.html) to 
`odeint()` in python
- [C.P. Dullemond of the University of Heidelberg's introduction to the N-body model](http://www.ita.uni-heidelberg.de/~dullemond/lectures/studtage_compastro_2018/Chapter_1.pdf)
- [University of Oslo AST1100 lecture notes on the general solution to the two-body problem](https://www.uio.no/studier/emner/matnat/astro/AST1100/h13/undervisningsmateriale/ast1100-fullstendig.pdf)
