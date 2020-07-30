# Solar system animation

This repo contains a python script to simulate solar systems. The script has only been tested with our own [Solar System](https://en.wikipedia.org/wiki/Solar_System).

![Solar system animation example](solarsystem-sim-example.gif)

## Requirements

The script requires the following modules:

- datetime
- spiceypy
- scipy
- numpy
- matplotlib

## How the script works

Enter the simulation start and end dates in the start_date_utc variables. The sec_per_iteration variable controls the precision of the calculations and the speed of the animation. The default value of 86 400 is equal to one day (60 sec/min \* 60 min/hour \* 24 hours/day), which means we perform one calculation of the orbits per day and that every frame of the animation represents one day.

Initialize a planetary system by doing something like `Sun = planetary_system(10, "Sun", sun_mass_in_kg)`.

Initialize planets as `Planet` instances by doing something like `earth = Planet(399, "Earth", Sun)`, where Sun is the planetary_system instance you constructed earlier. 399 is the [NASA NAIF ID](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html) for the planet.

TODO: write more text here... 😅
