# Change Log

##Closed issues
- Different seeds for CMB and noise.
- Code crashes if debug is not passed.
- Nside of different templates. 

##Changes made
- Changed curved synchrotron power law option to follow Kogut 2012 curvature model. (29.6.16)
- Change free-free emission to be a -2.17 power law rather than analytical for from Commander. (29.6.16)  

- Added smoothing with a Gaussian kernel after summation of all components. (29.6.16)
- Cleared up file saving code for clarity.(29.6.16)
- Now saves main_config information to the header of the output fits file. (29.6.16)

- Added multithreading to main.py to allow component computation to be done faster. (29.6.16)

**Merged pull requests:**

- made debug optional with default 'False' [\#6](https://github.com/bthorne93/PySM/pull/6) ([damonge](https://github.com/damonge))

