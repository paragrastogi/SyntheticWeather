# SyntheticWeather

This repository contains scripts to create synthetic weather time series from a short weather record of at least one year. The tool (scripts and how they might be implemented in ESP-r and Ladybug) is provisionally named `indra`. I haven't thought up a smarmy, contrived acronym yet but I'm working on it.

While the MATLAB/R scripts are well documented (in two conference papers and my thesis), no further work is being undertaken on that algorithm (call it the 'old' algorithm). I will, hopefully by November 2017, translate the 'old' scripts to Python faithfully. 

The Python scripts in the repository can be called with two methods: `arma` and `gp`. The method `arma` is an almost-completely-faithful translation of the MATLAB scripts<sup>[1]</sup>. The method `gp` is new and experimental - which means that it has not been published or tested extensively yet. This consitutes a significant change in the algorithm - I would not recommend using them just yet (true as of October 2017). 

Like the older work, the newer work will also, most likely, be released with a BSD-3 license. Using the older scripts requires only a valid MATLAB license and R (R is free to download and reuse). While you are free to use the scripts as you please, I am not liable for anything that happens as a result of using my scripts. Like if you accidentally release nuclear missiles in the USA, ruin the ski season in Switzerland, or cause a drought in Scotland.

This work is linked to Parag Rastogi's PhD thesis at the Ecole Polytechnique Federale de Lausanne, EPFL (search on <infoscience.epfl.ch> or <paragrastogi.com>) in Lausanne, Switzerland. Parag is the sole author of these scripts and the thesis was supervised by Professor Marilyne Andersen. The scripts are the intellectual property of EPFL, and their redistribution under the liberal BSD-3 license is with the approval of Prof. Marilyne Andersen.

If you have questions or concerns, or notice errors, please contact me at `contact[at]paragrastogi.com`.

Happy creating fake weather!

### Footnotes

[1] The difference lies in the 'simulation' of the SARMA model to produce synthetic 'de-mean-ed' series. I am not convinced that I should reproduce the 'custom noise' functionality used in the old scripts to simulate the SARMA models with bootstrapped residuals. For now, I am doing a 'conventional' simulation by using white noise.
