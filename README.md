# SyntheticWeather

This repository contains scripts to create synthetic weather time series from a short weather record of at least one year. The tool (scripts and how they might be implemented in ESP-r and Ladybug) is provisionally named `indra` <sup>(1)</sup>. _All the scripts here should be treated as **experimental** unless explicitly stated otherwise._

### Great, where do I begin?

Start by going into either the folder `m-files` (MATLAB files) or the folder `py-files` (Python files, if you are feeling adventurous). Each of them has a `README.md` or `Installation.md` file to get you started.

## The methods

The MATLAB/R scripts are based on the algorithms published in Parag's thesis. While these scripts are well documented (in two conference papers and the thesis), no further work is being undertaken on that algorithm (call it the 'old' algorithm). These scripts will, hopefully by November 2017, be faithfully translated to Python (update: translation about 60% complete).

The Python scripts in the repository *should not be used for now* (October 2017). If you are really curious, they can be called with two methods: `arma` and `gp`. The method `arma` will be an almost-completely-faithful translation of the MATLAB scripts<sup>(2)</sup>. The method `gp` is new and experimental - which means that it has not been published or tested extensively yet. This consitutes a significant change in the algorithm, and I would not recommend using it just yet (October 2017). 

## License, implementation, and compatibility

I intend to make this tool fully compatible with [ESP-r](https://github.com/ESP-rCommunity "ESP-r") and [Dragonfly](https://github.com/chriswmackey/Dragonfly "Dragonfly"). However, it exists on its own as a command-line tool.

Like the older work, the newer work will also, most likely, be released with a BSD-3 license. Using the older scripts requires only a valid MATLAB license and R (R is free to download and reuse). While you are free to use the scripts as you please, I am not liable for anything that happens as a result of using my scripts. Like if you accidentally release the nuclear missiles stored in the Dakotas, ruin the ski season in Switzerland, or cause a drought in Scotland.

This work is linked to Parag Rastogi's PhD thesis at the Ecole Polytechnique Federale de Lausanne, EPFL. You can find the full text of the thesis on [EPFL's scientific repository](https://infoscience.epfl.ch/record/220971?ln=en) or [my website](https://paragrastogi.com). 

## I'm panicking/clueless

If you have questions or concerns, or notice errors, please contact me at `contact[at]paragrastogi.com`.

Happy creating fake weather!


### Footnotes

(1) I haven't thought up a smarmy, contrived acronym yet but I'm working on it.

(2) The difference lies in the 'simulation' of the SARMA model to produce synthetic 'de-mean-ed' series. I am not convinced that I should reproduce the 'custom noise' functionality used in the old scripts to simulate the SARMA models with bootstrapped residuals. For now, I am doing a 'conventional' simulation by using white noise.

### Acknowledgements

1. This work began during Parag's PhD at the Ecole Polytechnique Federale de Lausanne. It was funded by the CCEM SECURE project and the EuroTech consortium. February 2012 - August 2016
2. Newer work and the Python translation is being undertaken by Parag as a visiting scientist at the Energy Systems Research Unit (ESRU), University of Strathclyde, Glasgow; and the RIKEN Institute for Advanced Intelligence Project (RIKEN-AIP), Tokyo. His stay at these institutions is being financed by the Swiss National Science Foundation (SNSF). 

Parag woud like to thank his hosts: Prof. Joe Clarke (Strathclyde, Glasgow) and Dr Mohammad Emtiyaz Khan and Prof. Masashi Sugiyama (RIKEN-AIP, Tokyo). The advice of Prof. Anthony Davison (EPFL, Lausanne) was crucial in creating the first models.

### Disclaimer

Parag is the sole author of these scripts and the thesis was supervised by Professor Marilyne Andersen. The scripts are the intellectual property of Parag and EPFL. Their redistribution under the liberal __ BSD-3 license __ is with the approval of EPFL, acting through Prof. Marilyne Andersen.
