<h1>SyntheticWeather</h1>

This repository contains scripts to create synthetic weather time series from a short weather record of at least one year. The tool (scripts and how they might be implemented in ESP-r and Ladybug) is provisionally named `indra` <sup>(1)</sup>. _All the scripts here should be treated as **experimental** unless explicitly stated otherwise._

### Great, where do I begin?

<strong>Please, please see the <a href='https://github.com/paragrastogi/SyntheticWeather/wiki'>wiki</a> first.</strong> I can't make you do it though, so yes, be your own boss. The wiki contains a step-by-step guide to installing and running mighty <b>indra</b>.

If you know your way around MATLAB or Python, go directly into either the folder `m-files` (MATLAB files) or the folder `py-files` (Python files). Most of the scripts explain themselves. **Sample Python commands** are <a href='https://github.com/paragrastogi/SyntheticWeather/wiki/Sample-Commands'>given here</a>.

## The methods

The MATLAB/R scripts are based on the algorithms published in Parag's thesis. While these scripts are well documented (in two conference papers and the thesis), I won't be working on these any more.

The Python scripts in the repository are translations of these original scripts. The method is an almost-completely-faithful translation of the MATLAB scripts<sup>(2)</sup>.

## License, implementation, and compatibility

This tool is distributed under the GPLv3 license. Please read what this means <a href='https://en.wikipedia.org/wiki/GNU_General_Public_License'>here</a>.

This tool is in the process of being made compatible with <a href='https://github.com/ESP-rCommunity'>ESP-r</a> and <a href='https://github.com/chriswmackey/Dragonfly'>Dragonfly</a>. However, it exists on its own as a command-line tool.

Using the older scripts requires only a valid MATLAB license and R (R is free to download and reuse). While you are free to use the scripts as you please, I am not liable for anything that happens as a result of using my scripts. Like if you accidentally release the nuclear missiles stored in the Dakotas, ruin the ski season in Switzerland, or cause a drought in Scotland.

This work is linked to my PhD thesis at the Ecole Polytechnique Federale de Lausanne, EPFL. You can find the full text of the thesis on [EPFL's scientific repository](https://infoscience.epfl.ch/record/220971?ln=en) or [my website](https://paragrastogi.com). Please cite it as:

Rastogi, Parag. 2016. ‘On the Sensitivity of Buildings to Climate: The Interaction of Weather and Building Envelopes in Determining Future Building Energy Consumption’. PhD, Lausanne, Switzerland: Ecole polytechnique fédérale de Lausanne. EPFL Infoscience. https://infoscience.epfl.ch/record/220971?ln=en.

@phdthesis{rastogi_sensitivity_2016,
	address = {Lausanne, Switzerland},
	type = {{PhD}},
	title = {On the sensitivity of buildings to climate: the interaction of weather and building envelopes in determining future building energy consumption},
	shorttitle = {Sensitivity of {Buildings} to {Climate}},
	url = {https://infoscience.epfl.ch/record/220971?ln=en},
	language = {EN},
	school = {Ecole polytechnique fédérale de Lausanne},
	author = {Rastogi, Parag},
	month = aug,
	year = {2016},
	note = {doi:10.5075/epfl-thesis-6881},
	file = {EPFL_TH6881.pdf:C\:\\Users\\prastogi\\AppData\\Roaming\\Zotero\\Zotero\\Profiles\\wjfrmd14.default\\zotero\\storage\\3295R3P6\\EPFL_TH6881.pdf:application/pdf}
}

## I'm panicking/clueless

If you have questions or concerns, or notice errors, please contact me at `contact[at]paragrastogi.com`.

Happy creating fake weather!


### Footnotes

(1) I haven't thought up a smarmy, contrived acronym yet but I'm working on it.

(2) The difference lies in the 'simulation' of the SARMA model to produce synthetic 'de-mean-ed' series. I am not convinced that I should reproduce the 'custom noise' functionality used in the old scripts to simulate the SARMA models with bootstrapped residuals. For now, I am doing a 'conventional' simulation by using white noise.

### Acknowledgements

1. This work began during my PhD at the Ecole Polytechnique Federale de Lausanne. It was funded by the CCEM SECURE project and the EuroTech consortium. February 2012 - August 2016
2. Newer work and the Python translation is being undertaken by me as a visiting scientist at the Energy Systems Research Unit (ESRU), University of Strathclyde, Glasgow; and the RIKEN Institute for Advanced Intelligence Project (RIKEN-AIP), Tokyo. His stay at these institutions is being financed by the Swiss National Science Foundation (SNSF).

I woud like to thank his hosts: Prof. Joe Clarke (Strathclyde, Glasgow) and Dr Mohammad Emtiyaz Khan and Prof. Masashi Sugiyama (RIKEN-AIP, Tokyo). The advice of Prof. Anthony Davison (EPFL, Lausanne) was crucial in creating the first models.

### Disclaimer

Parag is the sole author of these scripts. The original PhD thesis was supervised by Professor Marilyne Andersen. The scripts are the intellectual property of Parag and EPFL. Their redistribution under the liberal __GPLv3__ license is with the approval of EPFL, acting through Prof. Marilyne Andersen.
