<h1>Vali</h1>

Vali is a cross-platform, Python-based program to create synthetic weather time series from a weather record of at least one year.

<h2>Great, where do I begin?</h2>

<strong>Please, please see the <a href='https://github.com/paragrastogi/SyntheticWeather/wiki'>wiki</a> first.</strong> I can't make you do it though, so yes, be your own boss. The wiki contains a step-by-step guide to installing and running mighty <b>indra</b>.

If you know your way around Python, go directly to the sample commands below. This version of indra (<i>v2.78_esru</i>) has been

<h2>Bibliography</h2>
The program is based on the algorithms published in Rastogi (2016), Rastogi and Andersen (2015, 2016).

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
	note = {doi:10.5075/epfl-thesis-6881}
}

<ul>
<li>Rastogi, Parag. 2016. “On the Sensitivity of Buildings to Climate: The Interaction of Weather and Building Envelopes in Determining Future Building Energy Consumption.” PhD, Lausanne, Switzerland: Ecole polytechnique fédérale de Lausanne. EPFL Infoscience. https://infoscience.epfl.ch/record/220971?ln=en.
<li>Rastogi, Parag, and Marilyne Andersen. 2015. “Embedding Stochasticity in Building Simulation Through Synthetic Weather Files.” In Proceedings of BS 2015. Hyderabad, India. http://infoscience.epfl.ch/record/208743.
<li>———. 2016. “Incorporating Climate Change Predictions in the Analysis of Weather-Based Uncertainty.” In Proceedings of SimBuild 2016. Salt Lake City, UT, USA. http://infoscience.epfl.ch/record/208743.
</ul>


<h2>License, implementation, and compatibility</h2>

This tool is distributed under the GPLv3 license. Please read what this means <a href='https://en.wikipedia.org/wiki/GNU_General_Public_License'>here</a>.

<h2>I'm panicking/clueless</h2>

If you have questions or concerns, or notice errors, please contact me at `contact[at]paragrastogi.com`.

Happy creating fake weather!


<h2>Footnotes</h2>

(1) I haven't thought up a smarmy, contrived acronym yet but I'm working on it.

(2) The difference lies in the 'simulation' of the SARMA model to produce synthetic 'de-mean-ed' series. I am not convinced that I should reproduce the 'custom noise' functionality used in the old scripts to simulate the SARMA models with bootstrapped residuals. For now, I am doing a 'conventional' simulation by using white noise.

### Acknowledgements

1. This work began during my PhD at the Ecole Polytechnique Federale de Lausanne. It was funded by the CCEM SECURE project and the EuroTech consortium. February 2012 - August 2016
2. Newer work and the Python translation is being undertaken by me as a visiting scientist at the Energy Systems Research Unit (ESRU), University of Strathclyde, Glasgow; and the RIKEN Institute for Advanced Intelligence Project (RIKEN-AIP), Tokyo. His stay at these institutions is being financed by the Swiss National Science Foundation (SNSF).

I woud like to thank his hosts: Prof. Joe Clarke (Strathclyde, Glasgow) and Dr Mohammad Emtiyaz Khan and Prof. Masashi Sugiyama (RIKEN-AIP, Tokyo). The advice of Prof. Anthony Davison (EPFL, Lausanne) was crucial in creating the first models.

### Disclaimer

Parag is the sole author of these scripts. The original PhD thesis was supervised by Professor Marilyne Andersen. The scripts are the intellectual property of Parag and EPFL. Their redistribution under the liberal __GPLv3__ license is with the approval of EPFL, acting through Prof. Marilyne Andersen.
