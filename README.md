<!-- utf-8 -->

<h1>Indra and Vali</h1>

Indra is a cross-platform, Python-based program to create synthetic weather time series from a weather record of at least one year. Vali is the child of Indra created for use with ESP-r on 31 May 2018 (<i>v2.78_esru</i>).

See the <a href='https://github.com/paragrastogi/SyntheticWeather/wiki'><i>wiki</i></a> if you want to find out more about <b>indra</b>. It contains a step-by-step guide to installing and running mighty <b>indra</b>. If you know your way around, go directly to the sample commands below.
<br>
<strong>Indra is NOT a weather forecasting tool.</strong> It is designed to be used to create variations on weather patterns learned from a source file.
<br>
You can call <b>indra</b> using the python and shell scripts called <kbd>vali.py</kbd> and <kbd>vali.sh</kbd> respectively (they are the same script, just written in the two different languages).

<h2>Bibliography</h2>
The program is based on the algorithms published in Rastogi (2016), Rastogi and Andersen (2015, 2016).
<br>
<ul>
	<li> Full wiki for Indra: https://github.com/paragrastogi/SyntheticWeather/wiki
	<li>Rastogi, Parag. 2016. On the Sensitivity of Buildings to Climate: The Interaction of Weather and Building Envelopes in Determining Future Building Energy Consumption. PhD, Lausanne, Switzerland: Ecole polytechnique federale de Lausanne. EPFL Infoscience. https://infoscience.epfl.ch/record/220971?ln=en.
	<li>Rastogi, Parag, and Marilyne Andersen. 2015. Embedding Stochasticity in Building Simulation Through Synthetic Weather Files. In Proceedings of BS 2015. Hyderabad, India. http://infoscience.epfl.ch/record/208743.
	<li>Rastogi, Parag, and Marilyne Andersen. 2016. Incorporating Climate Change Predictions in the Analysis of Weather-Based Uncertainty. In Proceedings of SimBuild 2016. Salt Lake City, UT, USA. http://infoscience.epfl.ch/record/208743.
</ul>
@phdthesis{rastogi_sensitivity_2016,
	address = {Lausanne, Switzerland},
	type = {{PhD}},
	title = {On the sensitivity of buildings to climate: the interaction of weather and building envelopes in determining future building energy consumption},
	shorttitle = {Sensitivity of {Buildings} to {Climate}},
	url = {https://infoscience.epfl.ch/record/220971?ln=en},
	language = {EN},
	school = {Ecole polytechnique federale de Lausanne},
	author = {Rastogi, Parag},
	month = aug,
	year = {2016},
	note = {doi:10.5075/epfl-thesis-6881}
}
<br>

<h2>License, implementation, and compatibility</h2>

This tool is distributed under the GNU General Public License v3 (GPLv3). Please read what this means <a href='https://en.wikipedia.org/wiki/GNU_General_Public_License'>here</a>.

These scripts come with absolutely no warranties/guarantees of any kind. Happy creating fake weather!
