
A Jupyter Notebook based on a R kernel illustrating how the climate projections are constrained by observations using the Kriging for Climate Change (KCC) package. 

Content
-------

* [Requirements](#requirements)
* [Related Content](#related-content)
* [Contributors](#contributors)
* [License](#license)


Requirements
------------

- Clone the repository:

<code>git clone git@gitlab.com:saidqasmi/KCC_notebook.git</code>

- Install the KCC R package 

Follow the instructions available at: <a href="https://gitlab.com/saidqasmi/KCC">https://gitlab.com/saidqasmi/KCC</a> 

- If necessary, install the Jupyter Notebook interface

Follow the official documentation to install the interface: <a href="https://jupyter.readthedocs.io/en/latest/install/notebook-classic.html">https://jupyter.readthedocs.io/en/latest/install/notebook-classic.html</a>

Install the R kernel for Jupyter Notebook by running the following command in a R console:

<code>install.packages("IRkernel")</code>

Make the kernel available to Jupyter by running in a R console (if you have root access, set `user = FALSE`):
<code>IRkernel::installspec(user = TRUE)</code>


Related Content
---------------

Qasmi, S. (in revision): Past and future response of the North Atlantic warming hole to anthropogenic forcings _Environment research letters_

Ribes, A. et al. (2021): Making climate projections conditional on historical observations. _Science Advances_, 7, eabc0671, DOI: <a href="https://doi.org/10.1126/sciadv.abc0671">10.1126/sciadv.abc0671</a>

Qasmi, S. and Ribes, A. (in revision): Reducing uncertainty in local climate projections _Science Advances_ DOI:  <a href="https://doi.org/10.21203/rs.3.rs-364943/v1">10.21203/rs.3.rs-364943/v1</a>



Contributors
------------

- Saïd Qasmi (said.qasmi@meteo.fr)

License
-------

KCC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

