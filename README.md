LakeMichTribBacteria
====================

LakeMichTribBacteria is a research compendium that brings together the
data and analyses associated with *Sanfilippo et al. In Review.
Watershed-scale landuse is associated with seasonal and spatial
variation in Lake Michigan basin stream bacterial communities, Journal
of Great Lakes Research.*

Analysis scripts are provided in analysis/ and human-readable data in
extData/.

Note: All analyses originally conducted and package built using R
version 4.0.0 in Windows 10

### To get started using the package

    options(repos=structure(c(CRAN="http://cran.r-project.org")))
    install.packages("devtools")
    library(devtools)
    install_github("ScribnerLab/LakeMichTribBacteria")

### /analysis contents

1.  dbRDA.R: Distance-based redundancy analysis
2.  dbRDA\_correlationplot\_04102020.R: Variable correlation plot
3.  NMDSplots\_03042020.R: Non-metric multidimensional scaling
4.  heatmap.R: Code for generating community heatmap figure

#### /analysis/mothur

1.  16S\_V4\_Tiny\_data\_interactive\_mothur\_script\_PDS\_101019\_TM\_nodistseqscutoff\_10162019\_1112021.txt:
    Mothur script to process metabarcoding data. See
    <a href="https://mothur.org/" class="uri">https://mothur.org/</a>
    for more information
2.  Bacterial\_Mothur\_cons\_tax\_and\_shared\_file\_processing.R:
    Mothur post-processing script

### /extData contents

1.  stability.trim.contigs.good.unique.good.filter.unique.pick.precluster.pick.opti\_mcc.0.03.cons.taxonomy:
    The taxonomy file from Mothur (input for file 2. of
    /analysis/mothur)
2.  stability.trim.contigs.good.unique.good.filter.unique.pick.precluster.pick.opti\_mcc.0.03.subsample.shared:
    Subsample file (input for file 2. of /analysis/mothur)
3.  Bacterial\_OTU\_Matrix.txt
4.  2019Bacteria\_RiverCoordinates\_LkMich.csv: Site coordinates
5.  RiverNames\_01112021.xlsx: Site names
6.  bacteria\_state.xlsx: Key for site state based on naming scheme
7.  bacteria\_time.xlsx: Key for sampling time based on naming scheme
8.  bacteria\_rivers\_reduced\_landscape\_vars\_raw\_&\_transformed.xlsx:
    Landscape variables

#### Contact

Kim T. Scribner <br>
<a href="mailto:scribne3@msu.edu" class="email">scribne3@msu.edu</a>
<br>

Jared J. Homola <br>
<a href="mailto:jaredhomola20@gmail.com" class="email">jaredhomola20@gmail.com</a>
<br> www.jaredhomola.com

#### Copyright (c) 2019 ScribnerLab

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
“Software”), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
