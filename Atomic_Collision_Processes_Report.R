#'---
#'title: Atomic Physics Notes
#'author: ''
#'date: ''
#'output:
#'  bookdown::html_document2:
#'    css: "style.css"
#'    toc: true
#'    toc_float: true
#'    toc_collapsed: true
#'    toc_depth: 3
#'---

#' <style>
#' .row {
#'    position: absolute;
#'    left: 10px;
#' 
#' }
#' .tocify-subheader {
#'    text-indent: 40px;
#' }
#' .tocify-subheader .tocify-subheader {
#'    text-indent: 80px; 
#' }
#' .toc-content {
#'    position: absolute;
#'    left: 275px;
#'    min-width: 1000px; 
#'    width: 100%;
#' }
#'.tocify {
#'  max-width: 1600px;
#'  mix-width: 1600px;
#'}
#'.tocify-header {
#'  max-width: 1600px;
#'  mix-width: 1600px;
#' }
#' </style>

#' <h1 class="h1first"> Atomic processes in plasma </h1>
#' <p class="p1">This Jupyter notebook contains information relating to atomic data, atomic physics, atomic modelling and coding of said processes for collisions inside plasmas.</p>
#' <div id="divh2first"><h2 class="h2first">Charge exchange</h2></div>
#' <p class="p2">This section will outline what is involved in implementing charge exchange. An example code implmenting charge exchange can be found below. </p> 

#' <h3 class="h3first"> Physics </h3>
#' <p class="p3"> In Charge exchange an electron is transferred from one atom to another. An example is
#' $$ H + H^{+} \rightarrow H^{+} + H $$
#' If the 2 interacting atoms are one electron different that this means that the density of each atom is unchanged.   
#' </p>   
#' <h3 class="h3next"> Definitions </h3>
#' <p class="p3">The following subsections explains the meaning of various quantities used in the code </p>
#' <ul>
#'   <li> Weight ($w$) - A macroparticle representing a particular number of physical particles </li>
#'   <li> Deltaweight ($dw$) - Change in the weight of a particle \@ref(eq:dw-charge-exchange) </li>
#'   <li> Ion density ($n_{ions}$) - number density of ions </li>
#'   <li> Time step ($d t$) - Time step interval </li>
#'   <li> Charge exchange rate ($R_{CE}$) - Rate at which charge exchange rate occurs </li>
#'   <li> Fluid picture decribes the ions of the plasma (H<sup>+</sup>) </li>
#'   <li> Particle picture describes the neutrals of the plasma (H) </li>
#' </ul>
#' <h3 class="h3next"> Equations </h3>
#' \begin{equation}
#' dw = R_{CE} \ w \ n_{ions} \ d t (\#eq:dw-charge-exchange)
#' \end{equation}
#' <h3 class="h3next"> Results </h3>
#' <p class="p3">
#' The results of my momentum conservation are show in figure \@ref(fig:momentum-conservation) . This
#' came from running the python script [here](#code:charge-exchange-code)
#' </p>
#' <div id="divh3first">
#+ echo=FALSE,
py_run_file("./charge-exchange.py")
#+, momentum-conservation, out.width="75%", fig.cap="Momentum Conservation test results.", echo=FALSE
knitr::include_graphics("Momentum.png")
#' </div>
#' <div id="divh2next"><h2 class="h2next" > Ionise </h2></div>
#' <p class="p2"> This section will outline what is involved in implementing ionisation </p>
#' <h3 class="h3first"> Physics </h3>    
#' <p class="p3">
#' Electron impact ionisation is where an electron and atom/ion collide resulting in the atom/ion losing an electron, e.g.
#' \begin{equation}
#' e + H \rightarrow 2e + H^{+}  
#' \end{equation}
#' <h3 class="h3next"> Definitions </h3>
#' <ul>
#' <li> Weight ($w$)- A macroparticle representing a particular number of physical particles </li>
#' </ul>
#' <div id="divh2next"><h2 class="h2next" > Electron Capture </h2></div>
#' <h3 class="h3first"> Physics </h3>  
#' <p class="p3">
#' Electron capture is where an electron and atom/ion collide resulting in the atom/ion gaining an electron, e.g.
#' \begin{equation}
#'  e + H^{+} \rightarrow H 
#' \end{equation}    
#' </p>   
#'<h1 class="h1next">Atomic Databases and Atomic Data</h1>
#'<div id="divh2first"><h2 class="h2first">Atomic Databases</h2></div>
#'<p class="p2">
#'This section will outline various databases of atomic data    
#'<ul>
#'<li> <a href="https://www.adas.ac.uk/about.php"> ADAS</a> </li>   
#'<li> <a href="https://open.adas.ac.uk">OpenADAS </a> </li>
#'<li> <a href="https://www-amdis.iaea.org/ALADDIN/">ALADDIN</a></li>   
#'</ul>
#'</p>    
#'<div id="divh2first"><h2 class="h2next" >Atomic Data</h2></div>
#'<p class="p2">
#'This section will outline various atomic data, analytical formula which have been used in NESO.    
#'<ul>
#'<li> <a href="https://scientific-publications.ukaea.uk/wp-content/uploads/CLM-R137.pdf">Charge exchange data</a> </li>
#'<li> <a href="https://link.springer.com/article/10.1007/BF01392963">Electron impact ionisation rate data</a> </li>
#'</ul>    
#'</p>    
#'<h1 class="h1next">Codes</h1>
#'<div id="divh2first"><h2 class="h2first">Charge exchange code</h2></div>
#'<a id="code:charge-exchange-code"></a>
#' <p class="p2">Below is the python script which was used to test the charge exchange impementation. </p> 
#+ , echo=FALSE
system2("pygmentize", args=c("-O","full,style=emacs,linenos=1","-o","test.html","charge-exchange.py"))
#' <div id="divh2code">
#' <iframe src="./test.html" width=940, height=2150 style="overflow:hidden"></iframe>
#' </div>