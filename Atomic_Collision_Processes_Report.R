#'---
#'title: Atomic Physics Notes
#'author: ''
#'date: ''
#'output:
#'  bookdown::gitbook:
#'    css: "css/style-bookdown.css"
#'    config:
#'      split_by: chapter+section
#'      toc:
#'        collapse: section
#'        scroll_highlight: true
#'        before: null
#'        after: null
#'      toolbar:
#'        position: fixed
#'      edit : null
#'      search:
#'        engine: lunr # or fuse
#'        # options to control/tune search engine behavior (for
#'        # fuse.js, refer to https://fusejs.io/api/options.html)
#'        options: true
#'      fontsettings:
#'        theme: white
#'        family: sans
#'        size: 2
#'      sharing:
#'         facebook: false
#'         github: false
#'         twitter: false
#'         linkedin: false
#'         weibo: false
#'         instapaper: false
#'         vk: false
#'         whatsapp: false
#'         all: null
#'      info: true
#'---

#' # Atomic processes in plasma 
#' This Jupyter notebook contains information relating to atomic data, atomic physics, atomic modelling and coding of said processes for collisions inside plasmas.
#'
#' ## Charge exchange 
#' This section will outline what is involved in implementing charge exchange. An example code implmenting charge exchange can be found below. 
#'
#' ### Physics 
#' In Charge exchange an electron is transferred from one atom to another. An example is
#' $$ H + H^{+} \rightarrow H^{+} + H $$
#' If the 2 interacting atoms are one electron different that this means that the density of each atom is unchanged.   
#'  
#' ### Definitions 
#' The following subsections explains the meaning of various quantities used in the code 
#'
#' * Weight ($w$) - A macroparticle representing a particular number of physical particles 
#'
#' * Deltaweight ($dw$) - Change in the weight of a particle \@ref(eq:dw-charge-exchange) 
#'
#' * Ion density ($n_{ions}$) - number density of ions 
#'
#' * Time step ($d t$) - Time step interval 
#'
#' * Charge exchange rate ($R_{CE}$) - Rate at which charge exchange rate occurs 
#'
#' * Fluid picture decribes the ions of the plasma (H<sup>+</sup>) 
#'
#' * Particle picture describes the neutrals of the plasma (H) 
#' 
#' ### Equations 
#' \begin{equation}
#' dw = R_{CE} \ w \ n_{ions} \ d t (\#eq:dw-charge-exchange)
#' \end{equation}
#' The evolution of the velocity of the neutrals is
#' \begin{equation}
#' \dfrac{dv_n}{dt} = - \dfrac{dw}{w} v_n + \dfrac{dw}{w} v_i (\#eq:dvpdt-charge-exchange)
#' \end{equation}
#' The evolution of the velocity of the ions is
#' \begin{equation}
#' \dfrac{dv_i}{dt} = \dfrac{dw}{w} v_n + \dfrac{dw}{w} v_i (\#eq:dvidt-charge-exchange)
#' \end{equation}
#' The differential equations form a set out of coupled differential equations where the analytical solution is
#' \begin{equation}
#'  = (\#eq:analyical-solution-charge-exchange)
#' \end{equation}
#'
#' ### Results
#' The results of my momentum conservation are show in figure \@ref(fig:momentum-conservation) . This
#' came from running the python script [here](#code:charge-exchange-code)
#+ echo=FALSE,
py_run_file("./charge-exchange.py")
#+, momentum-conservation, out.width="75%", fig.cap="Momentum Conservation test results.", echo=FALSE
knitr::include_graphics("Momentum.png")
#' ## Ionise 
#' This section will outline what is involved in implementing ionisation 
#'
#' ### Physics  
#' Electron impact ionisation is where an electron and atom/ion collide resulting in the atom/ion losing an electron, e.g.
#' \begin{equation}
#' e + H \rightarrow 2e + H^{+}  
#' \end{equation}
#'
#' ### Definitions 
#' * Weight ($w$)- A macroparticle representing a particular number of physical particles 
#'
#' ## Electron Capture
#' This section will outline what is involved in implementing electron capture
#'
#' ### Physics 
#' Electron capture is where an electron and atom/ion collide resulting in the atom/ion gaining an electron, e.g.
#' \begin{equation}
#'  e + H^{+} \rightarrow H 
#' \end{equation}   
#'
#' # Atomic Databases and Atomic Data
#'
#' ## Atomic Databases 
#' This section will outline various databases of atomic data    
#'
#' * <a href="https://www.adas.ac.uk/about.php"> ADAS </a>   
#'
#' * <a href="https://open.adas.ac.uk">OpenADAS </a> 
#'
#' * <a href="https://www-amdis.iaea.org/ALADDIN/">ALADDIN </a>  
#'
#' ## Atomic Data 
#' This section will outline various atomic data, analytical formula which have been used in NESO.    
#'
#' *  <a href="https://scientific-publications.ukaea.uk/wp-content/uploads/CLM-R137.pdf">Charge exchange data</a>
#'
#' *  <a href="https://link.springer.com/article/10.1007/BF01392963">Electron impact ionisation rate data</a> 
#'
#' # Codes
#'
#' ## Charge exchange code
#' <a id="code:charge-exchange-code"></a>
#' Below is the python script which was used to test the charge exchange impementation.

#+ echo=FALSE
system2("pygmentize", args=c("-O","full,style=emacs,linenos=1","-o","test.html","./charge-exchange.py"))

#' <iframe src="../test.html" width=850, height=1500 style="overflow:hidden"></iframe>