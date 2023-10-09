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
#' If the 2 interacting atoms are one electron different that this means that the density of each atom is unchanged by the reaction.   
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
#' 
#' #### Single Particle added
#'
#' The rate of change of the weight with time is
#' \begin{equation}
#' \dfrac{dw}{dt} = R_{CE} \ w \ n_{ions} (\#eq:dw-charge-exchange)
#' \end{equation}
#' The evolution of the velocity of the particle is
#' \begin{equation}
#' \dfrac{dv_{p}}{dt} = - \dfrac{dw}{dt} \dfrac{v_{p}}{w} + \dfrac{dw}{dt} \dfrac{v_{f}}{w} (\#eq:dvpdt-charge-exchange)
#' \end{equation}
#' This evolution equation can be found by thinking about the average velocity of particle before and averge charge exchange 
#' has taken place. Assuming that the average velocity of a particle before a collision is $v_{p}$, the average velocity
#' of a particle after change exchange is
#' \begin{equation}
#' = \dfrac{(w-\dfrac{dw}{dt})v_{p} + \dfrac{dw}{dt} v_{f} }{w}
#' \end{equation}
#' which means that the average velocity of a particle changes by
#' \begin{equation}
#' = - \dfrac{dw}{dt} \dfrac{v_{p}}{w} + \dfrac{dw}{dt} \dfrac{v_{f}}{w}
#' \end{equation}
#' The evolution of the velocity of the fluid is
#' \begin{equation}
#' \dfrac{dv_{f}}{dt} = \dfrac{dw}{dt} \dfrac{v_{p}}{w} - \dfrac{dw}{dt} \dfrac{v_{f}}{w} (\#eq:dvfdt-charge-exchange)
#' \end{equation}
#' These differential equations form a set out of coupled differential equations, which when written in
#' matrix form, can be seen to have the analytical solution 
#' \begin{equation}
#' \begin{pmatrix} v_{p} (t) \\ v_{f} (t) \end{pmatrix} = e^{ R_{CE} \ n_{ions} \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix} t }\begin{pmatrix} v_{p} (t=0) \\ v_{f} (t=0) \end{pmatrix} (\#eq:analyical-solution-charge-exchange)
#' \end{equation}
#' 
#' #### Multiple particles added
#' 
#' The rate of change of the $n^{th}$ weight weight is 
#' \begin{equation}
#' \dfrac{dw_{n}}{dt} = R_{CE,n} \ w_{n} \ n_{ions} (\#eq:dw-charge-exchange-n)
#' \end{equation}
#' The evolution of the $n^{th}$ particles velocity is
#' \begin{equation}
#' \dfrac{dv_{pn}}{dt} = - \dfrac{dw_{n}}{dt} \dfrac{v_{pn}}{w_{n}} + \dfrac{dw_{n}}{dt} \dfrac{v_{f}}{w_{n}} (\#eq:dvpdt-charge-exchange-n)
#' \end{equation}
#' The evolution of the velocity of the fluid is
#' \begin{equation}
#' \dfrac{dv_{f}}{dt} = \sum_{n} \Bigg( \dfrac{dw_{n}}{dt} \dfrac{v_{pn}}{w_{n}} - \dfrac{dw_{n}}{dt} \dfrac{v_{f}}{w_{n}} \Bigg) (\#eq:dvfdt-charge-exchange-n)
#' \end{equation}
#' For 2 macroparticles and a fluid this equations to (writing out in matrix form)
#' \begin{equation}
#' \dfrac{d}{dt} \begin{pmatrix} v_{p1} \\ v_{p2} \\ v_{f} \end{pmatrix}  =
#' \begin{pmatrix} - \dfrac{dw_{1}}{dt} \dfrac{1}{w_{1}} & 0 &  \dfrac{dw_{1}}{dt} \dfrac{1}{w_{1}} \\
#'                 0 &  - \dfrac{dw_{2}}{dt} \dfrac{1}{w_{2}} & \dfrac{dw_{2}}{dt} \dfrac{1}{w_{2}} \\
#'                 \dfrac{dw_{1}}{dt} \dfrac{1}{w_{1}} & \dfrac{dw_{2}}{dt} \dfrac{1}{w_{2}} & - \dfrac{dw_{1}}{dt} \dfrac{1}{w_{1}} - \dfrac{dw_{2}}{dt} \dfrac{1}{w_{2}} \end{pmatrix}
#' \begin{pmatrix} v_{p1} \\ v_{p2} \\ v_{f} \end{pmatrix}
#' \end{equation}
#'
#' ### Results
#' Performing charge exchange between just 100 macroparticles and a fluid, and calculating their
#' respective momentums results in figure \@ref(fig:momentum-conservation). The figure also shows a comparison
#' between the numerical solution which includes sampling the fluid momentum at each time step, and the 
#' analytical solution which doesn't (see \@ref(eq:analyical-solution-charge-exchange)). This figure make use the python script shown [here](#code:charge-exchange-code) 
#+ echo=FALSE,
system2("python", args=c("./unittests.py","2>","/dev/null","&>","/dev/null"))
#+, momentum-conservation, out.width="75%", fig.cap="Momentum Conservation test results.", echo=FALSE
knitr::include_graphics("MomentumNumerical.png")
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
#' 
#' ### Atomic data/formulas used within NESO
#' 
#' This section will outline various atomic data, analytical formula which have been used in NESO.    
#'
#' *  <a href="https://scientific-publications.ukaea.uk/wp-content/uploads/CLM-R137.pdf">Charge exchange data</a>
#'
#' *  <a href="https://link.springer.com/article/10.1007/BF01392963">Electron impact ionisation rate data</a> 
#'
#' ### Working with the data/formulas
#'
#' Each atomic database serves their data in different file formats.
#'
#' # Codes
#'
#' ## Charge exchange code
#' <a id="code:charge-exchange-code"></a>
#' Below is the python script which was used to test the charge exchange impementation.

#+ echo=FALSE
system2("pygmentize", args=c("-O","full,style=emacs,linenos=1","-o","test.html","./chargeexchange.py"))

#' <iframe src="../test.html" width=850, height=1500 style="overflow:hidden"></iframe>