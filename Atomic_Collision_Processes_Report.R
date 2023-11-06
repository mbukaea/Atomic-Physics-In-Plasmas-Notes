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
#' * $N_{CE}$ - The number of charge exchange events which happen in timestep dt \@ref(eq:dw-charge-exchange) 
#'
#' * Ion density ($n_{ions}$) - number density of ions 
#'
#' * Time step ($dt$) - Time step interval 
#'
#' * V - Volume of the cell 
#'
#' * $v_{p}$ - Mean velocity of a particle within a macroparticles
#'
#' * $v_{f}$ - Mean velocity of an ion in the fluid
#'
#' * Charge exchange rate ($R_{CE}$) - Rate at which charge exchange rate occurs 
#'
#' * Fluid picture decribes the ions of the plasma (H<sup>+</sup>) 
#'
#' * Particle picture describes the neutrals of the plasma (H) 
#' 
#' ### Equations 
#' 
#' #### Single Macroparticle exchanging with fluid
#'
#' In timestep $dt$ the number of charge exchange events which occur $N_{CE}$ is
#' \begin{equation}
#' N_{CE} = R_{CE} \ w \ n_{ions} dt (\#eq:dw-charge-exchange)
#' \end{equation}
#' The change in the average velocity of particle within a macroparticle in this timestep $dt$ is
#' \begin{equation}
#' dv_{p} = - \dfrac{N_{CE} v_{p}}{w} + \dfrac{N_{CE} v_{f}}{w} (\#eq:dvpdt-charge-exchange)
#' \end{equation}
#' This equation can be found by thinking about the average velocity of particle before and averge charge exchange 
#' has taken place. Assuming that the average velocity of a particle before a collision is $v_{p}$, the average velocity
#' of a particle after change exchange is
#' \begin{equation}
#' = \dfrac{(w-N_{CE})v_{p} + N_{CE} v_{f} }{w}
#' \end{equation}
#' which means that the average velocity of a particle changes by
#' \begin{equation}
#' = - \dfrac{N_{CE} v_{p}}{w} + \dfrac{N_{CE} v_{f}}{w}
#' \end{equation}
#' By similar arguments you find that the change in the average velocity of the fluid is
#' \begin{equation}
#' dv_{f} = \dfrac{N_{CE} v_{p}}{N_{ions}} - \dfrac{N_{CE} v_{f}}{N_{ions}} (\#eq:dvfdt-charge-exchange)
#' \end{equation}
#' where $N_{ions}$ is the number of ions ($=n_{ions}V$ where V is the volume of the cell). These couple differential equations form a set out of coupled differential equations, which when written in
#' matrix form, can be seen to have the analytical solution  (assumming w is time independent)
#' \begin{equation}
#' \begin{pmatrix} v_{p} (t) \\ v_{f} (t) \end{pmatrix} = e^{ R_{CE} \ n_{ions} \begin{pmatrix} -1 & 1 \\ \frac{w}{V} & - \frac{w}{V} \end{pmatrix} t }\begin{pmatrix} v_{p} (t=0) \\ v_{f} (t=0) \end{pmatrix} (\#eq:analyical-solution-charge-exchange)
#' \end{equation}
#' This result may at first look to be dependent on the volume of the cell your in, but as w is constructed using the formula
#' \begin{equation}
#' \dfrac{V n_{neutrals}}{n_{macroparticles}},
#' \end{equation}
#' It can easily be seen that $\frac{w}{V}$ is independent ont the volume.
#'
#' #### Total Energy
#'
#' The average kinetic energy of a ion is
#' \begin{equation}
#' \dfrac{m}{2 \sigma \sqrt{2 \pi} } \int_{-\infty}^{\infty} v^{2} e^{-\dfrac{1}{2} \Big( \dfrac{v - \mu}{\sigma} \Big)^{2} } dv
#' \end{equation}
#' where $\mu$ is the mean velocity, and $\sigma$ is the particle thermal velocity defined by
#' \begin{equation}
#' \sigma \ = \ \sqrt{\dfrac{kT}{m}}
#' \end{equation}
#' This integral can be solved by making the substitution 
#' \begin{equation}
#' u = \dfrac{x - \mu}{\sqrt{2} \sigma}
#' \end{equation}
#' After substitution you get
#' \begin{equation}
#' \dfrac{m}{2 \sqrt{\pi}} \int_{-\infty}^{\infty} ( 2 \sigma^{2} u^{2} + 2 \sqrt{2} \mu \sigma u + \mu^{2} ) e^{-u^2} du 
#' \end{equation}
#' The integral over $u e^{-u^2}$ can be seen to be 0 by symmetry. The other integrals are known and are as follows 
#' \begin{equation}
#' \int_{-\infty}^{\infty} x^{2} e^{-x^2} \ = \ \dfrac{\sqrt{\pi}}{2}
#' \end{equation}
#' and 
#' \begin{equation}
#' \int_{-\infty}^{\infty} e^{-x^2} \ = \ \sqrt{\pi}
#' \end{equation}
#' leading to the average kinetic energy of a particle to be
#' \begin{equation}
#' \ = \ \dfrac{m}{2} ( \mu^2 + \sigma^2 )
#' \end{equation}
#' This can then be multiplied by the number of ions to get the total energy. The total energy is the macroparticles can be
#' calculated by simply performing the sum
#' \begin{equation}
#' \dfrac{1}{2} m_{neutral} \sum_{n=0}^{n_{macroparticles}} w_{n} v_{n}^{2}
#' \end{equation}
#'
#'
#' #### Multiple particles added
#' 
#' The rate of change of the $n^{th}$ weight weight is 
#' \begin{equation}
#' N_{CE,n} = R_{CE,n} \ w_{n} \ n_{ions} (\#eq:dw-charge-exchange-n)
#' \end{equation}
#' The evolution of the $n^{th}$ particles velocity is
#' \begin{equation}
#' \dfrac{dv_{pn}}{dt} = - N_{CE,n} \dfrac{v_{pn}}{w_{n}} + N_{CE,n} \dfrac{v_{f}}{w_{n}} (\#eq:dvpdt-charge-exchange-n)
#' \end{equation}
#' The evolution of the velocity of the fluid is
#' \begin{equation}
#' \dfrac{dv_{f}}{dt} = \sum_{n} \Bigg( N_{CE,n} \dfrac{v_{pn}}{N_{ions}} - N_{CE,n} \dfrac{v_{f}}{N_{ions}} \Bigg) (\#eq:dvfdt-charge-exchange-n)
#' \end{equation}
#' If we now sum equation \@ref(\#eq:dvpdt-charge-exchange-n) over n we get (replacing the sum of the velocities of the particles with the $n_{macroparticles} * v_{all \ macroparticles}$ and dividing by $n_{max}$)
#' \begin{equation}
#' \dfrac{d v_{all \ macroparticles} }{dt} \ = \ - R_{CE,n} n_{ions} v_{all \ macroparticles} +  R_{CE,n} n_{ions} v_{f}
#' \end{equation}
#' where $v_{all \ macroparticles}$ is the average velocity of all the macroparticles. Writing out $w$ explicity it can also be seen that \@ref(\#eq:dvfdt-charge-exchange-n)
#' can be written as follows (againt replacing the sum of the velocities of the particles with the $n_{macroparticles} * v_{all \ macroparticles}$
#' \begin{equation}
#' \dfrac{dv_{f}}{dt} = R_{CE,n} n_{neutrals} v_{all \ macroparticles} -  R_{CE,n} n_{neutrals} v_{f}
#' \end{equation}
#' The general solution of these two coupled differential equations is
#' \begin{equation}
#' v_{all \ macroparticles} \ = \ \dfrac{c_1 (ae^{-(a+b)t} + b) - ac_2 ( e^{-(a+b)t} - 1) }{a+b}
#' \end{equation}
#' and
#' \begin{equation}
#' v_{f} \ = \ \dfrac{c_2 (be^{-(a+b)t} + a) - bc_2 ( e^{-(a+b)t} - 1) }{a+b}
#' \end{equation}
#' where due to the initial conditions 
#' \begin{equation}
#' c_{1} \ = \ v_{all \ macroparticles} (t=0) \quad \textrm{and} \quad c_2 \ = \ v_{f} (t=0)
#' \end{equation}
#'
#' ### Algorithm
#' 
#' #### Initial State
#' 
#'
#' #### Fluid Sampling
#' When randomly sampling from a normal distribution because we are sampling a finite number of samples the mean and standard deviation that your
#' returned sample has will not be the mean and standard deviation you requested. This can be rectified by making the following transform on your outputted sample
#' \begin{equation}
#' new_{sample} \ = \ \Bigg( old_{sample} - \mu_{old\_sample} \Bigg) \Bigg( \dfrac{\sigma_{new\_sample}}{\sigma_{old\_sample}} \Bigg) + \mu_{new\_sample}
#' \end{equation}
#'
#' #### Charge Exchange
#'
#'
#' #### Fluid Temperature Update
#'
#'
#' ### Results
#' One example of performing charge exchange between 2000 macroparticles and a fluid, and calculating the
#' respective momentums for a single ion and single neutral results in figure \@ref(fig:momentum-conservation-single). A mean momentum
#' for a macroparticle momentum is calculated, and the red, green, and blue lines represent 1,2 and 3 standard deviations from that value.
#' This figure make use the python script shown [here](#code:charge-exchange-code) 
#+, momentum-conservation-bulk, out.width="75%", fig.cap="Evolution of the bulk momentum of the particles and fluid.", echo=FALSE
knitr::include_graphics("./Momentum_Numerical_Bulk.png")

#+, momentum-conservation-single, out.width="75%", fig.cap="Evolution of the momentum of a single ion and neutral.", echo=FALSE
knitr::include_graphics("./Momentum_Numerical_Single.png")

#+, momentum-exp, out.width="75%", fig.cap="Check of exponential factor in evolution of single atom properties compared to expected solution.", echo=FALSE
knitr::include_graphics("./Momentum_Numerical_Slope.png")
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
