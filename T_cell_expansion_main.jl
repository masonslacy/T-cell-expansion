## set-up (download packages)
using Pkg 
Pkg.add("Plots")
Pkg.add("LaTeXStrings")
Pkg.add("LaTeXStrings")
Pkg.add("Images")
Pkg.add("ImageView")
Pkg.add("Colors")
Pkg.add("SparseArrays")
Pkg.add("LinearAlgebra")
Pkg.add("ProgressMeter")
Pkg.add("JLD")
Pkg.add("Statistics")
Pkg.add("DifferentialEquations")
Pkg.add("Random")
Pkg.add("Distributions")
Pkg.add("StatsPlots")
Pkg.add("LsqFit")


## main script
using Plots
using Plots.PlotMeasures
using LaTeXStrings
using Images
using ImageView
using Colors
using SparseArrays
using LinearAlgebra
using ProgressMeter
using JLD
using Statistics
using DifferentialEquations
using Random
using Distributions
using StatsPlots
using LsqFit

include("T_cell_expansion_functions.jl")


# default parameters 
def_font = 13; # default font size
default(titlefont = (def_font, "times"), legendfont = (def_font, "times"), guidefont = (def_font, "times"), tickfont = (def_font, "times"), framestyle = :box, yminorgrid = true, xminorgrid = true, size = (600,400), linewidth = 2); # default plot settings

# plot colours
col_total = RGB(0, 0, 0); # black for total T cells
col_naive = RGB(44/255, 143/255, 163/255); # light blue for naive cells
col_active = RGB(237/255, 62/255, 62/255); # red for activated cells
col_IL2 = RGB(48/255, 194/255, 95/255); # green for IL-2


# choose simulation options
cell_IC = "uniform"; # initial condition for cells: "pinpoint", "uniform"
start_with_active = false; # true to start with activated cells (half of N0)
cyt_IC = "zero"; # "zero", "pinpoint", "uniform"
approx = "central"; # "finite", "central"
BC = "periodic"; # boundary conditions for PDEs: "noflux" or "periodic"

# switch on/off versions of the models
use_ODE_for_IL2 = true; # true to use ODE model for IL-2
supplement_IL2 = false; # true to supplement IL-2 (amounts supp_IL2_mass at times supp_IL2_time)


# biological parameters
T_diam = 5; # diameter of T cell (μm) (small estimate)

MSR_len = 70.4; # mean MSR length (μm)   (D. K. Y. Zhang, A. S. Cheung, and D. J. Mooney, “Activation and expansion of human T cells using artificial antigen-presenting cell scaffolds,” Nat Protoc, vol. 15, no. 3, pp. 773–798, Mar. 2020, doi: 10.1038/s41596-019-0249-0.)
MSR_diam = 5.1; # mean MSR diameter (μm)                                           ^
MSR_SA = 572*1000000; # surface area of MSR per gram (μm^2/μg)                     ^
m_MSR = π*MSR_diam*MSR_len/MSR_SA; # mass of individual MSR (μg)


τ = 1; # simulation timestep for PDE model (min)
τ_ABM = 1; # simulation timestep for cytokine in ABM model (min)
ε = 5*τ; # T cell step length (in x and y directions) (μm) 
α = 0.1; # proportion of ε for β
β = α*ε; # scaling used to influence chances of movement through the scaffold (higher represents greater chance to respond to scaffolds and move onto the micro-rods), ∈ [0,1]


# generate scaffold
APC_ms_conc = 333; # concentration of inputted APC-ms (μg/mL) 
APC_ms_vol = 100e-3; # input volume of APC-ms (mL) 
m_input = APC_ms_conc*APC_ms_vol; # input mass of MSRs (μg) assuming majority of APC-ms mass is MSRs       m_MSR*result_scale*n for n micro-rods
max_d = Inf; # maximum number of micro-rods overlapped at any lattice site in the domain
N_cell_2D = 1; # number of cells stacked vertically (number of cell diameters) to approximate this 3D domain in 2D
scale = 0.001; # scales down the domain area by the proportion scale ∈ [0,1]
linemethod = "original"; # method for generating the lines in scaffolds: "original", "XiaolinWu"
# probability density functions (for x and y) that govern the x and y positions of the centre of micro-rods (functions of the domain width and height)   -   check the plot of a multivariable PDF with code at the bottom of the script
xpos_PDF(max) = Uniform(1,max); # options:  Uniform(1,max), Normal(max/2,max/5), SkewNormal(max/2,max/5,-max/10), Arcsine(1,max), ...
ypos_PDF(max) = Uniform(1,max); # options cont        Exponential(max/5), Laplace(max/2), SymTriangularDist(max/2,max/5), TriangularDist(1,max,1)
# probability density function (of θ, in rad ∈ [0, π)) that governs the rotation of micro-rods    -    check plot with:   plot(x->pdf(rot_PDF,x),xlab="θ (rad)",ylab="Probability density",legend=false)
rot_PDF = Uniform(0,π); # options:  Uniform(0,π), Normal(π/2,π/20) 
blurring = 1; # side length (pixels, h μm) of the uniform kernel used to convolute (blur) the scaffold image (must be odd)
seed = NaN; # seed for random number generation (same seed will generate the same scaffold, NaN to use default seed)       for uniform distribution: 24 has small hole in the middle, 9319 has a big hole and cluster of micro-rods,  45, 69 have holes
ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

#display(heatmap(0:h:hei, 0:h:wid, reverse(ρ,dims=1), xlims=[0,hei], ylims=[0,wid], colorbar=false, c=cgrad(:grayC, rev=true), clims=(0,ρ_max), aspect_ratio=:equal, size=(300,300));) # view scaffold

h = T_diam; # distance between nodes in grid (μm)
mLtoh2 = h^2*N_cell_2D*T_diam/1e12; # ratio of mL to h^2 (mL per lattice site), approximation from 3D to 2D

# IL-2 loading
#IL2_stock = 0.2e3; # concentration of IL-2 in typical stock used to load micro-rods (μg/mL) (Zhang, et al. 2020)
IL2_M = 40e-9; # molar concentration of loaded IL-2 (mol/mL or kM) (Cheung, et al. 2018)
IL2_mass = 15.5e12; # molar mass of IL-2 (nDa, or ng/mol) (Zhang, et al. 2020)
IL2_vol = 10e-3; # volume of solution containing IL-2 loaded into micro-rods (mL)
MSR_loading = 500; # mass of MSRs used for loading IL-2 (in production of APC-ms) (μg)
IL2_ratio = (IL2_M*IL2_mass*IL2_vol/2)/MSR_loading; # ratio of loaded IL-2 (50% is retained in rods) to micro-rod mass (ng/μg)
I_loaded = IL2_ratio*m_input; # initial mass of loaded cytokine (retained inside MSRs with total mass m_input) (ng) 
I0_mass = I_loaded/result_scale; # initial mass of free cytokine if non-zero IC (ng)

k = 0.00023065269705856172; # secretion rate of cytokine (min^-1)   (as per model fit, see Fig 3 code)
K_S = 1.0341; # scaling to make IL-2 secretion accurate to the expected secretion (reduce error from scaffold discretisation)
IL2_IUtong = 18; # ratio of IU (international units) to ng of IL-2   (https://www.stemcell.com/international-units-conversion-data-for-cytokines)

# get IL-2 diffusivity from Stokes-Einstein-Sutherland equation: similar thing done in "Agent-Based Modelling Reveals the Role of the Tumor Microenvironment ... Checkpoint Blockade..." (2023)
k_B = 1.380649e-23; # Boltzmann constant (m^2 kg s^-2 K^-1)
T_media = 20+273.15; # temperature of the culture media (K) - room temperature (Zhang, et al. 2020)
η_media = 0.958e-3; # dynamic viscosity of the culture media (Pa s) - specific to RPMI-1640 + 10% FBS culture media (used in Zhang, et al. 2020), from https://fluidic.com/molecular-weight-to-hydrodynamic-radius-converter/
r_IL2 = 2.17e-9; # Stokes radius of IL-2 (m) - converted from a molecular mass of IL2_mass to hydrodynamic radius (https://fluidic.com/molecular-weight-to-hydrodynamic-radius-converter/)
D_I = k_B*T_media/(6*pi*η_media*r_IL2)*6e13; # diffusivity of IL-2 (μm^2/min) (Stokes-Einstein-Sutherland equation)
# we see a value around 6000  which is validated by "Competition for IL-2 between regulatory and effector T cells to chisel immune responses" (2012)


supp_IL2_mass = 30*[2,2,5,4,4,4,4]/IL2_IUtong; # mass of supplemented IL-2 (ng) (Zhang, et al. 2020)
supp_IL2_time = [6,8,9,10,11,12,13]*24*60; # times at which IL-2 supplements are given (min) (Zhang, et al. 2020)



mutable struct cell # T cells
    x::Float32 # x and y positions
    y::Float32
    activated::Bool # activation state (true for activated, false for naive)
end
T_loaded = 10*result_scale; # total number of loaded T cells
N0 = Int(floor(T_loaded/result_scale)); # initial number of T cells in scaled domain
halfN0 = Int(floor(T_loaded/result_scale/2)); # half of initial number of cells as an integer


ind(x, y) = [Int(round(y/h))+1, Int(round(x/h))+1]; # function to convert (x,y) coordinates to [i,j] indices (for indexing scaffold density)



T = 0.1*24*60; # total simulation time (min)
ntotal = 10; # total number of simulations for averaging/computing variance

T = Int(floor(T)); # make T into an integer for simplicity
Ts = Int(floor(T/τ)); # discretised number of timesteps for PDE model
Ts_ABM = Int(floor(T/τ_ABM)); # discretised number of timesteps for cytokine in ABM


# all probabilities and rates
ρ_avg_w = 333/m_MSR*mLtoh2; # average micro-rod concentraiton in 333 μg/mL scaffold
#min_w = 0.2; avg_w = 0.2; # defines linear w function
#w(ρ) = (avg_w-min_w)*ρ/ρ_avg_w+min_w; # density-dependent probability of waiting
#w_dash(ρ) = (avg_w-min_w)/ρ_avg_w; # derivative (w.r.t. ρ) of the probability of waiting
min_w = 0.15; max_w = 0.2; # minimum (no micro-rods) and maximum (inf micro-rods) wait probability
w(ρ) = (max_w-min_w)*ρ/(ρ_avg_w + ρ) + min_w; # density-dependent probability of waiting (Michaelis-Menten equation)
w_dash(ρ) = (max_w-min_w)*ρ_avg_w/(ρ_avg_w + ρ)^2; # derivative (w.r.t. c) of r_p (i.e. r_p'(c))  (1/min^2)


ra_max = 2*0.005; # maximum activation rate (1/min)
K_m_a = 1; # Michaelis constant (number of micro-rods at which the activation rate is half of rp_max)  (rods)
r_a(ρ) = ra_max*ρ*(MSR_len*MSR_diam/h^2)/(K_m_a + ρ*(MSR_len*MSR_diam/h^2)); # ρ-dependent activation rate (with ρ scaled to represent 1 if there is 1 rod at the lattice site) (Michaelis-Menten equation)  (1/min)
ra_dash(ρ) = ra_max*K_m_a*(MSR_len*MSR_diam/h^2)/(K_m_a + ρ*(MSR_len*MSR_diam/h^2))^2; # derivative (w.r.t. ρ) of r_a (i.e. r_a'(ρ))  (1/min^2)

rp_max = 1.5*0.0003188952794471368; # maximum proliferation rate (as per model fit, see Fig 3 code)  (1/min)
K_m_p = 0.6932514606417738*mLtoh2; # Michaelis constant - concentration of IL-2 at which the reaction rate is half of rp_max  (ng/h^2)     as per model fit (see Fig 3 code)
r_p(I) = rp_max*I/(K_m_p + I); # cytokine concentration (c, ng/h^2) dependent proliferation rate (Michaelis-Menten equation)  (1/min)
rp_dash(I) = rp_max*K_m_p/(K_m_p + I)^2; # derivative (w.r.t. c) of r_p (i.e. r_p'(c))  (1/min^2)

r_d = 0.5*0.2/(24*60); # constant death rate  (1/min)   -   0.4/day from "Mathematical Models of Tumor Immune System Dynamics" (pg 38)

#K_DI = 20; # scaling for IL-2 diffusion from T cell diffusion in the absence of micro-rods (i.e. D_n(0) or D_a(0))
K_λ = 2; # scaling for IL-2 decay from activated T cell concentration (1/Tcell/min),   ~2 for proportional to a*c,  ~0.00002 for proportional to R_na(c)a

m_n(ρ) = 1 - w(ρ) - τ*(r_d - r_a(ρ)); # probability of naive cell moving under the influence of micro-rod density ρ
m_a(ρ,I) = 1 - w(ρ) - τ*(r_d - r_p(I)); # probability of activated cell moving under the influence of micro-rod density ρ and IL-2 concentration I


s = true; # prevent multiple actions in one timestep (true) or allow it (false)


# consider PDE system:   ∂n/∂t = ∇•[D_n(ρ)∇n - χ_n(ρ)∇ρ*n] + R_nn(ρ)n + R_na(I)a
#                        ∂a/∂t = ∇•[D_a(ρ,I)∇a - (χ_aρ(ρ,I)∇ρ+χ_aI(ρ,I)∇I)a] + R_aa a + R_an(ρ)n
#                          ∂I/∂t = ∇•(D_I∇I) - λ(a,I) + S(ρ,t)

D_n(ρ) = 1/4*ε^2/τ*(1-w(ρ)); # scaffold-dependent diffusivity for n
χ_n(ρ) = ε^2/τ*(1/4*w_dash(ρ) + α/(2*ρ_bar)*(1-w(ρ))); # scaffold-dependent haptotactic sensitivity for n
R_nn(ρ) = -(r_a(ρ)+r_d); # scaffold-dependent reaction term driven by naive cells for n
R_na(I) = r_p(I); # cytokine-dependent reaction term driven by activated cells for n

D_a(ρ) = D_n(ρ); # scaffold- and cytokine-dependent diffusivity for a
χ_aρ(ρ) = χ_n(ρ); # density- and cytokine-dependent haptotactic sensitivity in the direction of ∇ρ for a
χ_aI(ρ) = 0; # scaffold- and cytokine-dependent haptotactic sensitivity in the direction of ∇c for a
R_aa(I) = -r_d; # cytokine-dependent (if s is true) reaction term driven by activated cells for a
R_an(ρ) = r_a(ρ); # scaffold-dependent reaction term driven by naive cells for a

#D_I = D_n(0)*K_DI; # constant diffusivity for IL-2 (h^2/min)  -  6000/h^2 from https://doi.org/10.3389/fimmu.2012.00268,   20x T cell diffusion from "Mathematical Models of Tumor Immune System Dynamics" (pg 38)
S(ρ,t) = ρ*k*m_MSR/m_input*I_loaded*exp(-k*t)#*K_S; # density- and time-dependent source term for cytokine (IL-2) (time is t=τ*n, where n is the discretised time)  (ng/h^2/min)
λ(a,I) = K_λ*a*I; # activated T cell- and cytokine-dependent or constant decay rate for IL-2  (1/min)    -    2/(60*24)*I from "Mathematical Models of Tumor Immune System Dynamics" (pg 38),  or K_λ*a*I (1/min?) to decay IL-2 proportionally to local activated T cells,  or K_λ*R_na(I)*a (1/min?) to represent consumption by activated T cells when they proliferate (THIS IS NON-LINEAR),  or 0 based off fitted data from (Cheung et al. 2018)



# parameters for the Theta method
atol = 1e-12; # FOM or GMRES absolute error tolerence
rtol = 1e-12; # FOM or GMRES relative error tolerence
maxiters_θ = 500; # FOM or GMRES maximum number of iterations

θ = 1/2; # theta selection for Theta method (forward-Euler: 0, backward-Euler: 1, Crank-Nicolson: 1/2)

lsolve = "backslash"; # method to solve linear system: "backslash", "FOM" or "GMRES"
pretype = "right"; # FOM or GMRES preconditioner type "none", "left" or "right"
precond = "Gauss-Seidel"; # FOM or GMRES preconditioner option "Jacobi" or "Gauss-Seidel"



# T cell initial condition
x_inj = wid÷2; y_inj = hei÷2; # injection position if using pinpoint IC
inj_ind = ind(x_inj,y_inj); # injection position in indices
if cell_IC == "pinpoint" # pinpoint injection at a single point in the domain (usually the centre)
    IC_u = zeros(N,M); # setting initial condition
    if start_with_active
        IC_u[inj_ind[1],inj_ind[2]] = halfN0; 
    else
        IC_u[inj_ind[1],inj_ind[2]] = N0; 
    end
elseif cell_IC == "uniform" # cells begin (approximately) uniformly spread out
    if start_with_active
        u0 = halfN0/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4); # concentration of cells at each point (conc = mass/volume)
    else
        u0 = N0/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4);
    end
    IC_u = zeros(N,M) .+ u0; # assume perfectly uniformly spread out (may not match with ABM)
end


# cytokine initial condition
IC_I = zeros(N,M); # initialise at zero
if cyt_IC == "pinpoint" # all cytokine starts at the same point in space
    xj = M÷2; yi = N÷2; # injection at (yi,xj) (indices)
    IC_I[yi,xj] = I0_mass; # setting initial condition
elseif cyt_IC == "uniform" # cytokine particles begin (perfectly) uniformly spread out
    I0 = I0_mass/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4); # concentration of particles at each point (conc = mass/volume)
    IC_I = IC_I .+ I0; # setting initial condition
end 
# otherwise left at zero



# initialise cytokine density if using a PDE for IL-2
I = zeros(N*M, Ts+1); I[:,1] = IC_I[:]; # cytokine density as a 2-dimensional array (for use in PDE solving)

t_PDE = 0:τ:T; # time vector

spI = sparse(1.0*LinearAlgebra.I, N*M, N*M); # generate sparse identity matrix for solving ODE systems (from PDE discretisation)
b_I = zeros(N*M, 1); # constant terms 

A_I = genA_I_IMEX(N, M, h, D_I, BC); # generate A matrix assuming using IMEX method
Ã_I = spI-τ*θ*A_I; b̃_I = spI+τ*(1-θ)*A_I; # generate matrices used in Theta method

alg = RK4(); # algorithm for solving ODEs





## code to run/test each model
# run ABM
n_ABM, a_ABM, I_ABM, n_tot_sims, a_tot_sims, I_tot_sims, runtime_ABM = simulate_ABM(); # outputs from ABM: average positions of naive and activated cells, average total mass/spatial concentration of IL-2, total number of naive and activated cells for each simulation, total mass of IL-2 for each simulation, computation time
T_ABM = n_ABM + a_ABM; # average positions of total T cells from ABM

# run PDE model 
n_PDE, a_PDE, I_PDE, runtime_PDE = simulate_PDE(); # outputs from PDE model: spatial concentrations of naive and activated cells, total mass/spatial concentration of IL-2, computation time
T_PDE = n_PDE + a_PDE; # spatial concentrations of total T cells from PDE model

# run ODE model 
n_ODE, a_ODE, I_ODE, runtime_ODE = simulate_ODE(); # outputs from ODE model: total number of naive and activated T cells and total mass of IL-2 (scaled up for the whole dish), computation time
T_ODE = n_ODE + a_ODE; # mass of total T cells from ODE model




## Fig 3 - model fitting for r_p and secretion rate
# ====          ===          ===          =============         =============
# ====  ===============  =======  =============================  ============
# ====  ===============  =======  =============================  ============
# ====       ==========  =======  ===     ================      =============
# ====  ===============  =======  ======  =====================  ============
# ====  ===============  =======  ======  =====================  ============
# ====  ===========          ===          =============         =============

# plots 3a and b
# data digitized (using WebPlotDigitizer) from "Low interleukin-2 concentration favors generation of early memory T cells over effector phenotypes during chimeric antigen receptor T-cell expansion" (Kaartinen, et al. 2017)
scale_t = 24*60; # scale for time (convert to minutes)
scale_IL2 = 1/18; # scale for IL2 amount (convert to IU to ng) (https://www.stemcell.com/international-units-conversion-data-for-cytokines)

IL2 = [0, 5, 20, 100, 300]*scale_IL2; # IL-2 concentration (ng/mL)

data0_t = [6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 10, 10, 10, 10, 10, 13, 13, 13, 13, 15]*scale_t; # data for 0 IU/mL IL-2  (min)
data0_y = [1.49, 2.23, 2.39, 4.01, 3.02, 7.34, 1.52, 3.51, 2.2, 7.59, 8.11, 3.63, 17.82, 8.82, 17.52, 13.18, 17.82, 7.21, 6.2, 2.6]; #  (fold expansion)
data5_t = [6, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 13, 13, 13, 13, 13, 15, 15, 15, 17, 17]*scale_t; # data for 5 IU/mL IL-2
data5_y = [1.6, 1.11, 2.55, 2.55, 5.16, 5.34, 9.12, 1.4, 1.49, 2.78, 3.69, 12.75, 20.04, 12.75, 11.72, 3.02, 21.07, 6.52, 13.86, 19.05, 14.1, 7.84, 42.57, 5.25, 4.44, 11.15, 2.92, 2.64];
data20_t = [3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 17, 17, 20, 20, 20, 20, 20, 20]*scale_t; # data for 20 IU/mL IL-2
data20_y = [1.17, 1.29, 1.39, 1.46, 3.18, 10.08, 1.82, 3.39, 4.53, 5.13, 5.96, 7.66, 12.01, 24.21, 1.12, 2.62, 5.68, 7.48, 8.07, 9.86, 14.71, 17.53, 21.41, 37.14, 40.04, 15.9, 23.73, 30.48, 34.54, 3.82, 11.2, 40.14, 61.42, 101.33, 17.64, 21.55, 36.45, 40.29, 50.47, 68.14, 71.64, 96.74, 99.19, 329.78, 107.19, 13.1, 31.45, 44.64, 134.28, 70.93, 21.66, 21.39, 344.11, 38.51, 57.48, 168.63, 482.45, 24.64, 42.73, 128.52, 21.74, 555.72, 173.56];
data100_t = [3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 17, 20, 20, 20, 20, 20, 20]*scale_t; # data for 100 IU/mL IL-2
data100_y = [1.11, 1.22, 1.35, 2.39, 2.74, 1.4, 1.89, 2.47, 3.03, 3.35, 4.84, 5.92, 9.16, 21.91, 24.23, 4.23, 5.01, 47.4, 56.06, 1.4, 16.2, 20.83, 36.85, 8.56, 2.65, 25.91, 25.91, 70.9, 143.43, 42.15, 16.75, 86.71, 95.89, 109.67, 22.66, 61.99, 125.42, 134.12, 200.62, 229.43, 354.89, 379.52, 464.16, 548.94, 627.79, 221.86, 300.08, 441.38, 849.11, 596.98, 419.71, 1207.72, 354.89, 742.46, 767.8, 1837.03, 1502.05, 2284.73, 587.04, 878.08, 1073.9, 4545.26];
data300_t = [3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 20, 20, 20, 20, 20, 20, 20]*scale_t; # data for 300 IU/mL IL-2
data300_y = [2.63, 1.73, 1.66, 1.13, 1.13, 18.99, 9.29, 2.63, 1.09, 3.18, 15.39, 5.85, 4.36, 2.13, 40.49, 30.16, 24.44, 20.66, 15.39, 5.97, 5.16, 4.55, 11.23, 1.29, 2.86, 140, 102.13, 67.06, 37.22, 18.99, 16.05, 3.61, 86.31, 59.11, 1122.61, 597.39, 484.1, 236.83, 208.75, 123.4, 52.11, 33.51, 392.29, 331.55, 1273.57, 1032.04, 572.78, 484.1, 872.25, 376.13, 9588.16, 2022.72, 1506.89, 1273.57, 584.96, 64980.85, 6566.77, 4892.13, 4690.65, 948.79, 3212.55, 2022.72];

data = [[data0_t,data0_y], [data5_t,data5_y], [data20_t,data20_y], [data100_t,data100_y], [data300_t,data300_y]]; # all data collated

T = 20*scale_t; # final time (min)

@. model(t, p) = exp.(p[1]*t); # model to fit to data, with unknown p = [rp] 
p0 = [0.0]; # initial guess for p = [rp]

rp = zeros(length(IL2)); # doubling rate (probability of reproducing) of T cells at corresponding IL-2 concentration
for i in eachindex(rp)
    fit = curve_fit(model, data[i][1], data[i][2], p0); # model fit
    rp[i] = fit.param[1]; # fitted rp value
end

@. rpmodel(C, p) = p[1]*C/(p[2] + C); # model (Michaelis-Menten equation) to fit to concentration (C) vs rp, with unknown p = [rp_max, K_m] where rp_max is the limiting rate, K_m is the Michaelis constant (concentration of substrate at which the reaction rate is half of rp_max)
p0 = [0.00035, 0.002]; # initial guess for p = [rp_max, K_m]
fit = curve_fit(rpmodel, IL2, rp, p0); # model fit
rp_max = fit.param[1]; # fitted rp_max value (limiting rate)
K_m = fit.param[2]; # fitted K_m value (Michaelis constant)

# plot parameters
def_font = 13; # default font size
default(titlefont = (def_font, "times"), legendfont = (def_font, "times"), guidefont = (def_font, "times"), tickfont = (def_font, "times"), framestyle = :box, yminorgrid = true, xminorgrid = true, size = (600,400), linewidth = 2); # default plot settings
col_line = "#4292c6"; # line colour

plot_x = 550; plot_y = 300; # default size for plots (x and y sizes)

rp_scale = 1e-4; # factor this scale out of proliferation rates for plotting

fit = curve_fit(model, data[3][1], data[3][2], p0); # model fit
rp[3] = fit.param[1]; # fitted rp value

plot(data[3][1]/scale_t, data[3][2], yaxis=:log, seriestype=:scatter, legend=false, c=:black, size=(plot_x,plot_y))
plot!(xlabel="𝑡 (days)", ylabel="Fold expansion, N_T(t)/N_T(0)", title="1.11 ng/mL IL-2")
display(plot!((1:T)/scale_t, exp.(rp[3]*(1:T)), yaxis=:log, c=col_line))

plot(IL2, rp/rp_scale, seriestype=:scatter, legend=false, c=:black, size=(plot_x,plot_y))
C = range(0, IL2[end], 100); # IL2 concentrations for plotting
plot!(xlabel="IL-2 concentration (ng/mL)", ylabel="Proliferation rate, r_p")
display(plot!(C, rp_max*C./(K_m .+ C)/rp_scale, c=col_line))



# data from (Cheung et al. 2018), digitised with WebPlotDigitizer, used to get fit for IL-2 secretion rate in Fig 3c
L0 = 1000; # initial loaded amount of IL-2 (ng)
I0 = 0; # initial amount of cytokine free in the dish (ng)

data_t = [0, 1, 3, 5, 7, 9, 11, 13, 15]*24*60; # time (minutes)
data_I = [0, 235.08, 589.37, 946.89, 992.17, 1005.26, 1005.47, 1004.07, 1004.29]/L0; # cumulative IL-2 (proportion of total)
data_std = [0, 0, 648.92-data_I[3], 977.47-data_I[4], 0, 0, 0, 0, 0]/L0; # standard deviation of cumulative IL-2 (proportion of total), only have data for third and fourth data points
paperfit_x = [0, 0.2615046837014796, 0.47598205492942336, 0.7140444440215594, 0.9521068331136967, 1.2373392579342233, 1.5698130801105619, 1.9021441790320492, 2.3053016931735466, 2.70831648406019, 3.1819436052845633, 3.6789416594908566, 4.223038387798108, 4.861403825934716, 5.617551630137443, 6.373342626203039, 7.128990899013779, 7.813598671722235, 8.851910350766182, 9.84298063245431, 10.857564570379218, 11.966203133251184, 12.862576535345417, 14.159966602758365, 15.009027245869346]*24*60; # model fit in paper (Cheung, et al. 2018)
paperfit_y = [0, 65.84894170706522, 136.22934675566262, 206.61331988563109, 276.99729301559955, 347.38840230831084, 420.05594751590183, 488.1848932192181, 558.5938429187858, 624.4641931140786, 685.8066480492109, 740.3447718093028, 792.6207319799997, 842.6416647240446, 888.1418383706713, 922.2955132566118, 951.9105886382773, 972.4377607672799, 990.7491543647151, 1006.7841120472697, 1020.5533380590587, 1025.259637387784, 1029.933823984167, 1036.9379677159998, 1037.066418645366]/L0;

@. model(t, k) = 1-exp.(-k[1]*t); # model to fit to data, with unknown k (secretion rate)  (normalised loaded cytokine)
k0 = [0.0]; # initial guess for k
fit = curve_fit(model, data_t, data_I, k0); # model fit for k
k_guess = fit.param[1]; # fitted k value

# fitting both k and L0 to update k fit based on uncertainty of L0 guess in previous fit
data_I = [0, 235.08, 589.37, 946.89, 992.17, 1005.26, 1005.47, 1004.07, 1004.29];
data_std = [0, 0, 648.92-data_I[3], 977.47-data_I[4], 0, 0, 0, 0, 0];
paperfit_y = [0, 65.84894170706522, 136.22934675566262, 206.61331988563109, 276.99729301559955, 347.38840230831084, 420.05594751590183, 488.1848932192181, 558.5938429187858, 624.4641931140786, 685.8066480492109, 740.3447718093028, 792.6207319799997, 842.6416647240446, 888.1418383706713, 922.2955132566118, 951.9105886382773, 972.4377607672799, 990.7491543647151, 1006.7841120472697, 1020.5533380590587, 1025.259637387784, 1029.933823984167, 1036.9379677159998, 1037.066418645366];

@. model(t, p) = p[2]*(1-exp.(-p[1]*t)); # model to fit to data, with unknown p = [k, L0] 
p0 = [k_guess, 1000.0]; # initial guess for p = [k, L0]
fit = curve_fit(model, data_t, data_I, p0); # model fit for k
k = fit.param[1]; # fitted k value
L0 = fit.param[2]; # fitted L0 value

covar = estimate_covar(fit); # covariance of model fit
CI_k = confidence_interval(fit)[1]; # confidence interval for k
CI_L0 = confidence_interval(fit)[2]; # confidence interval for L0

t_half = 10*24*60; # IL-2 half life (min)
λ_c = log(2)/t_half; # IL-2 decay rate (min^-1)
model_λ(t) = k*L0/(λ_c-k)*(exp.(-k*t) .- exp.(-λ_c*t));

# plot parameters
def_font = 13; # default font size
default(titlefont = (def_font, "times"), legendfont = (def_font, "times"), guidefont = (def_font, "times"), tickfont = (def_font, "times"), framestyle = :box, yminorgrid = true, xminorgrid = true, size = (600,400), linewidth = 2); # default plot settings
col_line = "#4292c6"; # line colour

plot_x = 550; plot_y = 300; # default size for plots (x and y sizes)
plot(data_t/24/60, data_I, yerror=data_std, seriestype=:scatter, c=:black, legend=false, size=(plot_x,plot_y))
plot!(xlabel="𝑡 (days)", ylabel="Mass of IL-2 (ng)")
display(plot!(0:0.1:15, model((0:0.1:15)*24*60, [k,L0]), c=col_line))




## Fig 4 - stochastic realisation with single cell tracking
# ====          ===          ===          =================    ==============
# ====  ===============  =======  ========================  =  ==============
# ====  ===============  =======  =======================  ==  ==============
# ====       ==========  =======  ===     ==============  ===  ==============
# ====  ===============  =======  ======  =============  ====  ==============
# ====  ===============  =======  ======  =============          ============
# ====  ===========          ===          ===================  ==============

T = 28*24*60; # total simulation time (min)
Ts_ABM = T; # number of timesteps for IL-2 solving

N0_track = 10; # initial number of cells
num_gens = 1; # number of cell generations to plot
animate_traj = false; # true to animate
place_cells = true; # true to decide where cells are placed

cell_xy = [[140,310], [340,130]]; # (x,y) for each cell if their starting position is chosen

num_sims = 1; # number of simulations (to average the activation time)

mutable struct cell_track # T cells for tracking
    x::Float32 # x and y positions
    y::Float32
    activated::Bool # activation state (true for activated, false for naive)
    ID::Int32 # unique ID for tracking
end

seed_ABM = Int(floor(rand()*10000)); # random seed for ABM simulation 
seed_ρ = 9319; # seed for scaffold generation 
ρ_1, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

# simulate indiviudal realisation
n_ABMxV, a_ABMxV, NI_ABM, cells_track_1, act_time_1, act_prop_1, avg_runtime_1 = ABM_track_cells();

N_ABM = sum(n_ABMxV,dims=[1,2])[:]; # total number of cells - integral of n_ABM * lattice site volumes (output of ABM) over the whole domain
A_ABM = sum(a_ABMxV,dims=[1,2])[:];
NT_ABM = N_ABM + A_ABM;


# plots 
ρ_1 = ρ_1*m_MSR/h^2*10^8; # scale ρ for plotting (values now represent mass density of micro-rods μg/μm^2)
ρ_max = maximum(ρ_1); # maximum value of scaled-up ρ

spat_plot_x = 300; spat_plot_y = spat_plot_x; # default size for spatial plots (x and y sizes)
pop_plot_x = 550; pop_plot_y = spat_plot_y; # default size for population plots (x and y sizes)
traj_plot_x = 500; traj_plot_y = traj_plot_x; # size for trajectory plots (x and y sizes)

T_colbar_scaling = 1e-3; # factor this scale out of the colorbar for displaying above, for T cells
I_colbar_scaling = 1e-8; # factor this scale out of the colorbar for displaying above, for IL-2

t_plots = [0, 7, 14, 21]; # times at which plots are generated for each population (days)
t_plots = Int.(floor.(t_plots*24*60)).+1; # convert to min and ensure whole numbers, +1 to correspond to array indices

IL2_scale = 1e-3; # factor this scale out of IL-2 mass for plotting

# spatial plots
n_plot = reverse(n_ABMxV,dims=1); # reversed matrices for plotting
a_plot = reverse(a_ABMxV,dims=1); 
titles = ["𝑡 = 0 days", "𝑡 = 7 days", "𝑡 = 14 days", "𝑡 = 21 days"];
for n = eachindex(t_plots)
    t_n = t_plots[n];
    # generate information for spatial plot
    naive_xind = zeros(0); naive_yind = zeros(0); # x and y indices for naive and activated cells
    active_xind = zeros(0); active_yind = zeros(0);
    for i = 1:N # for each node
        for j = 1:M 
            if n_plot[i,j,t_n] != 0
                append!(naive_yind,i);
                append!(naive_xind,j);
            end
            if a_plot[i,j,t_n] != 0
                append!(active_yind,i);
                append!(active_xind,j);
            end
        end
    end
    
    plt_spatial = heatmap(0:h:hei, 0:h:wid, reverse(ρ_1*m_MSR/h^2*10^8,dims=1), xlims=[0,hei], ylims=[0,wid], size=(spat_plot_x,spat_plot_y), colorbar=false, c=cgrad(:grayC, rev=true), clims=(0,maximum(ρ_1*m_MSR/h^2*10^8)), aspect_ratio=:equal); # plot scaffold
    for tick = [100, 200, 300, 400]
        plot!(plt_spatial, [tick,tick], [0,8], c=:black, legend=false, linewidth=2)
        plot!(plt_spatial, [0,8], [tick,tick], c=:black, legend=false, linewidth=2)
    end
    for edge = [0, wid]
        plot!(plt_spatial, [edge,edge], [0,wid], c=:black, legend=false, linewidth=2)
        plot!(plt_spatial, [0,wid], [edge,edge], c=:black, legend=false, linewidth=2)
    end
    
    plot!(plt_spatial, naive_xind*h, naive_yind*h, seriestype=:scatter, mc=col_naive, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), label="")
    plot!(plt_spatial, active_xind*h, active_yind*h, seriestype=:scatter, mc=col_active, alpha=0.7, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), label="")

    plot!(plt_spatial, xlabel="𝑥 (μm)", ylabel="𝑦 (μm)", title=titles[n])
    display(plt_spatial)
end

# population plots 
plt_pop = plot(legend=false, size=(pop_plot_x,pop_plot_y));
plt_IL2 = plot(legend=false, size=(pop_plot_x,pop_plot_y));

plot!(plt_pop, (0:T)/60/24, N_ABM, c=col_naive)
plot!(plt_pop, (0:T)/60/24, A_ABM, c=col_active)
plot!(plt_pop, (0:T)/60/24, NT_ABM, c=col_total)

plot!(plt_IL2, (0:T)/60/24, NI_ABM/IL2_scale, c=col_IL2)

# display vertical lines to indicate time points in spatial plots 
vline_width = 1.5; 
for t_n in t_plots/(24*60) # for each time point (in days)
    vline!(plt_pop, [t_n], c=col_total, linestyle=:dash, linewidth=vline_width)
    vline!(plt_IL2, [t_n], c=col_total, linestyle=:dash, linewidth=vline_width)
end

plot!(plt_pop, xlabel="𝑡 (days)", ylabel="Number of T cells")
plot!(plt_IL2, xlabel="𝑡 (days)", ylabel="Mass of IL-2 (pg)")

display(plt_pop)
display(plt_IL2)



# move micro-rods into cluster 
xpos_PDF(max) = Normal(max/4,max/10); # options:  Uniform(1,max), Normal(max/2,max/5), SkewNormal(max/2,max/5,-max/10), Arcsine(1,max), ...
ypos_PDF(max) = Normal(3max/4,max/10); # options cont        Exponential(max/5), Laplace(max/2), SymTriangularDist(max/2,max/5), TriangularDist(1,max,1)

seed_ABM = Int(floor(rand()*10000)); # random seed for ABM simulation 
seed_ρ = 9319; # seed for scaffold generation 
ρ_2, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

cell_xy = [[250,250], [300,100]]; # (x,y) for each cell if their starting position is chosen

# simulate indiviudal realisation
n_ABMxV, a_ABMxV, NI_ABM, cells_track_2, act_time_2, act_prop_2, avg_runtime_2 = ABM_track_cells();


xpos_PDF(max) = Uniform(1,max); # reset micro-rod position PDFs to default
ypos_PDF(max) = Uniform(1,max); 

ρ_2 = ρ_2*m_MSR/h^2*10^8; # scale ρ for plotting (values now represent mass density of micro-rods μg/μm^2)


# trajectory plots at t=14 days
traj_plot_x = 400; traj_plot_y = traj_plot_x; # size for trajectory plots (x and y sizes)

line_step = 5; # plot the cell trajectory as a line which skips line_step timesteps

traj_width = 3; # linewidth of trajectory lines
traj_alpha = 0.5; # alpha value (opacity) of trajectory lines
marker_size = 8; # size of markers for starting naive cells and death
act_size = marker_size/1.5; # size of marker for marking activation location
marker_line = 3; # linewidth of marker line

cells_plot_1 = [cells_track_1[1],cells_track_1[2]]; # plot only the first two cells
cells_plot_2 = [cells_track_2[1],cells_track_2[2]];
for i in eachindex(cells_plot_1)
    cells_plot_1[i][4][2,:] = abs.(cells_plot_1[i][4][2,:].-hei); # flip y coords for plotting
    cells_plot_2[i][4][2,:] = abs.(cells_plot_2[i][4][2,:].-hei);
end

plt_ρ = [ρ_1, ρ_2];
scaf_colours = ["#1b9e77", "#7570b3"];
ρ_max_colbar = maximum(ρ_2)*0.85;
plt_results = [cells_plot_1, cells_plot_2];
for i = eachindex(plt_ρ)
    cells_plot = plt_results[i];

    ρ = plt_ρ[i];
    ρ = ρ; # scale ρ for plotting (values now represent mass density of micro-rods μg/μm^2)
    ρ_max = maximum(ρ); # maximum value of scaled-up ρ

    plt = heatmap(0:h:hei, 0:h:wid, reverse(ρ,dims=1), xlims=[0,hei], ylims=[0,wid], size=(traj_plot_x,traj_plot_y), colorbar=false, c=cgrad(:grayC, rev=true), clims=(0,ρ_max_colbar), aspect_ratio=:equal, xlabel="𝑥 (μm)", ylabel="𝑦 (μm)"); # plot scaffold
    for edge = [0, wid]
        plot!(plt, [edge,edge], [0,wid], c=scaf_colours[i], legend=false, linewidth=8)
        plot!(plt, [0,wid], [edge,edge], c=scaf_colours[i], legend=false, linewidth=8)
    end
    for tick = [100, 200, 300, 400]
        plot!(plt, [tick,tick], [0,8], c=:black, legend=false, linewidth=2)
        plot!(plt, [0,8], [tick,tick], c=:black, legend=false, linewidth=2)
    end

    # plot cell trajectory in blue (naive) or red (activated) with colour proportional to 1/generation
    for cell in reverse(cells_plot) # [time of birth, time of death (NaN if not), time of activation (NaN if not), (x,y) positions over time, generation]
        birth_i = cell[1]; # time index of birth
        if T+1 >= birth_i # if cell was born by this time
            death_i = cell[2]; # time index of death (NaN if it never died)
            act_i = cell[3]; # time index of activation (NaN if it never got activated)
            cell_traj = cell[4]; # current cell's trajectory
            cell_gen = cell[5]^2; # current cell's generation squared (for plotting)
            if isnan(death_i)
                end_i = T+1;
            else
                end_i = min(T+1,death_i);
            end
            traj_path = birth_i:line_step:end_i; 
            if end_i == death_i
                traj_path = [traj_path;death_i];
            end
            cell_act = false; # keep track of when cell becomes activated
            for traj_i = length(traj_path)-1:-1:1 # for each plot time from the current time/death time back to their birth time
                # don't draw line if they cross through the boundary
                if maximum(abs.([diff(cell_traj[1,traj_path[traj_i]:traj_path[traj_i+1]]);diff(cell_traj[2,traj_path[traj_i]:traj_path[traj_i+1]])])) <= ε
                    if traj_path[traj_i] >= act_i
                        # plot trajectory with a red line (activated) 
                        plot!(cell_traj[1,traj_path[traj_i]:traj_path[traj_i+1]], cell_traj[2,traj_path[traj_i]:traj_path[traj_i+1]], c=col_active, alpha=traj_alpha, linewidth=traj_width, label="")
                    else
                        # plot trajectory with a blue line (naive) 
                        if !cell_act && !isnan(act_i) # if this is the time range in which the cell becomes activated, plot the segment when they were activated and then when they were naive
                            plot!(cell_traj[1,[traj_path[traj_i],act_i]], cell_traj[2,[traj_path[traj_i],act_i]], c=col_active, alpha=traj_alpha, linewidth=traj_width, label="")
                            plot!(cell_traj[1,[act_i,traj_path[traj_i+1]]], cell_traj[2,[act_i,traj_path[traj_i+1]]], c=col_naive, alpha=traj_alpha, linewidth=traj_width, label="")
                            cell_act = true
                        else
                            plot!(cell_traj[1,traj_path[traj_i]:traj_path[traj_i+1]], cell_traj[2,traj_path[traj_i]:traj_path[traj_i+1]], c=col_naive, alpha=traj_alpha, linewidth=traj_width, label="")
                        end
                    end
                end
            end
            if T+1 >= death_i # if cell died by this time
                # plot the end with a black cross where it died
                plot!([cell_traj[1,death_i]], [cell_traj[2,death_i]], seriestype=:scatter, markershape=:x, c=:black, markersize=marker_size, markerstrokewidth=marker_line, alpha=1/cell_gen, label="") 
                #=if isnan(act_i)
                    # plot the end with a blue cross (naive) where it died
                    plot!([cell_traj[1,death_i]], [cell_traj[2,death_i]], seriestype=:scatter, markershape=:x, c=col_naive, alpha=1/cell_gen, label="") 
                else
                    # plot the end with a red cross (activated) where it died
                    plot!([cell_traj[1,death_i]], [cell_traj[2,death_i]], seriestype=:scatter, markershape=:x, c=col_active, alpha=1/cell_gen, label="") 
                end=#
            end
            # plot the start of the cell trajectory with a blue circle (naive cell)
            plot!([cell_traj[1,birth_i]], [cell_traj[2,birth_i]], seriestype=:scatter, c=col_naive, markersize=marker_size, markerstrokewidth=marker_line, alpha=1/cell_gen, label="") # plot the start of the cell trajectory with a blue circle (naive cell)
            # plot the location of activation if activated 
            if !isnan(act_i) && T+1 >= act_i
                plot!([cell_traj[1,act_i]], [cell_traj[2,act_i]], seriestype=:scatter, c=col_active, markersize=act_size, markerstrokewidth=marker_line, alpha=1/cell_gen, label="") # plot the start of the cell trajectory with a blue circle (naive cell)
            end
        end
    end

    display(plt)
end



# average over a few simulations to see activation times
place_cells = false; # true to decide where cells are placed
num_sims = 5; # number of simulations (to average the activation time)

seed_ρ = 9319; # seed for scaffold generation 
ρ_1, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

# simulate indiviudal realisation
n_ABMxV, a_ABMxV, NI_ABM, cells_track_1, act_time_1, act_prop_1, avg_runtime_1 = ABM_track_cells();


# move micro-rods into cluster 
xpos_PDF(max) = Normal(max/4,max/10); # options:  Uniform(1,max), Normal(max/2,max/5), SkewNormal(max/2,max/5,-max/10), Arcsine(1,max), ...
ypos_PDF(max) = Normal(3max/4,max/10); # options cont        Exponential(max/5), Laplace(max/2), SymTriangularDist(max/2,max/5), TriangularDist(1,max,1)

ρ_2, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

# simulate indiviudal realisation
n_ABMxV, a_ABMxV, NI_ABM, cells_track_2, act_time_2, act_prop_2, avg_runtime_2 = ABM_track_cells();

xpos_PDF(max) = Uniform(1,max); # reset micro-rod position PDFs to default
ypos_PDF(max) = Uniform(1,max); 


# activation plot
plt_gens = 4;
plt_act_time_1 = act_time_1[1:plt_gens];
plt_act_time_2 = act_time_2[1:plt_gens];

plot_x = 250; plot_y = 400;
linewidth = 2.5;
timescale = 60*24;

plt_act = plot(ylims=[0,6], size=[plot_x,plot_y])
plot_x1 = [[0.6],[2.6],[4.6],[6.6]];
plot_x2 = [[1.4],[3.4],[5.4],[7.4]];

boxplot!(plot_x1, plt_act_time_1/timescale, color=scaf_colours[1], xticks=(1:2:7,1:4), label="")
boxplot!(plot_x2, plt_act_time_2/timescale, color=scaf_colours[2], label="")
plot!(xlabel="T cell generation", ylabel="Time to activation (days)")
display(plt_act)




## animations
# animate cell movement
pres_font = 25; # font size for presentation
dot_size = 7; # size of dots representing T cells

fps = 60; # animation frames per second
pltframe = 200; # plot every 'pltframe' frames

plt_empty = plot([], [], label=false, showaxis=false, axis=([], false)); # empty plot
anim = Animation();
@showprogress 1 "Generating expansion animation..." for n = 1:pltframe:Int(floor(T))+1
    # generate information for spatial plot
    naive_xind = zeros(0); naive_yind = zeros(0); # x and y indices for naive and activated cells
    active_xind = zeros(0); active_yind = zeros(0);
    for i = 1:N # for each node
        for j = 1:M 
            if n_plot[i,j,n] != 0
                append!(naive_yind,i);
                append!(naive_xind,j);
            end
            if a_plot[i,j,n] != 0
                append!(active_yind,i);
                append!(active_xind,j);
            end
        end
    end

    plt_spatial = heatmap(0:h:hei, 0:h:wid, reverse(ρ_1,dims=1), xlims=[0,hei], ylims=[0,wid], size=(spat_plot_x,spat_plot_y), colorbar=false, c=cgrad(:grayC, rev=true), clims=(0,maximum(ρ_1)), aspect_ratio=:equal, tickfont = (pres_font, "times")); # plot scaffold
    for tick = [100, 200, 300, 400]
        plot!(plt_spatial, [tick,tick], [0,8], c=:black, legend=false, linewidth=2)
        plot!(plt_spatial, [0,8], [tick,tick], c=:black, legend=false, linewidth=2)
    end
    for edge = [0, wid]
        plot!(plt_spatial, [edge,edge], [0,wid], c=:black, legend=false, linewidth=2)
        plot!(plt_spatial, [0,wid], [edge,edge], c=:black, legend=false, linewidth=2)
    end
    plot!(plt_spatial, naive_xind*h, naive_yind*h, seriestype=:scatter, mc=col_naive, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), label="", markersize=dot_size)
    plot!(plt_spatial, active_xind*h, active_yind*h, seriestype=:scatter, mc=col_active, alpha=0.7, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), label="", markersize=dot_size)


    plt_pop = plot(legend=false, xtickfontcolor=:white, xtickfont = 1, size=(pop_plot_x,pop_plot_y), leftmargin=30mm, ytickfont = (pres_font, "times"));
    plt_IL2 = plot(legend=false, size=(pop_plot_x,pop_plot_y), leftmargin=30mm, tickfont = (pres_font, "times"), yticks = 0:0.1:0.3);

    plot!(plt_pop, (0:T)/60/24, N_ABM, c=col_naive)
    plot!(plt_pop, (0:T)/60/24, A_ABM, c=col_active)
    plot!(plt_pop, (0:T)/60/24, NT_ABM, c=col_total)
    vline!([(n-1)/60/24], linestyle=:dash, c=col_total, linewidth=2, label="")

    plot!(plt_IL2, (0:T)/60/24, NI_ABM/IL2_scale, c=col_IL2)
    vline!([(n-1)/60/24], linestyle=:dash, c=col_total, linewidth=2, label="")


    plt = plot(plt_spatial, plt_pop, plt_IL2, size = (3*pop_plot_x,2*pop_plot_y), layout = (@layout([a grid(2,1)])), bottommargin=5mm)

    frame(anim, plt)
end
display(gif(anim, fps=fps))


# animate cell trajectories
fps = 60; # animation frames per second
pltframe = 200; # plot every 'pltframe' frames

line_step = 5; # plot the cell trajectory as a line which skips line_step timesteps

# first scaffold
cells_plot = [cells_track_1[1],cells_track_1[2]];
for i in eachindex(cells_plot)
    cells_plot[i][4][2,:] = abs.(cells_plot[i][4][2,:].-hei); # flip y coords for plotting
end

anim = Animation();
max_time = max(cells_track_1[1][2],cells_track_1[2][2],cells_track_2[1][2],cells_track_2[2][2])+12*60; 
if isnan(max_time) || max_time > T
    max_time = T;
end
@showprogress 1 "Generating trajectory animation 1..." for i = 1:pltframe:max_time
    plt = heatmap(0:h:hei, 0:h:wid, reverse(ρ_1,dims=1), xlims=[0,hei], ylims=[0,wid], size=(spat_plot_x,spat_plot_y), colorbar=false, c=cgrad(:grayC, rev=true), clims=(0,maximum(ρ_2)), aspect_ratio=:equal, xlab="x (μm)", ylab="y (μm)"); # plot scaffold
    
    # plot cell trajectory in blue (naive) or red (activated) with colour proportional to 1/generation
    for cell in reverse(cells_plot) # [time of birth, time of death (NaN if not), time of activation (NaN if not), (x,y) positions over time, generation]
        birth_i = cell[1]; # time index of birth
        if i >= birth_i # if cell was born by this time
            death_i = cell[2]; # time index of death (NaN if it never died)
            act_i = cell[3]; # time index of activation (NaN if it never got activated)
            cell_traj = cell[4]; # current cell's trajectory
            cell_gen = cell[5]^2; # current cell's generation squared (for plotting)
            if isnan(death_i)
                end_i = i;
            else
                end_i = min(i,death_i);
            end
            traj_path = birth_i:line_step:end_i; 
            if end_i == death_i
                traj_path = [traj_path;death_i];
            end
            cell_act = false; # keep track of when cell becomes activated
            for traj_i = length(traj_path)-1:-1:1 # for each plot time from the current time/death time back to their birth time
                # don't draw line if they cross through the boundary
                if maximum(abs.([diff(cell_traj[1,traj_path[traj_i]:traj_path[traj_i+1]]);diff(cell_traj[2,traj_path[traj_i]:traj_path[traj_i+1]])])) <= ε
                    if traj_path[traj_i] >= act_i
                        # plot trajectory with a red line (activated) 
                        plot!(cell_traj[1,traj_path[traj_i]:traj_path[traj_i+1]], cell_traj[2,traj_path[traj_i]:traj_path[traj_i+1]], c=col_active, alpha=traj_alpha, linewidth=traj_width, label="")
                    else
                        # plot trajectory with a blue line (naive) 
                        if !cell_act && !isnan(act_i) && traj_path[traj_i]>=act_i>=traj_path[traj_i+1] # if this is the time range in which the cell becomes activated, plot the segment when they were activated and then when they were naive
                            plot!(cell_traj[1,[traj_path[traj_i],act_i]], cell_traj[2,[traj_path[traj_i],act_i]], c=col_active, alpha=traj_alpha, linewidth=traj_width, label="")
                            plot!(cell_traj[1,[act_i,traj_path[traj_i+1]]], cell_traj[2,[act_i,traj_path[traj_i+1]]], c=col_naive, alpha=traj_alpha, linewidth=traj_width, label="")
                            cell_act = true
                        else
                            plot!(cell_traj[1,traj_path[traj_i]:traj_path[traj_i+1]], cell_traj[2,traj_path[traj_i]:traj_path[traj_i+1]], c=col_naive, alpha=traj_alpha, linewidth=traj_width, label="")
                        end
                    end
                end
            end
            if i >= death_i # if cell died by this time
                # plot the end with a black cross where it died
                plot!([cell_traj[1,death_i]], [cell_traj[2,death_i]], seriestype=:scatter, markershape=:x, c=:black, alpha=1/cell_gen, label="") 
                #=if isnan(act_i)
                    # plot the end with a blue cross (naive) where it died
                    plot!([cell_traj[1,death_i]], [cell_traj[2,death_i]], seriestype=:scatter, markershape=:x, c=col_naive, alpha=1/cell_gen, label="") 
                else
                    # plot the end with a red cross (activated) where it died
                    plot!([cell_traj[1,death_i]], [cell_traj[2,death_i]], seriestype=:scatter, markershape=:x, c=col_active, alpha=1/cell_gen, label="") 
                end=#
            end
            # plot the start of the cell trajectory with a blue circle (naive cell)
            plot!([cell_traj[1,birth_i]], [cell_traj[2,birth_i]], seriestype=:scatter, c=col_naive, alpha=1/cell_gen, label="") # plot the start of the cell trajectory with a blue circle (naive cell)
            # plot the current location of the cell
            if i < death_i
                if i >= act_i
                    plot!([cell_traj[1,i]], [cell_traj[2,i]], seriestype=:scatter, c=col_active, markersize=act_size, alpha=1/cell_gen, label="") # plot the start of the cell trajectory with a blue circle (naive cell)
                else
                    plot!([cell_traj[1,i]], [cell_traj[2,i]], seriestype=:scatter, c=col_naive, markersize=act_size, alpha=1/cell_gen, label="") # plot the start of the cell trajectory with a blue circle (naive cell)
                end
            end
        end
    end

    frame(anim, plot!(plt, title="𝑡 = $(round((i-1)/24/60,digits=1)) days"))
end
display(gif(anim, fps=fps))

# second scaffold
cells_plot = [cells_track_2[1],cells_track_2[2]];
for i in eachindex(cells_plot)
    cells_plot[i][4][2,:] = abs.(cells_plot[i][4][2,:].-hei); # flip y coords for plotting
end

anim = Animation();
max_time = max(cells_track_1[1][2],cells_track_1[2][2],cells_track_2[1][2],cells_track_2[2][2])+12*60; 
if isnan(max_time) || max_time > T
    max_time = T;
end
@showprogress 1 "Generating trajectory animation 2..." for i = 1:pltframe:max_time
    plt = heatmap(0:h:hei, 0:h:wid, reverse(ρ_2,dims=1), xlims=[0,hei], ylims=[0,wid], size=(spat_plot_x,spat_plot_y), colorbar=false, c=cgrad(:grayC, rev=true), clims=(0,maximum(ρ_2)), aspect_ratio=:equal, xlab="x (μm)", ylab="y (μm)"); # plot scaffold
    
    # plot cell trajectory in blue (naive) or red (activated) with colour proportional to 1/generation
    for cell in reverse(cells_plot) # [time of birth, time of death (NaN if not), time of activation (NaN if not), (x,y) positions over time, generation]
        birth_i = cell[1]; # time index of birth
        if i >= birth_i # if cell was born by this time
            death_i = cell[2]; # time index of death (NaN if it never died)
            act_i = cell[3]; # time index of activation (NaN if it never got activated)
            cell_traj = cell[4]; # current cell's trajectory
            cell_gen = cell[5]^2; # current cell's generation squared (for plotting)
            if isnan(death_i)
                end_i = i;
            else
                end_i = min(i,death_i);
            end
            traj_path = birth_i:line_step:end_i; 
            if end_i == death_i
                traj_path = [traj_path;death_i];
            end
            cell_act = false; # keep track of when cell becomes activated
            for traj_i = length(traj_path)-1:-1:1 # for each plot time from the current time/death time back to their birth time
                # don't draw line if they cross through the boundary
                if maximum(abs.([diff(cell_traj[1,traj_path[traj_i]:traj_path[traj_i+1]]);diff(cell_traj[2,traj_path[traj_i]:traj_path[traj_i+1]])])) <= ε
                    if traj_path[traj_i] >= act_i
                        # plot trajectory with a red line (activated) 
                        plot!(cell_traj[1,traj_path[traj_i]:traj_path[traj_i+1]], cell_traj[2,traj_path[traj_i]:traj_path[traj_i+1]], c=col_active, alpha=traj_alpha, linewidth=traj_width, label="")
                    else
                        # plot trajectory with a blue line (naive) 
                        if !cell_act && !isnan(act_i) && traj_path[traj_i]>=act_i>=traj_path[traj_i+1] # if this is the time range in which the cell becomes activated, plot the segment when they were activated and then when they were naive
                            plot!(cell_traj[1,[traj_path[traj_i],act_i]], cell_traj[2,[traj_path[traj_i],act_i]], c=col_active, alpha=traj_alpha, linewidth=traj_width, label="")
                            plot!(cell_traj[1,[act_i,traj_path[traj_i+1]]], cell_traj[2,[act_i,traj_path[traj_i+1]]], c=col_naive, alpha=traj_alpha, linewidth=traj_width, label="")
                            cell_act = true
                        else
                            plot!(cell_traj[1,traj_path[traj_i]:traj_path[traj_i+1]], cell_traj[2,traj_path[traj_i]:traj_path[traj_i+1]], c=col_naive, alpha=traj_alpha, linewidth=traj_width, label="")
                        end
                    end
                end
            end
            if i >= death_i # if cell died by this time
                # plot the end with a black cross where it died
                plot!([cell_traj[1,death_i]], [cell_traj[2,death_i]], seriestype=:scatter, markershape=:x, c=:black, alpha=1/cell_gen, label="") 
                #=if isnan(act_i)
                    # plot the end with a blue cross (naive) where it died
                    plot!([cell_traj[1,death_i]], [cell_traj[2,death_i]], seriestype=:scatter, markershape=:x, c=col_naive, alpha=1/cell_gen, label="") 
                else
                    # plot the end with a red cross (activated) where it died
                    plot!([cell_traj[1,death_i]], [cell_traj[2,death_i]], seriestype=:scatter, markershape=:x, c=col_active, alpha=1/cell_gen, label="") 
                end=#
            end
            # plot the start of the cell trajectory with a blue circle (naive cell)
            plot!([cell_traj[1,birth_i]], [cell_traj[2,birth_i]], seriestype=:scatter, c=col_naive, alpha=1/cell_gen, label="") # plot the start of the cell trajectory with a blue circle (naive cell)
            # plot the current location of the cell
            if i < death_i
                if i >= act_i
                    plot!([cell_traj[1,i]], [cell_traj[2,i]], seriestype=:scatter, c=col_active, markersize=act_size, alpha=1/cell_gen, label="") # plot the start of the cell trajectory with a blue circle (naive cell)
                else
                    plot!([cell_traj[1,i]], [cell_traj[2,i]], seriestype=:scatter, c=col_naive, markersize=act_size, alpha=1/cell_gen, label="") # plot the start of the cell trajectory with a blue circle (naive cell)
                end
            end
        end
    end

    frame(anim, plot!(plt, title="𝑡 = $(round((i-1)/24/60,digits=1)) days"))
end
display(gif(anim, fps=fps))




## test for significance 
using HypothesisTests
using DataFrames 

function make_dataframe(times, domain_label)
    rows = []
    for (gen_idx, gen_times) in enumerate(times)
        for t in gen_times
            push!(rows, (scaffold=domain_label, generation=gen_idx, activation_time=t))
        end
    end
    return DataFrame(rows)
end

df1 = make_dataframe(plt_act_time_1, "Homogeneous")
df2 = make_dataframe(plt_act_time_2, "Heterogeneous")

df = vcat(df1, df2)

using CategoricalArrays
df.scaffold = categorical(df.scaffold)
df.generation = categorical(df.generation)


# test normality with Shapiro Wilk test
group_data = df[(df.scaffold .== "Heterogeneous") .&& (df.generation .== 3), :activation_time]

sw_test = ShapiroWilkTest(group_data)
println(sw_test)
println("p-value = ", pvalue(sw_test))


# tests with Mann–Whitney U test
# test difference between scaffolds for each generation 
for g in unique(df.generation)
    data_subset = df[df.generation .== g, :]
    d1 = data_subset[data_subset.scaffold .== "Homogeneous", :activation_time]
    d2 = data_subset[data_subset.scaffold .== "Heterogeneous", :activation_time]
    
    println("Generation $g")
    println(MannWhitneyUTest(d1, d2))
end

# test difference between activaiton times in scaffolds for generations 2 and 3
for g in unique(df.scaffold)
    data_subset = df[df.scaffold .== g, :]
    d1 = data_subset[data_subset.generation .== 2, :activation_time]
    d2 = data_subset[data_subset.generation .== 3, :activation_time]
    
    println("Scaffold $g")
    println(MannWhitneyUTest(d1, d2))
end




## Fig 5 - comparison between averaged ABM and PDE model
# ====          ===          ===          =============          ============
# ====  ===============  =======  =====================  ====================
# ====  ===============  =======  =====================  ====================
# ====       ==========  =======  ===     =============         ============
# ====  ===============  =======  ======  =====================  ============
# ====  ===============  =======  ======  =====================  ============
# ====  ===========          ===          =============         =============

T = Int(floor(2*24*60)); # total simulation time (mins) as Int            -       run for >= 14 days to see plots at 7 and 14 days
ntotal = 100; # number of stochastic simulations to average over

seed_ρ = 9319; # seed for scaffold generation 
ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

# run ABM
n_ABM, a_ABM, I_ABM, n_tot_sims, a_tot_sims, I_tot_sims, runtime_ABM = simulate_ABM(); # outputs from ABM: average positions of naive and activated cells, average total mass/spatial concentration of IL-2, total number of naive and activated cells for each simulation, total mass of IL-2 for each simulation, computation time
T_ABM = n_ABM + a_ABM; # average positions of total T cells from ABM

# run PDE model 
n_PDE, a_PDE, I_PDE, runtime_PDE = simulate_PDE(); # outputs from PDE model: spatial concentrations of naive and activated cells, total mass/spatial concentration of IL-2, computation time
T_PDE = n_PDE + a_PDE; # spatial concentrations of total T cells from PDE model



# plot parameters
spat_plot_x = 300; spat_plot_y = spat_plot_x; # default size for spatial plots (x and y sizes)
pop_plot_x = spat_plot_x*2; pop_plot_y = spat_plot_y; # default size for population plots (x and y sizes)

T_colbar_scaling = 1e-4; # factor this scale out of the colorbar for displaying above, for T cells
I_colbar_scaling = 1e-9; # factor this scale out of the colorbar for displaying above, for IL-2

pop_scaling = 1e6; # factor this scale out of the total T cell population

t_plots = [min(7,T/24/60),min(14,T/24/60)]; # times at which plots are generated for each population (days)
t_plots = Int.(floor.(t_plots*24*60)); # convert to min and ensure whole numbers

linewidths = 3; # population plot linewidths


# arrays for plotting (specific time points, reversed dimensions, rescaled for colourbar and change spatial units to μm^2)
n_ABM_plot = reverse(n_ABM[:,:,t_plots.+1],dims=1)/T_colbar_scaling/h^2;
a_ABM_plot = reverse(a_ABM[:,:,t_plots.+1],dims=1)/T_colbar_scaling/h^2;

n_PDE_plot = reverse(n_PDE[:,:,t_plots.÷τ.+1],dims=1)/T_colbar_scaling/h^2;
a_PDE_plot = reverse(a_PDE[:,:,t_plots.÷τ.+1],dims=1)/T_colbar_scaling/h^2;
#I_PDE_plot = reverse(I_PDE[:,:,t_plots.+1],dims=1)/I_colbar_scaling/h^2; 

# mass arrays (with standard deviation for ABM results)
N_ABM = mean(n_tot_sims*result_scale/pop_scaling, dims=2)[:]; SD_N_ABM = std(n_tot_sims*result_scale/pop_scaling, dims=2)[:];
A_ABM = mean(a_tot_sims*result_scale/pop_scaling, dims=2)[:]; SD_A_ABM = std(a_tot_sims*result_scale/pop_scaling, dims=2)[:];
NT_ABM = N_ABM + A_ABM; SD_NT_ABM = std((n_tot_sims+a_tot_sims)*result_scale/pop_scaling, dims=2)[:];
NI_ABM = mean(I_tot_sims*result_scale, dims=2)[:]; SD_NI_ABM = std(I_tot_sims*result_scale, dims=2)[:];

N_PDE = (sum(n_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(n_PDE[1,2:M-1,:],dims=1)+sum(n_PDE[N,2:M-1,:],dims=1)+sum(n_PDE[2:N-1,1,:],dims=1)+sum(n_PDE[2:N-1,M,:],dims=1))[:]/2+(n_PDE[1,1,:]+n_PDE[1,M,:]+n_PDE[N,1,:]+n_PDE[N,M,:])/4)*result_scale/pop_scaling;
A_PDE = (sum(a_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(a_PDE[1,2:M-1,:],dims=1)+sum(a_PDE[N,2:M-1,:],dims=1)+sum(a_PDE[2:N-1,1,:],dims=1)+sum(a_PDE[2:N-1,M,:],dims=1))[:]/2+(a_PDE[1,1,:]+a_PDE[1,M,:]+a_PDE[N,1,:]+a_PDE[N,M,:])/4)*result_scale/pop_scaling;
NT_PDE = N_PDE + A_PDE;
NI_PDE = I_PDE*result_scale;

n_max = maximum(n_PDE_plot); # maximum value shown for naive T cell concentration
a_max = maximum(a_PDE_plot); # maximum value shown for activated T cell concentration


# spatial plots 
for n in eachindex(t_plots) # for each time point
    plt_n_ABM = heatmap(0:h:hei, 0:h:wid, n_ABM_plot[:,:,n], color=:blues, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), clims=(0,n_max), colorbar=false, size = (spat_plot_x,spat_plot_y), xlabel="𝑥 (μm)", ylabel="𝑦 (μm)", title="𝑡 = $(round(t_plots[n]/24/60,digits=1)) days");
    plt_a_ABM = heatmap(0:h:hei, 0:h:wid, a_ABM_plot[:,:,n], color=:reds, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), clims=(0,a_max), colorbar=false, size = (spat_plot_x,spat_plot_y), xlabel="𝑥 (μm)", ylabel="𝑦 (μm)", title="𝑡 = $(round(t_plots[n]/24/60,digits=1)) days");

    plt_n_PDE = heatmap(0:h:hei, 0:h:wid, n_PDE_plot[:,:,n], color=:blues, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), clims=(0,n_max), colorbar=false, size = (spat_plot_x,spat_plot_y), xlabel="𝑥 (μm)", ylabel="𝑦 (μm)", title="𝑡 = $(round(t_plots[n]/24/60,digits=1)) days");
    plt_a_PDE = heatmap(0:h:hei, 0:h:wid, a_PDE_plot[:,:,n], color=:reds, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), clims=(0,a_max), colorbar=false, size = (spat_plot_x,spat_plot_y), xlabel="𝑥 (μm)", ylabel="𝑦 (μm)", title="𝑡 = $(round(t_plots[n]/24/60,digits=1)) days");

    # add tick marks
    for plt = [plt_n_ABM, plt_a_ABM, plt_n_PDE, plt_a_PDE]
        for tick = [100,200,300,400]
            plot!(plt, [tick,tick], [0,8], c=:black, legend=false, linewidth=1.5)
            plot!(plt, [0,8], [tick,tick], c=:black, legend=false, linewidth=1.5)
        end
    end

    display(plt_n_ABM)
    display(plt_a_ABM)

    display(plt_n_PDE)
    display(plt_a_PDE)
end


# population plots
plt_pop = plot(legend=false, ylims=(0,maximum([N_ABM+SD_N_ABM; A_ABM+SD_A_ABM; NT_ABM+SD_NT_ABM])), size=(pop_plot_x,pop_plot_y), xlabel="𝑡 (days)", ylabel="Number of T cells (millions)");
plt_IL2 = plot(legend=false, ylims=(0,maximum(NI_ABM+SD_NI_ABM)), size=(pop_plot_x,pop_plot_y), xlabel="𝑡 (days)", ylabel="Mass of IL-2 (ng)");

plot!(plt_pop, (0:Int(floor(T)))/60/24, N_ABM, ribbon=SD_N_ABM, fillalpha = 0.5, c=col_naive, linewidth=linewidths)
plot!(plt_pop, (0:τ:Int(floor(T)))/60/24, N_PDE, c=col_naive, linestyle=:dash, linewidth=linewidths)

plot!(plt_pop, (0:Int(floor(T)))/60/24, A_ABM, ribbon=SD_A_ABM, fillalpha = 0.5, c=col_active, linewidth=linewidths)
plot!(plt_pop, (0:τ:Int(floor(T)))/60/24, A_PDE, c=col_active, linestyle=:dash, linewidth=linewidths)

plot!(plt_pop, (0:Int(floor(T)))/60/24, NT_ABM, ribbon=SD_NT_ABM, fillalpha = 0.5, c=col_total, linewidth=linewidths)
plot!(plt_pop, (0:τ:Int(floor(T)))/60/24, NT_PDE, c=col_total, linestyle=:dash, linewidth=linewidths)

plot!(plt_IL2, (0:Int(floor(T)))/60/24, NI_ABM, ribbon=SD_NI_ABM, fillalpha = 0.5, c=col_IL2, linewidth=linewidths)
plot!(plt_IL2, (0:τ:Int(floor(T)))/60/24, NI_PDE, c=col_IL2, linestyle=:dash, linewidth=linewidths)

# display vertical lines to indicate time points in spatial plots 
vline_width = 1.5; 
for t_n in t_plots/(24*60) # for each time point (in days)
    vline!(plt_pop, [t_n], c=col_total, linestyle=:dash, linewidth=vline_width)
    vline!(plt_IL2, [t_n], c=col_total, linestyle=:dash, linewidth=vline_width)
end

display(plt_pop)
display(plt_IL2)




## animation of average T cell concentrations and expansion 
pres_font = 25; # font size for presentation
dot_size = 7; # size of dots representing T cells

fps = 60; # animation frames per second
pltframe = 100; # plot every 'pltframe' frames

# arrays for plotting (reversed dimensions, rescaled for colourbar and change spatial units to μm^2)
n_ABM_plot = reverse(n_ABM,dims=1)/T_colbar_scaling/h^2;
a_ABM_plot = reverse(a_ABM,dims=1)/T_colbar_scaling/h^2;
n_PDE_plot = reverse(n_PDE,dims=1)/T_colbar_scaling/h^2;
a_PDE_plot = reverse(a_PDE,dims=1)/T_colbar_scaling/h^2;

n_max = maximum(n_PDE_plot); # maximum value shown for naive T cell concentration
a_max = maximum(a_PDE_plot); # maximum value shown for activated T cell concentration
tot_max = maximum(n_ABM_plot+a_ABM_plot); # maximum value shown for total T cell concentration

n_ABM_anim = Animation();
n_PDE_anim = Animation();
a_ABM_anim = Animation();
a_PDE_anim = Animation();
#IL2_anim = Animation();
tot_ABM_anim = Animation();
pop_anim = Animation();
IL2pop_anim = Animation();
@showprogress 1 "Generating animation..." for n = 1:pltframe:Int(floor(T))+1
    plt_n_ABM = heatmap(0:h:hei, 0:h:wid, n_ABM_plot[:,:,n], color=:blues, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), clims=(0,n_max), colorbar=false, size = (spat_plot_x,spat_plot_y), xlabel="𝑥 (μm)", ylabel="𝑦 (μm)", title="𝑡 = $(round((n-1)/24/60,digits=1)) days");
    plt_a_ABM = heatmap(0:h:hei, 0:h:wid, a_ABM_plot[:,:,n], color=:reds, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), clims=(0,a_max), colorbar=false, size = (spat_plot_x,spat_plot_y), xlabel="𝑥 (μm)", ylabel="𝑦 (μm)", title="𝑡 = $(round((n-1)/24/60,digits=1)) days");

    plt_n_PDE = heatmap(0:h:hei, 0:h:wid, n_PDE_plot[:,:,n], color=:blues, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), clims=(0,n_max), colorbar=false, size = (spat_plot_x,spat_plot_y), xlabel="𝑥 (μm)", ylabel="𝑦 (μm)", title="𝑡 = $(round((n-1)/24/60,digits=1)) days");
    plt_a_PDE = heatmap(0:h:hei, 0:h:wid, a_PDE_plot[:,:,n], color=:reds, aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), clims=(0,a_max), colorbar=false, size = (spat_plot_x,spat_plot_y), xlabel="𝑥 (μm)", ylabel="𝑦 (μm)", title="𝑡 = $(round((n-1)/24/60,digits=1)) days");

    plt_tot_ABM = heatmap(0:h:hei, 0:h:wid, n_ABM_plot[:,:,n]+a_ABM_plot[:,:,n], color=cgrad(:Greys_3, rev=true), aspect_ratio=:equal, xlimits=(0,wid), ylimits=(0,hei), clims=(0,tot_max), colorbar=false, size = (spat_plot_x,spat_plot_y), xlabel="𝑥 (μm)", ylabel="𝑦 (μm)", title="𝑡 = $(round((n-1)/24/60,digits=1)) days");

    # add tick marks
    for plt = [plt_n_ABM, plt_a_ABM, plt_n_PDE, plt_a_PDE, plt_tot_ABM]
        for tick = [100,200,300,400]
            plot!(plt, [tick,tick], [0,8], c=:black, legend=false, linewidth=1.5)
            plot!(plt, [0,8], [tick,tick], c=:black, legend=false, linewidth=1.5)
        end
    end

    frame(n_ABM_anim, plt_n_ABM)
    frame(n_PDE_anim, plt_n_PDE)
    frame(a_ABM_anim, plt_a_ABM)
    frame(a_PDE_anim, plt_a_PDE)
    #frame(IL2_anim, plt_I)

    frame(tot_ABM_anim, plt_tot_ABM)


    # population plots 
    plt_pop = plot(legend=false, ylims=(0,maximum([N_ABM+SD_N_ABM; A_ABM+SD_A_ABM; NT_ABM+SD_NT_ABM])), size=(pop_plot_x,pop_plot_y), xtickfontcolor=:white, xlabel="𝑡 (days)", ylabel="Number of T cells (millions)");
    plt_IL2 = plot(legend=false, ylims=(0,maximum(NI_ABM+SD_NI_ABM)), size=(pop_plot_x,pop_plot_y), xlabel="𝑡 (days)", ylabel="Mass of IL-2 (ng)");

    plot!(plt_pop, (0:Int(floor(T)))/60/24, N_ABM, ribbon=SD_N_ABM, fillalpha = 0.5, c=col_naive, linewidth=linewidths)
    plot!(plt_pop, (0:Int(floor(T)))/60/24, N_PDE, c=col_naive, linestyle=:dash, linewidth=linewidths)
    plot!(plt_pop, (0:Int(floor(T)))/60/24, A_ABM, ribbon=SD_A_ABM, fillalpha = 0.5, c=col_active, linewidth=linewidths)
    plot!(plt_pop, (0:Int(floor(T)))/60/24, A_PDE, c=col_active, linestyle=:dash, linewidth=linewidths)
    plot!(plt_pop, (0:Int(floor(T)))/60/24, NT_ABM, ribbon=SD_NT_ABM, fillalpha = 0.5, c=col_total, linewidth=linewidths)
    plot!(plt_pop, (0:Int(floor(T)))/60/24, NT_PDE, c=col_total, linestyle=:dash, linewidth=linewidths)
    vline!([(n-1)/60/24], linestyle=:dash, c=col_total, linewidth=2, label="")

    plot!(plt_IL2, (0:Int(floor(T)))/60/24, NI_ABM, ribbon=SD_NI_ABM, fillalpha = 0.5, c=col_IL2, linewidth=linewidths)
    plot!(plt_IL2, (0:Int(floor(T)))/60/24, NI_PDE, c=col_IL2, linestyle=:dash, linewidth=linewidths)
    vline!([(n-1)/60/24], linestyle=:dash, c=col_total, linewidth=2, label="")


    frame(pop_anim, plt_pop)
    frame(IL2pop_anim, plt_IL2)
end

display(gif(n_ABM_anim, fps=fps))
display(gif(n_PDE_anim, fps=fps))
display(gif(a_ABM_anim, fps=fps))
display(gif(a_PDE_anim, fps=fps))

display(gif(tot_ABM_anim, fps=fps))

display(gif(pop_anim, fps=fps))
display(gif(IL2pop_anim, fps=fps))





## Fig 6 - comparison between ABM, PDE and ODE models, and convergence of PDE to ODE model under blurred scaffolds
# ====          ===          ===          ==============         ============
# ====  ===============  =======  =====================  ====================
# ====  ===============  =======  =====================  ====================
# ====       ==========  =======  ===     =============         =============
# ====  ===============  =======  ======  =============  ======  ============
# ====  ===============  =======  ======  =============  ======  ============
# ====  ===========          ===          ==============         ============

use_ODE_for_IL2 = true;
seed_ρ = 9319;

T = Int(floor(0.5*24*60)); # total simulation time (mins) as Int            -       run for >= 14 days to see plots at 7 and 14 days
Ts = T;
ntotal = 100; # number of stochastic simulations to average over

conc_arr = [33, 333];# [10, 33, 333, 1000]; # concentrations of solutions containing micro-rods (μg/mL)

ABM_result = zeros(T+1,ntotal,length(conc_arr),3); # [time, simulation, concentration, naive/activated/IL-2]
PDE_result = zeros(T+1,length(conc_arr),3); # [time, concentration, naive/activated/IL-2]
ODE_result = zeros(T+1,length(conc_arr),3); # [time, concentration, naive/activated/IL-2]

for scaff = eachindex(conc_arr)
    APC_ms_conc = conc_arr[scaff];
    m_input = APC_ms_conc*APC_ms_vol; # input mass of MSRs (μg) assuming majority of APC-ms mass is MSRs       m_MSR*result_scale*n for n micro-rods
    ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

    N_rods = mean(ρ)*N*M*result_scale; # APC_ms_conc*APC_ms_vol/m_MSR; # total number of rods in the domain (integral of ρ(x,y))

    I_loaded = IL2_ratio*m_input; # initial mass of loaded cytokine (retained inside MSRs with total mass m_input) (ng) 
    I0_mass = I_loaded/result_scale; # initial mass of free cytokine if non-zero IC (ng)

    # mathematical parameters
    T_loaded = 10*result_scale; # total number of loaded T cells
    N0 = Int(floor(T_loaded/result_scale)); # initial number of T cells in scaled domain
    halfN0 = Int(floor(T_loaded/result_scale/2)); # half of initial number of cells as an integer

    print("Scaffold $(scaff) of $(length(conc_arr))")

    # run ABM
    n_ABM, a_ABM, I_ABM, n_tot_sims, a_tot_sims, I_tot_sims, runtime_ABM = simulate_ABM(); # outputs from ABM: average positions of naive and activated cells, average total mass/spatial concentration of IL-2, total number of naive and activated cells for each simulation, total mass of IL-2 for each simulation, computation time
    T_ABM = n_ABM + a_ABM; # average positions of total T cells from ABM
    ABM_result[:,:,scaff,1] = n_tot_sims;
    ABM_result[:,:,scaff,2] = a_tot_sims;
    ABM_result[:,:,scaff,3] = I_tot_sims;

    # run PDE model 
    n_PDE, a_PDE, I_PDE, runtime_PDE = simulate_PDE(); # outputs from PDE model: spatial concentrations of naive and activated cells, total mass/spatial concentration of IL-2, computation time
    T_PDE = n_PDE + a_PDE; # spatial concentrations of total T cells from PDE model
    PDE_result[:,scaff,1] = (sum(n_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(n_PDE[1,2:M-1,:],dims=1)+sum(n_PDE[N,2:M-1,:],dims=1)+sum(n_PDE[2:N-1,1,:],dims=1)+sum(n_PDE[2:N-1,M,:],dims=1))[:]/2+(n_PDE[1,1,:]+n_PDE[1,M,:]+n_PDE[N,1,:]+n_PDE[N,M,:])/4)*result_scale;
    PDE_result[:,scaff,2] = (sum(a_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(a_PDE[1,2:M-1,:],dims=1)+sum(a_PDE[N,2:M-1,:],dims=1)+sum(a_PDE[2:N-1,1,:],dims=1)+sum(a_PDE[2:N-1,M,:],dims=1))[:]/2+(a_PDE[1,1,:]+a_PDE[1,M,:]+a_PDE[N,1,:]+a_PDE[N,M,:])/4)*result_scale;
    PDE_result[:,scaff,3] = I_PDE*result_scale;

    # run ODE model 
    n_ODE, a_ODE, I_ODE, runtime_ODE = simulate_ODE(); # outputs from ODE model: total number of naive and activated T cells and total mass of IL-2 (scaled up for the whole dish), computation time
    T_ODE = n_ODE + a_ODE; # mass of total T cells from ODE model
    ODE_result[:,scaff,1] = n_ODE;
    ODE_result[:,scaff,2] = a_ODE;
    ODE_result[:,scaff,3] = I_ODE;
end
ABM_result = ABM_result*result_scale;

# set parameters back to default 
APC_ms_conc = 333; # concentration of inputted APC-ms (μg/mL) 
APC_ms_vol = 100e-3; # input volume of APC-ms (mL) 
m_input = APC_ms_conc*APC_ms_vol; # input mass of MSRs (μg) assuming majority of APC-ms mass is MSRs       m_MSR*result_scale*n for n micro-rods
I_loaded = IL2_ratio*m_input; # initial mass of loaded cytokine (retained inside MSRs with total mass m_input) (ng) 
I0_mass = I_loaded/result_scale; # initial mass of free cytokine if non-zero IC (ng)


pop_scale = 10^6; # value to factor out of the y axis for plotting
ABM_result[:,:,:,1:2] = ABM_result[:,:,:,1:2]/pop_scale; 
PDE_result[:,:,1:2] = PDE_result[:,:,1:2]/pop_scale; 
ODE_result[:,:,1:2] = ODE_result[:,:,1:2]/pop_scale; 

mean_ABM = mean(ABM_result, dims=2)[:,1,:,:]; # mean and standard deviation for all cell types and scaffolds over time
SD_ABM = std(ABM_result, dims=2)[:,1,:,:]; 
SD_NT_ABM_scaff = std(ABM_result[:,:,:,1]+ABM_result[:,:,:,2], dims=2)[:,1,:]; # standard deviation for total T cells for different scaffolds


# plot parameters 
spat_plot_x = 250; spat_plot_y = spat_plot_x; # default size for spatial plots (x and y sizes)
pop_plot_x = spat_plot_x*1.4; pop_plot_y = spat_plot_y; # default size for population plots (x and y sizes)

ODE_linestyle = :dashdot; # linestyle for ODE/MFA results

linewidths = 3; # line widths for population plots

# y axis limits for each population
naive_lims = [0, maximum(mean_ABM[:,:,1]+SD_ABM[:,:,1])]; # cells
active_lims = [0, maximum(mean_ABM[:,:,2]+SD_ABM[:,:,2])]; # cells
total_lims = [0, maximum(mean_ABM[:,:,1]+mean_ABM[:,:,2]+SD_NT_ABM_scaff)]; # cells
IL2_lims = [0, maximum(mean_ABM[:,:,3]+SD_ABM[:,:,3])]; # ng


# for each scaffold concentration:
for scaff = eachindex(conc_arr)
    # extract each result
    N_ABM = mean_ABM[:,scaff,1]; # average total amounts from ABM
    A_ABM = mean_ABM[:,scaff,2];
    NT_ABM = N_ABM + A_ABM;
    NI_ABM = mean_ABM[:,scaff,3];

    SD_N_ABM = SD_ABM[:,scaff,1]; # standard deviation of total amounts from ABM
    SD_A_ABM = SD_ABM[:,scaff,2];
    SD_NT_ABM = SD_NT_ABM_scaff[:,scaff];
    SD_NI_ABM = SD_ABM[:,scaff,3];

    N_PDE = PDE_result[:,scaff,1]; # total amounts from PDE
    A_PDE = PDE_result[:,scaff,2];
    NT_PDE = N_PDE + A_PDE;
    NI_PDE = PDE_result[:,scaff,3];

    N_ODE = ODE_result[:,scaff,1]; # total amounts from ODE (independent of scaffold)
    A_ODE = ODE_result[:,scaff,2];
    NT_ODE = N_ODE + A_ODE;
    NI_ODE = ODE_result[:,scaff,3];

    # re-generate scaffold 
    ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(conc_arr[scaff], APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

    # plot results
    plt_scaff = heatmap(0:h:hei, 0:h:wid, reverse(ρ,dims=1), xlims=[0,hei], ylims=[0,wid], size=(spat_plot_x,spat_plot_y), colorbar=false, c=cgrad(:grayC, rev=true), clims=(0,ρ_max), aspect_ratio=:equal); # plot scaffold
    for tick = [100, 200, 300, 400]
        plot!(plt_scaff, [tick,tick], [0,8], c=:black, legend=false, linewidth=1.5)
        plot!(plt_scaff, [0,8], [tick,tick], c=:black, legend=false, linewidth=1.5)
    end
    for edge = [0, wid]
        plot!(plt_scaff, [edge,edge], [0,wid], c=:black, legend=false, linewidth=1.5)
        plot!(plt_scaff, [0,wid], [edge,edge], c=:black, legend=false, linewidth=1.5)
    end

    plt_naive = plot(legend=false, size=(pop_plot_x,pop_plot_y), ylims=naive_lims);
    plt_active = plot(legend=false, size=(pop_plot_x,pop_plot_y), ylims=active_lims);
    plt_total = plot(legend=false, size=(pop_plot_x,pop_plot_y), ylims=total_lims);
    plt_IL2 = plot(legend=false, size=(pop_plot_x,pop_plot_y), ylims=IL2_lims);

    plot!(plt_naive, (0:Int(floor(T)))/60/24, N_ABM, ribbon=SD_N_ABM, fillalpha = 0.5, c=col_naive, linewidth=linewidths)
    plot!(plt_active, (0:Int(floor(T)))/60/24, A_ABM, ribbon=SD_A_ABM, fillalpha = 0.5, c=col_active, linewidth=linewidths)
    plot!(plt_total, (0:Int(floor(T)))/60/24, NT_ABM, ribbon=SD_NT_ABM, fillalpha = 0.5, c=col_total, linewidth=linewidths)
    plot!(plt_IL2, (0:Int(floor(T)))/60/24, NI_ABM, ribbon=SD_NI_ABM, fillalpha = 0.5, c=col_IL2, linewidth=linewidths)

    plot!(plt_naive, (0:Int(floor(T)))/60/24, N_PDE, c=col_naive, linestyle=:dash, linewidth=linewidths)
    plot!(plt_active, (0:Int(floor(T)))/60/24, A_PDE, c=col_active, linestyle=:dash, linewidth=linewidths)
    plot!(plt_total, (0:Int(floor(T)))/60/24, NT_PDE, c=col_total, linestyle=:dash, linewidth=linewidths)
    plot!(plt_IL2, (0:Int(floor(T)))/60/24, NI_PDE, c=col_IL2, linestyle=:dash, linewidth=linewidths)

    plot!(plt_naive, (0:Int(floor(T)))/60/24, N_ODE, c=col_naive, linestyle=ODE_linestyle, linewidth=linewidths)
    plot!(plt_active, (0:Int(floor(T)))/60/24, A_ODE, c=col_active, linestyle=ODE_linestyle, linewidth=linewidths)
    plot!(plt_total, (0:Int(floor(T)))/60/24, NT_ODE, c=col_total, linestyle=ODE_linestyle, linewidth=linewidths)
    plot!(plt_IL2, (0:Int(floor(T)))/60/24, NI_ODE, c=col_IL2, linestyle=ODE_linestyle, linewidth=linewidths)

    display(plt_scaff)
    display(plt_naive)
    display(plt_active)
    display(plt_total)
    display(plt_IL2)
end



## generate data for the end of Fig 6 (6e,f):  solve PDE model for incrementally blurred scaffolds (keeps the same seed for scaffold generation)
ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

# run ODE model 
n_ODE, a_ODE, I_ODE, runtime_ODE = simulate_ODE(); # outputs from ODE model: total number of naive and activated T cells and total mass of IL-2 (scaled up for the whole dish), computation time
T_ODE = n_ODE + a_ODE; # mass of total T cells from ODE model

blurring_arr = [1,5,20,50,N];# [1:2:21;25:4:45;53:8:N;N] #1:2:N; # array containing side length (pixels) of the uniform kernels used to convolute (blur) the scaffold image

T_blur = zeros(Ts+1, length(blurring_arr)); # array containing the total T cells over time for different blurring scales
I_blur = zeros(Ts+1, length(blurring_arr)); # array containing the mass of cytokine over time for different blurring scales

for blur_n = eachindex(blurring_arr) # for each blurring scale
    ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring_arr[blur_n], seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)
    
    # run PDE model 
    print("Scaffold $(blur_n) of $(length(blurring_arr))\n")
    n_PDE, a_PDE, I_PDE, runtime_PDE = simulate_PDE(); # outputs from PDE model: spatial concentrations of naive and activated cells, total mass/spatial concentration of IL-2, computation time
    T_PDE = n_PDE + a_PDE; # spatial concentrations of total T cells from PDE model

    # calculate mass of total T cells and cytokine and store in arrays
    n_mass_PDE = (sum(n_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(n_PDE[1,2:M-1,:],dims=1)+sum(n_PDE[N,2:M-1,:],dims=1)+sum(n_PDE[2:N-1,1,:],dims=1)+sum(n_PDE[2:N-1,M,:],dims=1))[:]/2+(n_PDE[1,1,:]+n_PDE[1,M,:]+n_PDE[N,1,:]+n_PDE[N,M,:])/4)*result_scale; # naive cells from PDE model
    a_mass_PDE = (sum(a_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(a_PDE[1,2:M-1,:],dims=1)+sum(a_PDE[N,2:M-1,:],dims=1)+sum(a_PDE[2:N-1,1,:],dims=1)+sum(a_PDE[2:N-1,M,:],dims=1))[:]/2+(a_PDE[1,1,:]+a_PDE[1,M,:]+a_PDE[N,1,:]+a_PDE[N,M,:])/4)*result_scale; # activated cells from PDE model
    T_blur[:,blur_n] = n_mass_PDE + a_mass_PDE; # total cells from PDE model

    if use_ODE_for_IL2
        I_blur[:,blur_n] = I_PDE*result_scale;
    else
        I_blur[:,blur_n] = (sum(I_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(I_PDE[1,2:M-1,:],dims=1)+sum(I_PDE[N,2:M-1,:],dims=1)+sum(I_PDE[2:N-1,1,:],dims=1)+sum(I_PDE[2:N-1,M,:],dims=1))[:]/2+(I_PDE[1,1,:]+I_PDE[1,M,:]+I_PDE[N,1,:]+I_PDE[N,M,:])/4)*result_scale; # cytokine from PDE model
    end
end


## plotting parameters
pop_plot_x = 450; pop_plot_y = 300; # default size for population plots (x and y sizes)

# compute RMSE between CL and MFA 
RMSE_arr = zeros(length(blurring_arr)); # contains the RMSE between the CL and MFA for increasingly blurred scaffolds
for scaf = eachindex(blurring_arr)
    RMSE_arr[scaf] = rmse(T_ODE, T_blur[:,scaf])/pop_scale;
end

# plot RMSE 
plt_RMSE = plot(legend=false, size=(pop_plot_x,pop_plot_y));
plot!(plt_RMSE, blurring_arr*h, RMSE_arr, c=col_total);

display(plt_RMSE)




## Fig 7 - impact of IL-2 loading on T cell expansion, using PDE model
# ====          ===          ===          =============          ============
# ====  ===============  =======  ============================  =============
# ====  ===============  =======  ===========================  ==============
# ====       ==========  =======  ===     ==================  ===============
# ====  ===============  =======  ======  =================  ================
# ====  ===============  =======  ======  ================  =================
# ====  ===========          ===          ===============  ==================

T = Int(floor(3*24*60));
Ts = T;

seed_ρ = 30; 
ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

# run PDE model for default IL-2 secretion
n_PDE, a_PDE, I_PDE, runtime_PDE = simulate_PDE(); # outputs from PDE model: spatial concentrations of naive and activated cells, total mass/spatial concentration of IL-2, computation time

n_PDE_sec = n_PDE; # solutions when IL-2 is secreted
a_PDE_sec = a_PDE;
I_PDE_sec = I_PDE;


# simulation without secretion and supplementation only at the start
cyt_IC = "uniform"; # start with some IL-2
S(ρ,t) = 0; # zero IL-2 secretion from micro-rods  (ng/h^2/min)

# cytokine initial condition
I0 = I0_mass/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4); # concentration of particles at each point (conc = mass/volume)
IC_I = zeros(N,M) .+ I0; # setting initial condition

# run PDE model
n_PDE, a_PDE, I_PDE, runtime_PDE = simulate_PDE(); # outputs from PDE model: spatial concentrations of naive and activated cells, total mass/spatial concentration of IL-2, computation time

n_PDE_supp1 = n_PDE; # solutions when IL-2 is supplemented at the start
a_PDE_supp1 = a_PDE;
I_PDE_supp1 = I_PDE;


# simulation without secretion and supplementation at the days indicated by supp_times
supplement_IL2 = true; # true to supplement IL-2
supp_IL2_time = Int.(floor.([0,2,4,6]*24*60)); # times at which IL-2 is supplemented
I0_mass = I_loaded/result_scale/length(supp_IL2_time); # initial mass of free cytokine if non-zero IC (ng)

# cytokine initial condition
I0 = I0_mass/((N-2)*(M-2)+(2*(N-2)+2*(M-2))/2+4/4); # concentration of particles at each point (conc = mass/volume), divided by the number of supplementations
IC_I = zeros(N,M) .+ I0; # setting initial condition

supp_IL2_mass = [I0_mass for i in eachindex(supp_IL2_time)];

# run PDE model
n_PDE, a_PDE, I_PDE, runtime_PDE = simulate_PDE(); # outputs from PDE model: spatial concentrations of naive and activated cells, total mass/spatial concentration of IL-2, computation time

n_PDE_supp2 = n_PDE; # solutions when IL-2 is supplemented
a_PDE_supp2 = a_PDE;
I_PDE_supp2 = I_PDE;


# reset parameters 
cyt_IC = "zero";
I0_mass = I_loaded/result_scale; # initial mass of free cytokine if non-zero IC (ng)
supplement_IL2 = false; # true to supplement IL-2
S(ρ,t) = ρ*k*m_MSR/m_input*I_loaded*exp(-k*t)#*K_S; # density- and time-dependent source term for cytokine (IL-2) (time is t=τ*n, where n is the discretised time)  (ng/h^2/min)
IC_I = zeros(N,M); # set back to zero




pop_scale = 1e6; # factor this scale out of the T cell population for plotting
# total cells / mass of IL-2
N_PDE_sec = (sum(n_PDE_sec[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(n_PDE_sec[1,2:M-1,:],dims=1)+sum(n_PDE_sec[N,2:M-1,:],dims=1)+sum(n_PDE_sec[2:N-1,1,:],dims=1)+sum(n_PDE_sec[2:N-1,M,:],dims=1))[:]/2+(n_PDE_sec[1,1,:]+n_PDE_sec[1,M,:]+n_PDE_sec[N,1,:]+n_PDE_sec[N,M,:])/4)*result_scale/pop_scale;
A_PDE_sec = (sum(a_PDE_sec[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(a_PDE_sec[1,2:M-1,:],dims=1)+sum(a_PDE_sec[N,2:M-1,:],dims=1)+sum(a_PDE_sec[2:N-1,1,:],dims=1)+sum(a_PDE_sec[2:N-1,M,:],dims=1))[:]/2+(a_PDE_sec[1,1,:]+a_PDE_sec[1,M,:]+a_PDE_sec[N,1,:]+a_PDE_sec[N,M,:])/4)*result_scale/pop_scale;
NT_PDE_sec = N_PDE_sec + A_PDE_sec;
NI_PDE_sec = I_PDE_sec*result_scale;
N_PDE_supp1 = (sum(n_PDE_supp1[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(n_PDE_supp1[1,2:M-1,:],dims=1)+sum(n_PDE_supp1[N,2:M-1,:],dims=1)+sum(n_PDE_supp1[2:N-1,1,:],dims=1)+sum(n_PDE_supp1[2:N-1,M,:],dims=1))[:]/2+(n_PDE_supp1[1,1,:]+n_PDE_supp1[1,M,:]+n_PDE_supp1[N,1,:]+n_PDE_supp1[N,M,:])/4)*result_scale/pop_scale;
A_PDE_supp1 = (sum(a_PDE_supp1[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(a_PDE_supp1[1,2:M-1,:],dims=1)+sum(a_PDE_supp1[N,2:M-1,:],dims=1)+sum(a_PDE_supp1[2:N-1,1,:],dims=1)+sum(a_PDE_supp1[2:N-1,M,:],dims=1))[:]/2+(a_PDE_supp1[1,1,:]+a_PDE_supp1[1,M,:]+a_PDE_supp1[N,1,:]+a_PDE_supp1[N,M,:])/4)*result_scale/pop_scale;
NT_PDE_supp1 = N_PDE_supp1 + A_PDE_supp1;
NI_PDE_supp1 = I_PDE_supp1*result_scale;
N_PDE_supp2 = (sum(n_PDE_supp2[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(n_PDE_supp2[1,2:M-1,:],dims=1)+sum(n_PDE_supp2[N,2:M-1,:],dims=1)+sum(n_PDE_supp2[2:N-1,1,:],dims=1)+sum(n_PDE_supp2[2:N-1,M,:],dims=1))[:]/2+(n_PDE_supp2[1,1,:]+n_PDE_supp2[1,M,:]+n_PDE_supp2[N,1,:]+n_PDE_supp2[N,M,:])/4)*result_scale/pop_scale;
A_PDE_supp2 = (sum(a_PDE_supp2[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(a_PDE_supp2[1,2:M-1,:],dims=1)+sum(a_PDE_supp2[N,2:M-1,:],dims=1)+sum(a_PDE_supp2[2:N-1,1,:],dims=1)+sum(a_PDE_supp2[2:N-1,M,:],dims=1))[:]/2+(a_PDE_supp2[1,1,:]+a_PDE_supp2[1,M,:]+a_PDE_supp2[N,1,:]+a_PDE_supp2[N,M,:])/4)*result_scale/pop_scale;
NT_PDE_supp2 = N_PDE_supp2 + A_PDE_supp2;
NI_PDE_supp2 = I_PDE_supp2*result_scale;



# plot parameters
spat_plot_x = 300; spat_plot_y = spat_plot_x; # default size for spatial plots (x and y sizes)
pop_plot_x = 600; pop_plot_y = 300; # default size for population plots (x and y sizes)

supp1_linestyle = :dash; # linestyle for supplementation 1 results
supp2_linestyle = :dashdot; # linestyle for supplementation 2 results

linewidths = 3; # line widths for population plots

col_IL2_supp = RGB(35/255, 148/255, 71/255); # dark green for IL-2 supplements


# plots
plt_IL2 = plot(legend=false, size=(pop_plot_x,pop_plot_y), xlabel="𝑡 (days)", ylabel="Total mass of IL-2 (ng)");
plot!(plt_IL2, (0:T)/24/60, NI_PDE_supp1, linestyle=supp1_linestyle, c=col_IL2_supp, linewidth=linewidths)
plot!(plt_IL2, (0:T)/24/60, NI_PDE_supp2, linestyle=supp2_linestyle, c=col_IL2_supp, linewidth=linewidths)
plot!(plt_IL2, (0:T)/24/60, NI_PDE_sec, c=col_IL2, linewidth=linewidths)

display(plt_IL2)

plt_T = plot(legend=false, size=(pop_plot_x,pop_plot_y), xlabel="𝑡 (days)", ylabel="Number of T cells (millions)");
plot!(plt_T, (0:T)/24/60, NT_PDE_sec, c=col_total, linewidth=linewidths)
plot!(plt_T, (0:T)/24/60, NT_PDE_supp1, linestyle=supp1_linestyle, c=col_total, linewidth=linewidths)
plot!(plt_T, (0:T)/24/60, NT_PDE_supp2, linestyle=supp2_linestyle, c=col_total, linewidth=linewidths)

display(plt_T)




## simulations with no IL-2 decay
T = 2*24*60;

K_λ = 0; # scaling for IL-2 decay from activated T cell concentration (1/Tcell/min),   ~2 for proportional to a*c,  ~0.00002 for proportional to R_na(c)a
λ(a,I) = K_λ*a*I; # activated T cell- and cytokine-dependent or constant decay rate for IL-2  (1/min)    -    2/(60*24)*I from "Mathematical Models of Tumor Immune System Dynamics" (pg 38),  or K_λ*a*I (1/min?) to decay IL-2 proportionally to local activated T cells,  or K_λ*R_na(I)*a (1/min?) to represent consumption by activated T cells when they proliferate (THIS IS NON-LINEAR),  or 0 based off fitted data from (Cheung et al. 2018)

# run PDE model 
n_PDE, a_PDE, I_PDE, runtime_PDE = simulate_PDE(); # outputs from PDE model: spatial concentrations of naive and activated cells, total mass/spatial concentration of IL-2, computation time
T_PDE = n_PDE + a_PDE; # spatial concentrations of total T cells from PDE model

# run ODE model 
n_ODE, a_ODE, I_ODE, runtime_ODE = simulate_ODE(); # outputs from ODE model: total number of naive and activated T cells and total mass of IL-2 (scaled up for the whole dish), computation time
T_ODE = n_ODE + a_ODE; # mass of total T cells from ODE model


# reset to default
K_λ = 2; # scaling for IL-2 decay from activated T cell concentration (1/Tcell/min),   ~2 for proportional to a*c,  ~0.00002 for proportional to R_na(c)a
λ(a,I) = K_λ*a*I; # activated T cell- and cytokine-dependent or constant decay rate for IL-2  (1/min)    -    2/(60*24)*I from "Mathematical Models of Tumor Immune System Dynamics" (pg 38),  or K_λ*a*I (1/min?) to decay IL-2 proportionally to local activated T cells,  or K_λ*R_na(I)*a (1/min?) to represent consumption by activated T cells when they proliferate (THIS IS NON-LINEAR),  or 0 based off fitted data from (Cheung et al. 2018)


pop_scaling = 1e6; # factor this scale out of the total T cell population

linewidths = 3; # population plot linewidths

# mass arrays (with standard deviation for ABM results)
N_ABM = mean(n_tot_sims*result_scale/pop_scaling, dims=2)[:]; SD_N_ABM = std(n_tot_sims*result_scale/pop_scaling, dims=2)[:];
A_ABM = mean(a_tot_sims*result_scale/pop_scaling, dims=2)[:]; SD_A_ABM = std(a_tot_sims*result_scale/pop_scaling, dims=2)[:];
NT_ABM = N_ABM + A_ABM; SD_NT_ABM = std((n_tot_sims+a_tot_sims)*result_scale/pop_scaling, dims=2)[:];
NI_ABM = mean(I_tot_sims*result_scale, dims=2)[:]; SD_NI_ABM = std(I_tot_sims*result_scale, dims=2)[:];

N_PDE = (sum(n_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(n_PDE[1,2:M-1,:],dims=1)+sum(n_PDE[N,2:M-1,:],dims=1)+sum(n_PDE[2:N-1,1,:],dims=1)+sum(n_PDE[2:N-1,M,:],dims=1))[:]/2+(n_PDE[1,1,:]+n_PDE[1,M,:]+n_PDE[N,1,:]+n_PDE[N,M,:])/4)*result_scale/pop_scaling;
A_PDE = (sum(a_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(a_PDE[1,2:M-1,:],dims=1)+sum(a_PDE[N,2:M-1,:],dims=1)+sum(a_PDE[2:N-1,1,:],dims=1)+sum(a_PDE[2:N-1,M,:],dims=1))[:]/2+(a_PDE[1,1,:]+a_PDE[1,M,:]+a_PDE[N,1,:]+a_PDE[N,M,:])/4)*result_scale/pop_scaling;
NT_PDE = N_PDE + A_PDE;
NI_PDE = I_PDE*result_scale;

N_ODE = n_ODE/pop_scaling;
A_ODE = a_ODE/pop_scaling;
NT_ODE = N_ODE + A_ODE;
NI_ODE = I_ODE;


pop_plot_x = 550; pop_plot_y = 300; # default size for population plots (x and y sizes)

plt_pop = plot(legend=false, ylims=(0,NT_ODE[end]), size=(pop_plot_x,pop_plot_y), xlabel="𝑡 (days)", ylabel="Number of T cells (millions)");
plt_IL2 = plot(legend=false, ylims=(0,NI_ODE[end]), size=(pop_plot_x,pop_plot_y), xlabel="𝑡 (days)", ylabel="Total mass of IL-2 (ng)");

plot!(plt_pop, (0:Int(floor(T)))/60/24, N_PDE, c=col_naive, linewidth=linewidths)
plot!(plt_pop, (0:Int(floor(T)))/60/24, N_ODE, c=col_naive, linestyle=:dash, linewidth=linewidths)

plot!(plt_pop, (0:Int(floor(T)))/60/24, A_PDE, c=col_active, linewidth=linewidths)
plot!(plt_pop, (0:Int(floor(T)))/60/24, A_ODE, c=col_active, linestyle=:dash, linewidth=linewidths)

plot!(plt_pop, (0:Int(floor(T)))/60/24, NT_PDE, c=col_total, linewidth=linewidths)
plot!(plt_pop, (0:Int(floor(T)))/60/24, NT_ODE, c=col_total, linestyle=:dash, linewidth=linewidths)

plot!(plt_IL2, (0:Int(floor(T)))/60/24, NI_PDE, c=col_IL2, linewidth=linewidths)
plot!(plt_IL2, (0:Int(floor(T)))/60/24, NI_ODE, c=col_IL2, linestyle=:dash, linewidth=linewidths)

display(plt_pop)
display(plt_IL2)






## Fig 8 - impact of scaffold heterogeneity on T cell expansion, using PDE and ODE models
# ====          ===          ===          ==============        =============
# ====  ===============  =======  =====================  ======  ============
# ====  ===============  =======  =====================  ======  ============
# ====       ==========  =======  ===     ==============        =============
# ====  ===============  =======  ======  =============  ======  ============
# ====  ===============  =======  ======  =============  ======  ============
# ====  ===========          ===          ==============        =============

T = Int(floor(0.5*24*60));
Ts = T;

seed_ρ = 1611;
ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, xpos_PDF, ypos_PDF, rot_PDF, blurring, seed_ρ); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

# run ODE model 
n_ODE, a_ODE, I_ODE, runtime_ODE = simulate_ODE(); # outputs from ODE model: total number of naive and activated T cells and total mass of IL-2 (scaled up for the whole dish), computation time
T_ODE = n_ODE + a_ODE; # mass of total T cells from ODE model


seeds = [1611]#; 1:9]; # seed for scaffold generation, averages results over multiple seeds

# arrays that contain probability density functions (for x and y) that govern the x and y positions of the centre of micro-rods
xpos_PDF_arr = [Uniform(1,M), Uniform(1,M), Normal(M/2,M/5), Normal(M/2,M/10), Exponential(M/5), Exponential(M/7), SymTriangularDist(M/2,M/4), Uniform(M/4,3*M/4)];
ypos_PDF_arr = [Uniform(1,N), Normal(1,N/10), Normal(N/2,N/5), Normal(N/2,N/10), Uniform(1,N), Exponential(N/7), SymTriangularDist(N/2,N/3), Normal(N/2,N/10)];
# names of distributions for x and y (for plotting)
x_PDF_names = ["Uniform", "Uniform", "Normal (high std dev)", "Normal (low std dev)", "Exponential", "Exponential", "Triangular", "Uniform"];
y_PDF_names = ["Uniform", "Edge Normal", "Normal (high std dev)", "Normal (low std dev)", "Uniform", "Exponential", "Triangular", "Middle Normal"];

# re-order by increasing RMSE on these seeds
scaf_RMSE = zeros(length(xpos_PDF_arr));
for scaf_n = eachindex(xpos_PDF_arr) # for each blurring scale
    RMSE_sum = 0;
    for seed = seeds
        ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, x->xpos_PDF_arr[scaf_n], x->ypos_PDF_arr[scaf_n], rot_PDF, blurring, seed); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)
        RMSE_sum = RMSE_sum + rmse(ρ, zeros(N,M).+mean(ρ));
    end
    scaf_RMSE[scaf_n] = RMSE_sum/length(seeds); # record the RMSE between the average scaffold and its mean-field approximation
end
scaf_RMSE_copy = copy(scaf_RMSE);

scaf_ind = Int.(zeros(length(scaf_RMSE)));
for i = eachindex(scaf_RMSE_copy)
    scaf_ind[i] = findall(x->x==minimum(scaf_RMSE_copy), scaf_RMSE_copy)[1];
    scaf_RMSE_copy[scaf_ind[i]] = Inf;
end

xpos_PDF_arr = xpos_PDF_arr[scaf_ind]; # reorder arrays
ypos_PDF_arr = ypos_PDF_arr[scaf_ind];
x_PDF_names = x_PDF_names[scaf_ind];
y_PDF_names = y_PDF_names[scaf_ind];
scaf_RMSE = scaf_RMSE[scaf_ind];



T_scaf = zeros(Ts+1, length(xpos_PDF_arr)); # array containing the total T cells over time for different scaffolds
I_scaf = zeros(Ts+1, length(xpos_PDF_arr)); # array containing the mass of cytokine over time for different scaffolds


for scaf_n = eachindex(xpos_PDF_arr) # for each scaffold
    T_sum = zeros(Ts+1);
    I_sum = zeros(Ts+1);
    for seed = seeds
        ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, x->xpos_PDF_arr[scaf_n], x->ypos_PDF_arr[scaf_n], rot_PDF, blurring, seed); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)
        
        N_rods = mean(ρ)*N*M*result_scale; # APC_ms_conc*APC_ms_vol/m_MSR; # total number of rods in the domain (integral of ρ(x,y))

        # run PDE model 
        print("Scaffold $(scaf_n) of $(length(xpos_PDF_arr))")
        n_PDE, a_PDE, I_PDE, runtime_PDE = simulate_PDE(); # outputs from PDE model: spatial concentrations of naive and activated cells, total mass/spatial concentration of IL-2, computation time

        # calculate mass of total T cells and cytokine and store in arrays
        n_mass_PDE = (sum(n_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(n_PDE[1,2:M-1,:],dims=1)+sum(n_PDE[N,2:M-1,:],dims=1)+sum(n_PDE[2:N-1,1,:],dims=1)+sum(n_PDE[2:N-1,M,:],dims=1))[:]/2+(n_PDE[1,1,:]+n_PDE[1,M,:]+n_PDE[N,1,:]+n_PDE[N,M,:])/4)*result_scale; # naive cells from PDE model
        a_mass_PDE = (sum(a_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(a_PDE[1,2:M-1,:],dims=1)+sum(a_PDE[N,2:M-1,:],dims=1)+sum(a_PDE[2:N-1,1,:],dims=1)+sum(a_PDE[2:N-1,M,:],dims=1))[:]/2+(a_PDE[1,1,:]+a_PDE[1,M,:]+a_PDE[N,1,:]+a_PDE[N,M,:])/4)*result_scale; # activated cells from PDE model
        T_sum = T_sum + n_mass_PDE + a_mass_PDE; # total cells from PDE model

        if use_ODE_for_IL2
            I_sum = I_sum + I_PDE*result_scale;
        else
            I_sum = I_sum + (sum(I_PDE[2:N-1,2:M-1,:],dims=[1,2])[:]+(sum(I_PDE[1,2:M-1,:],dims=1)+sum(I_PDE[N,2:M-1,:],dims=1)+sum(I_PDE[2:N-1,1,:],dims=1)+sum(I_PDE[2:N-1,M,:],dims=1))[:]/2+(I_PDE[1,1,:]+I_PDE[1,M,:]+I_PDE[N,1,:]+I_PDE[N,M,:])/4)*result_scale; # cytokine from PDE model
        end
    end
    T_scaf[:,scaf_n] = T_sum/length(seeds); # average total T cells for this scaffold distribution
    I_scaf[:,scaf_n] = I_sum/length(seeds); # average IL-2 for this scaffold distribution
end


pop_scale = 10^6; # value to factor out of the y axis for plotting

T_ODE = T_ODE/pop_scale; # renaming and scaling
NI_ODE = I_ODE; 
T_PDE_result = T_scaf/pop_scale;
NI_PDE_result = I_scaf;
PDF_arr = [xpos_PDF_arr ypos_PDF_arr];


# plot parameters
col_total = RGB(0, 0, 0); # black for total T cells

spat_plot_x = 250; spat_plot_y = spat_plot_x; # default size for spatial plots (x and y sizes)
pop_plot_x = 500; pop_plot_y = 300; # default size for population plots (x and y sizes)

ODE_linestyle = :dashdot; # linestyle for ODE/MFA results

linewidths = 3; # population plot linewidths

scaff_colors = palette(:viridis,size(PDF_arr,1)); # colours from palette

ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, x->PDF_arr[end,1], x->PDF_arr[end,2], rot_PDF, blurring, seeds[1]); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)
ρ = ρ*m_MSR/h^2*10^8; # scale ρ for plotting (values now represent mass density of micro-rods μg/μm^2)
ρ_max_plot = maximum(ρ)/2; # maximum value for plotting

# population plot
plt_pop = plot(legend=false, size=(pop_plot_x,pop_plot_y));
plot!(plt_pop, (0:τ:Int(floor(T)))/60/24, T_ODE, c=col_total, linestyle=ODE_linestyle, linewidth=linewidths)

# fold expansion plot
plt_expansion = plot(legend=false, size=(pop_plot_x,pop_plot_y), seriestype=:scatter);
scatter!(plt_expansion, [0], [maximum(T_ODE)/T_ODE[1]], c=col_total)

scaf_RMSE = zeros(size(PDF_arr,1)); # measure of "non-uniformity": RMSE between each scaffold and a mean-field approximation of the scaffold
pop_RMSE = zeros(size(PDF_arr,1)); # RMSE between population growths using CL and MFA

for scaf = 1:size(PDF_arr,1) # for each scaffold
    # re-create the scaffold 
    ρ, rot_dens, N, M, wid, hei, result_scale, ρ_max, V, ρ_avg, N_rods, ρ_bar, p_d, p_u, p_l, p_r = gen_scaffold(APC_ms_conc, APC_ms_vol, max_d, N_cell_2D, scale, linemethod, x->PDF_arr[scaf,1], x->PDF_arr[scaf,2], rot_PDF, blurring, seeds[1]); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)
    ρ = ρ*m_MSR/h^2*10^8; # scale ρ for plotting (values now represent mass density of micro-rods μg/μm^2)
    
    plt_scaff = heatmap(0:h:hei, 0:h:wid, reverse(ρ,dims=1), xlims=[0,hei], ylims=[0,wid], size=(spat_plot_x,spat_plot_y), colorbar=false, c=cgrad(:grayC, rev=true), clims=(0,ρ_max_plot), aspect_ratio=:equal, x_foreground_color_border=scaff_colors[scaf], y_foreground_color_border=scaff_colors[scaf]); # plot scaffold
    for edge = [0, wid]
        plot!(plt_scaff, [edge,edge], [0,wid], c=scaff_colors[scaf], legend=false, linewidth=7.5)
        plot!(plt_scaff, [0,wid], [edge,edge], c=scaff_colors[scaf], legend=false, linewidth=7.5)
    end
    for tick = [100, 200, 300, 400]
        plot!(plt_scaff, [tick,tick], [0,8], c=:black, legend=false, linewidth=1.5)
        plot!(plt_scaff, [0,8], [tick,tick], c=:black, legend=false, linewidth=1.5)
    end
    display(plt_scaff)
    
    # calculate RMSE between current scaffold and uniform scaffold
    scaf_RMSE[scaf] = rmse(ρ, zeros(N,M).+mean(ρ)); # record the RMSE between this scaffold and its mean-field approximation

    # calculate RMSE between population growth in current scaffold and that of the MFA 
    pop_RMSE[scaf] = rmse(T_ODE, T_PDE_result[:,scaf]);

    # add to population plot
    plot!(plt_pop, (0:τ:Int(floor(T)))/60/24, T_PDE_result[:,scaf], c=scaff_colors[scaf], linewidth=linewidths)

    # add to fold expansion plot 
    scatter!(plt_expansion, [scaf_RMSE[scaf]], [maximum(T_PDE_result[:,scaf])/T_PDE_result[:,scaf][1]], c=scaff_colors[scaf])
end

display(plt_pop)
display(plt_expansion)






## Fig 9 - further model fitting (appendix)
# ====          ===          ===          ==============        =============
# ====  ===============  =======  =====================  ======  ============
# ====  ===============  =======  =====================  ======  ============
# ====       ==========  =======  ===     ==============         ============
# ====  ===============  =======  ======  =====================  ============
# ====  ===============  =======  ======  =====================  ============
# ====  ===========          ===          =============         =============

# data digitized (using WebPlotDigitizer) from "Low interleukin-2 concentration favors generation of early memory T cells over effector phenotypes during chimeric antigen receptor T-cell expansion" (Kaartinen, et al. 2017)
scale_t = 24*60; # scale for time (convert to minutes)
scale_IL2 = 1/18; # scale for IL2 amount (convert to IU to ng) (https://www.stemcell.com/international-units-conversion-data-for-cytokines)

IL2 = [0, 5, 100, 300]*scale_IL2; # IL-2 concentration (ng/mL)

data0_t = [6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 10, 10, 10, 10, 10, 13, 13, 13, 13, 15]*scale_t; # data for 0 IU/mL IL-2  (min)
data0_y = [1.49, 2.23, 2.39, 4.01, 3.02, 7.34, 1.52, 3.51, 2.2, 7.59, 8.11, 3.63, 17.82, 8.82, 17.52, 13.18, 17.82, 7.21, 6.2, 2.6]; #  (fold expansion)
data5_t = [6, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 13, 13, 13, 13, 13, 15, 15, 15, 17, 17]*scale_t; # data for 5 IU/mL IL-2
data5_y = [1.6, 1.11, 2.55, 2.55, 5.16, 5.34, 9.12, 1.4, 1.49, 2.78, 3.69, 12.75, 20.04, 12.75, 11.72, 3.02, 21.07, 6.52, 13.86, 19.05, 14.1, 7.84, 42.57, 5.25, 4.44, 11.15, 2.92, 2.64];
data100_t = [3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 17, 20, 20, 20, 20, 20, 20]*scale_t; # data for 100 IU/mL IL-2
data100_y = [1.11, 1.22, 1.35, 2.39, 2.74, 1.4, 1.89, 2.47, 3.03, 3.35, 4.84, 5.92, 9.16, 21.91, 24.23, 4.23, 5.01, 47.4, 56.06, 1.4, 16.2, 20.83, 36.85, 8.56, 2.65, 25.91, 25.91, 70.9, 143.43, 42.15, 16.75, 86.71, 95.89, 109.67, 22.66, 61.99, 125.42, 134.12, 200.62, 229.43, 354.89, 379.52, 464.16, 548.94, 627.79, 221.86, 300.08, 441.38, 849.11, 596.98, 419.71, 1207.72, 354.89, 742.46, 767.8, 1837.03, 1502.05, 2284.73, 587.04, 878.08, 1073.9, 4545.26];
data300_t = [3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 20, 20, 20, 20, 20, 20, 20]*scale_t; # data for 300 IU/mL IL-2
data300_y = [2.63, 1.73, 1.66, 1.13, 1.13, 18.99, 9.29, 2.63, 1.09, 3.18, 15.39, 5.85, 4.36, 2.13, 40.49, 30.16, 24.44, 20.66, 15.39, 5.97, 5.16, 4.55, 11.23, 1.29, 2.86, 140, 102.13, 67.06, 37.22, 18.99, 16.05, 3.61, 86.31, 59.11, 1122.61, 597.39, 484.1, 236.83, 208.75, 123.4, 52.11, 33.51, 392.29, 331.55, 1273.57, 1032.04, 572.78, 484.1, 872.25, 376.13, 9588.16, 2022.72, 1506.89, 1273.57, 584.96, 64980.85, 6566.77, 4892.13, 4690.65, 948.79, 3212.55, 2022.72];

data = [[data0_t,data0_y], [data5_t,data5_y], [data100_t,data100_y], [data300_t,data300_y]]; # all data collated

T = 20*scale_t; # final time (min)

@. model(t, p) = exp.(p[1]*t); # model to fit to data, with unknown p = [rp] 
p0 = [0.0]; # initial guess for p = [rp]

rp = zeros(length(IL2)); # doubling rate (probability of reproducing) of T cells at corresponding IL-2 concentration
for i in eachindex(rp)
    fit = curve_fit(model, data[i][1], data[i][2], p0); # model fit
    rp[i] = fit.param[1]; # fitted rp value
end

@. rpmodel(C, p) = p[1]*C/(p[2] + C); # model (Michaelis-Menten equation) to fit to concentration (C) vs rp, with unknown p = [rp_max, K_m] where rp_max is the limiting rate, K_m is the Michaelis constant (concentration of substrate at which the reaction rate is half of rp_max)
p0 = [0.00035, 0.002]; # initial guess for p = [rp_max, K_m]
fit = curve_fit(rpmodel, IL2, rp, p0); # model fit
rp_max = fit.param[1]; # fitted rp_max value (limiting rate)
K_m = fit.param[2]; # fitted K_m value (Michaelis constant)

# plot parameters
def_font = 13; # default font size
default(titlefont = (def_font, "times"), legendfont = (def_font, "times"), guidefont = (def_font, "times"), tickfont = (def_font, "times"), framestyle = :box, yminorgrid = true, xminorgrid = true, size = (600,400), linewidth = 2); # default plot settings
col_line = "#4292c6"; # line colour

plot_x = 550; plot_y = 300; # default size for plots (x and y sizes)

rp_scale = 1e-4; # factor this scale out of proliferation rates for plotting

titles = ["0 ng/mL IL-2", "0.28 ng/mL IL-2", "5.55 ng/mL IL-2", "16.66 ng/mL IL-2"];
for i in eachindex(rp)
    fit = curve_fit(model, data[i][1], data[i][2], p0); # model fit
    rp[i] = fit.param[1]; # fitted rp value

    plot(data[i][1]/scale_t, data[i][2], yaxis=:log, seriestype=:scatter, legend=false, c=:black, size=(plot_x,plot_y))
    plot!(xlabel="𝑡 (days)", ylabel="Fold expansion, N_T(t)/N_T(0)", title=titles[i])
    display(plot!((1:T)/scale_t, exp.(rp[i]*(1:T)), yaxis=:log, c=col_line))
end