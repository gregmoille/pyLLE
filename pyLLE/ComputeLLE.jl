using HDF5
using Base
using FFTW
using LinearAlgebra

function Loadh5Param(dir)
    h5file = dir * "ParamLLEJulia.h5"
    res = Dict()
    sim = Dict()
    sim_name = ["res", "sim"]
    par = [["Qc", "R", "ng", "Qi", "gamma","dispfile"], ["Pin", "Tscan", "domega_init", "domega_end", "domega", "f_pmp", "mu_sim_center", "ind_pmp", "Dint", "ind_pump_sweep","f_center", "phi_pmp", "D1", "DKSinit_imag","DKSinit_real" ]]
    cnt = 1
    for sim_par = [res, sim]
        for it = par[cnt]
            sim_par[it] = h5read(h5file, sim_name[cnt] *"/" * it)
        end
        cnt += 1
    end

    return res, sim
end


function SaveResultsToFile(dir, S)
    h5file = dir * "ResultsJulia.h5"
    print(h5file)
    print("\n")
    h5open(h5file, "w") do file
        g = create_group(file, "Results") # create a group
        for ii in S
            print(ii[1])
            print("\n")
            g[ii[1]*"Real"] = real(ii[2])              # create a scalar dataset inside the group
            g[ii[1]*"Imag"] = imag(ii[2])
        end
        # attrs(g)["Description"] = "This group contains only a single dataset"
    end
end

# ----------------------------------------------------
# -- STARTING MAIN SCRIPT --
# ----------------------------------------------------
tmp_dir = ARGS[1] # "/var/folders/_v/2zfybdtj589c1fwk78q_q5qr0000gn/T/tmpqmatmei5" #
tol = parse(Float64,ARGS[2]) #1e-3 #
maxiter = parse(Int,ARGS[3]) #10 #
stepFactor =  parse(Float64,ARGS[4]) #1#

# tmp_dir = "/tmp/tmpbyi6lbwx"
res, sim = Loadh5Param(tmp_dir)
logfile =  open(tmp_dir * "Log.log" ,"w")

# ----------------------------------------------------
# -- Collect Data ---
# ----------------------------------------------------
c0 = 299792458
ħ = 6.634e-34/2/pi

# -- res parameters --
ng = res["ng"][1]
R = res["R"][1]
γ = res["gamma"][1]
L = 2*pi*R

Q0 = res["Qi"][1]
Qc = res["Qc"][1]
fpmp = sim["f_pmp"]
Ppmp = sim["Pin"]
φpmp = sim["phi_pmp"]
fcenter= sim["f_center"][1]

DKSinit_real = sim["DKSinit_real"]
DKSinit_imag = sim["DKSinit_imag"]

DKSinit = DKSinit_real .+ 1im .* DKSinit_imag

D1= sim["D1"][1]
FSR = D1/2π
ω0 = 2π*fpmp
ωcenter = 2π*fcenter

# -- Setup the different Parameters --
# ----------------------------------------------------
tR = 1/FSR
T = 1*tR
κext = ω0[1]/Qc*tR
κ0 = ω0[1]/Q0*tR
α = κ0 + κext
# -- Retrieve sim parameters --
# ----------------------------------------------------
# debug = Bool(sim["debug"][1])
δω_init = sim["domega_init"][1]
δω_end = sim["domega_end"][1]
ind_sweep = sim["ind_pump_sweep"] .+ 1 #because of 0 and 1 index start difference between Julia nd pyhton
t_end = sim["Tscan"][1]
Dint = sim["Dint"]


δω = sim["domega"]
ind_pmp = [ii+1 for ii in sim["ind_pmp"]] #because of 0 and 1 index start difference between Julia
μ_sim = sim["mu_sim_center"]
μ = collect(μ_sim[1]:μ_sim[2])
μ0 = Int(1 + (length(μ)-1)/2)  # center of the domain

# -- Losses --
# ----------------------------------------------------
dω = collect(μ_sim[1]:μ_sim[end])*2*pi*FSR
dω_pmp1 = collect(μ_sim[1]-ind_pmp[1]-1:μ_sim[end]-ind_pmp[1]-1)*2*pi*FSR
ω1 = ω0[1] .+ dω_pmp1

# Make sure there is no residue in the center of the domain
Dint=  (Dint .- Dint[μ0])
# dφ_Kerr=κext*Pin*γ*L/α^2

Ptot = 0
for ii=1:length(fpmp)
    global Ptot += Ppmp[ii]
end
# dt=0.2/(sqrt(Ptot))
dt = 1
t_end = t_end*tR
t_ramp = t_end

# -- Sim length --
# ----------------------------------------------------
Nt = round(t_ramp/tR/dt)
t1 = 0

# -- Angle through the resonator --
# ----------------------------------------------------
θ = Array(range(0,stop=2π,length=length(μ)))
# -- Time scale ---
# Here the normaizatiton of time is such that
# t=0 when the pump is align with the cold cavity
# δω=0 → t=0-
# hence we need to rescale if the detuning is not symmetric
# ----------------------------------------------------
δωtot = (abs(δω_end) + abs(δω_init))
δωcent_perc = -1*sign(δω_end+δω_init)*(abs(δω_end + δω_init)/2)/δωtot
t_sim = Array(range(-t_ramp/2 + δωcent_perc*t_ramp,
             stop = t_ramp/2 + δωcent_perc*t_ramp,
                   length=Int(Nt)))

# -- Frequency sweep ---
# index = 1 will always be the main pump, i.e. from where the Dint is defined
# δω_init = δω_init * tR
# δω_end = δω_end * tR
# ----------------------------------------------------
xx = collect(1:Nt)

δω_all = [δω[ii]*ones(Int(Nt),1) for ii in collect(1:length(fpmp))]
# define which thing we are actually sweeping
for ii in ind_sweep
    δω_all[ii][:] = δω_init.+ xx/Nt * (δω_end - δω_init)
end

# -- Seting up Simulation --
# ----------------------------------------------------
# -- Allocate FFT to go faster --
ifft_plan = FFTW.plan_ifft(zeros(length(μ), 1))
fft_plan = FFTW.plan_fft(zeros(length(μ), 1))

# -- Define the noise
# ----------------------------------------------------
function Noise()
    Ephoton=ħ.*ω1
    phase=2*pi*(rand(length(μ),1))
    array=rand(length(μ),1)
    Enoise=array.*sqrt.(Ephoton/2).*exp.(1im*phase).*length(μ)
    return FFTW.ifftshift(ifft_plan*(Enoise)) # not that ifftshift matters much here but let's do it right
end

# -- Pumps input --
# Need to detune them from the center of the domain (first index)
# Need to do it rhough FFT as it seems in the temporal domain
# introducing a shift through exp(μj x θ) introduce bunch of noise
# ---------------------------------------------------------------
Ain = [1im*zeros(length(μ),1) for ii in collect(1:length(fpmp))]
Ein = [1im*zeros(length(μ),1) for ii in collect(1:length(fpmp))]
for ii in collect(1:length(fpmp))
    Ein[ii][μ0+ind_pmp[ii]] = sqrt(Ppmp[ii])*length(μ)
    Ain[ii] = ifft_plan*(FFTW.fftshift(Ein[ii])) .* exp(-1im.*φpmp[ii])
end
# -- Initial State --
# ---------------------------------------------------------------
u0 = 1im * zeros(length(μ),1)
u0[:, 1] = DKSinit

# -- Output Dict --
# ---------------------------------------------------------------
S = Dict()
num_probe=5000
S["u_probe"] = 1im*zeros(num_probe, length(u0))
S["driving_force"] = 1im*zeros(num_probe,length(u0))
S["detuning"] = zeros(num_probe,)
S["t_sim"] = zeros(num_probe,)
S["kappa_ext"] = [κext]
S["kappa_0"] = [κ0]
S["alpha"] = [α]

# -- Misc --
# ---------------------------------------------------------------

Dint_shift = FFTW.ifftshift(Dint)
# ---------------------------------------------------------------
# ------------------------------------------
# -- Define the Functions
# ------------------------------------------
# ---------------------------------------------------------------

function ProgressBar_CallBack(it, Nt, S, u0, param)
    # -- Update the Progress Bar function every % --
    if it/Nt >= param["probe_pbar"]
        param["cnt_pbar"] = param["cnt_pbar"] + 1;
        #print(string(Int(round(param["probe_pbar"]*100))) * "\n");
        logfile =  open(tmp_dir * "log.log" ,"a")
        write(logfile,string(Int(round(param["probe_pbar"]*100))) * "\n")
        close(logfile)
        probe_pbar2 = param["probe_pbar"] + 0.01;
        param["probe_pbar"]=probe_pbar2;
    end
    return param
end

function SaveStatus_CallBack(it, Nt, S, u0, param)
    # -- Save the data for a data set of 1000 --
    if (it*num_probe/Nt > param["probe"])
        param["probe"] += 1;
        S["u_probe"][param["probe"],:]=u0[:,1] #* exp(1im*tsim[Int(it)]*Δω_pmp[Int(it)]/tR);
        S["detuning"][param["probe"]] = δω_all[ind_sweep[1]][it];
        S["t_sim"][param["probe"]] = t_sim[it]
        S["driving_force"][param["probe"],:] = Fdrive(it)
      end
    return param
end

function SaveData(S, num_probe, ω0, dω)
    # -- need to compute a full spectral κext for
    # extraction of the in-cavity power --
    logfile =  open(tmp_dir * "log.log" ,"a")
    write(logfile,string(Int(round(100))) * "\n")
    close(logfile)

    SaveResultsToFile(tmp_dir, S)
end

function Fdrive(it)
    # -- needed to declare it here as we will use it
    #  later to save it in the dict -
    # ------------------------------------------------
    Force = 0 .* exp.(1im.*θ)

    for ii=1:length(fpmp)
        # detuning for each pump
        # which is i(δωj + Dint[μ = μj]) - iμjθ with μ0 being the center of the
        # domain
        # The phase shift iμjθ has already been taken into account previously in
        # Ain (due to noise, it needs to be initilized throug FFT)
        # ---------------------------------------------------------------------
        σ = (δω_all[ii][it] .+ Dint[μ0+ind_pmp[ii]-1] .- δω_all[ind_sweep[1]][it] ) .* t_sim[it]
        # σ = (δω_all[ii][it] .+ Dint[μ0+ind_pmp[ii]-1]) .* t_sim[it]
        Force .= Force .- 1im .* Ain[ii] .*exp(1im .*σ)
    end
    return Force #.- 1im*Noise()
end


tol = 1e-2
maxiter = 10
success = false


L½prop = 1im .*zeros(length(μ),1)
L½prop_cpl = 1im .*zeros(length(μ),1)
A_L½prop = 1im .*zeros(length(μ),1)
NL½prop_0 = 1im .*zeros(length(μ),1)
CPL½prop_0 = 1im .*zeros(length(μ),1)
A½prop = 1im .*zeros(length(μ),1)
Aprop = 1im .*zeros(length(μ),1)
NL½prop_1 = 1im .*zeros(length(μ),1)
NLprop = 1im .*zeros(length(μ),1)
Force = 1im .*zeros(length(μ),1)
retL = 1im .*zeros(length(μ),1)
Eout = 1im .*zeros(length(μ),1)
Aout = 1im .*zeros(length(μ),1)
retNL = 1im .*zeros(length(μ),1)
retcpl = 1im .*zeros(length(μ),1)

function SSFM½step(A0, it, param)
    # ----------------------------------------------------------------------------------------
    # Standard split step fourier method
    # ∂A/∂t (t, τ)= [L + NL]A + Fdrive
    #
    # We can first split the equation in two such that:
    # and ∂A_l/∂t = L A_l
    # and ∂A_Nl/∂t = NL A_nl
    # where L and NL are both operators
    #
    # In the case of the LLE, the linear equation has an obivous solution in the frequency
    # domain (time domain more complicated with the derivatife).:
    # ^A^_l(t, ω)= exp(^L^ t)
    #
    # where the ^^ describe the Fourier Transform of the function
    #
    # After a smal propagation δt the nonlinear term can be expressed as
    # A_nl(t +δt , τ) = exp(NL x δt) A
    #
    # the linear term has an analytical solution in the frequency domain :
    # thus we could easily apply this sollution to our problem such that
    # ^A^(t + δt, ω) = exp(^L^ x δt)^A^_nl(t, ω)
    #
    # From here it is easy to get the valye on A with the driving force, which we will
    # increment for each δt
    #
    # A(t + δt, ω) = FFT-1(^A^(t + δt, ω)) + Fdrive x δt
    #
    #
    # However this direct SSF is not the best stable and implementing a 1/2 SSF is much better 
    # for more information, please checkout G. Agrawal's book (there is litterally the algorithm 
    # in it)


    # ----------------------------------------------------------------------------------------
    # G. Moille - 09/08/2022 - NIST/UMD - gmoille@umd.edu
    # ----------------------------------------------------------------------------------------
    # -- Define the Linear, Non-Linear and drive Force ---
    # Purpose is for easy modification to ad NL or other terms
    function FFT_Lin(it) 
        return -α/2 .+ 1im .* (Dint_shift .- δω_all[ind_sweep[1]][it] .* 2 )*tR
    end

    function NL(uu, it)
        return -1im*( γ*L* abs.(uu).^2 )
    end

    # --- Linear propagation ---
    
    A0 = A0 .+ Fdrive(Int(it)) .*sqrt(κext) .* dt;        
    
    L½prop .= exp.(FFT_Lin(Int(it)) .* dt/2);
    A_L½prop .= ifft_plan* (fft_plan*(A0) .* L½prop);
    NL½prop_0 .= NL(A0, it);
    A½prop .= A0;
    Aprop .= 1im .* zeros(size(A0));
    
    # --- iterative procedur to find ^LN^(t + δt, ω) ---
    success = false;
    for ii in collect(1:maxiter)
        err = 0.0
        success = false
        NL½prop_1 .= NL(A½prop, it);
        NLprop .= (NL½prop_0 .+ NL½prop_1) .* dt/2;
        Aprop .= ifft_plan*( fft_plan*(A_L½prop .* exp.(NLprop) ) .* L½prop )

        # --- check if we concerge or not ---
        err = LinearAlgebra.norm(Aprop-A½prop,2)/LinearAlgebra.norm(A½prop,2)
        if (err < tol)
            success = true
            break
        else
            success = false
            A½prop .= Aprop
        end
    end

    if success
        u0 = Aprop
        return u0
    else
        print("Convergence Error")
    end
end

function MainSolver(Nt, S, u0)
    # Here we are solvint the most general form of the χ(3)-only-LLE
    # For full explenation on how to derive this equation (usually from CMT)
    # please refere Taheri, H. et al  (2017). "Optical lattice trap for Kerr
    # solitons" The European Physical Journal D, 71(6)
    #
    # ∂A/∂t = (-α/2 - iγL |A|² + i TF⁻¹[Dint])A +
    #         + Σ2 √κextⱼ Fⱼ exp(i σⱼt+ i Dint(η=ηⱼ) - i(ηⱼ-η₀)θ)
    #
    # where Dint is computer at η0 (i.e. Dint(η0) = 0)
    # σⱼ is the detuning between the pump and the closest mode ηⱼ freqeuncy
    # and Fⱼ = -i√(Pinⱼ/ħωⱼ) exp(-φⱼ)
    # ----------------------------------------------------------------------------------------
    # G. Moille - 04/23/2020 - NIST/UMD - gmoille@umd.edu
    # last modification : 05/05/2020 -- Cinquo de Mayo!!
    # ----------------------------------------------------------------------------------------
    param = Dict()

    param["tol"] = 1e-3;
    param["maxiter"] = 6;
    param["step_factor"] = 0.1;
    param["probe_pbar"] = parse(Float64,"0.01");
    param["cnt_pbar"] = 0;
    param["probe"] = 0;

    loops = collect(1:1:Nt)

    for it in loops
        # -- Solve the Split Step Fourier Step --
        # u0 = SSFM(u0, it, param)
        u0 = SSFM½step(u0, it, param)
        
        # -- Update the Progress bar --
        param = ProgressBar_CallBack(Int(it), Nt, S, u0, param)
        # -- Save the Data in the dict --
        param = SaveStatus_CallBack(Int(it), Nt, S, u0, param)
    end
    SaveData(S, num_probe, ω0, dω)

end

# -- Start the solver
# ----------------------------------------------------
logfile =  open(tmp_dir * "log.log" ,"w")
write(logfile,string(Int(round(0))) * "\n")
close(logfile)
MainSolver(Nt, S, u0)
