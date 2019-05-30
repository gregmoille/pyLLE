using HDF5
using Base
using ProgressMeter
using FFTW
using LinearAlgebra

function Loadh5Param(dir)
    h5file = dir * "ParamLLEJulia.h5"
    res = Dict()
    sim = Dict()
    sim_name = ["res", "sim"]
    par = [["Qc", "R", "ng", "Qi", "gamma","dispfile"], ["dphi","Pin", "Tscan", "domega_init", "domega_end", "f_pmp", "mu_sim", "debug", "ind_aux"]]
    cnt = 1
    for sim_par = [res, sim]
        for it = par[cnt]
            sim_par[it] = h5read(h5file, sim_name[cnt] *"/" * it)
        end
        cnt += 1
    end

    return res, sim
end


function SaveResults(dir, S)
    h5file = dir * "ResultsJulia.h5"
    # print(h5file)
    # print("\n")
    h5open(h5file, "w") do file
        g = g_create(file, "Results") # create a group
        for ii in S
            g[ii[1]*"Real"] = real(ii[2])              # create a scalar dataset inside the group
            g[ii[1]*"Imag"] = imag(ii[2])
        end
        attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
    end
        # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

# ----------------------------------------------------
# -- STARTING MAIN SCRIPT --
# ----------------------------------------------------
# FFTW.set_num_threads(Sys.CPU_CORES)

tmp_dir = ARGS[1]
tol = parse(Float64,ARGS[2])
#maxiter = Int(float(ARGS[3]))
maxiter = parse(Int,ARGS[3])
#stepFactor = float(ARGS[4])
stepFactor = parse(Float64,ARGS[4])
res, sim = Loadh5Param(tmp_dir)
logfile =  open(tmp_dir * "Log.log" ,"w")

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

# -- sim parameters --
debug = Bool(sim["debug"][1])
δω_init = sim["domega_init"][1]
δω_end = sim["domega_end"][1]
t_end = sim["Tscan"][1]
fpmp = sim["f_pmp"][1]
Pin = sim["Pin"][1]
if length(sim["f_pmp"])>1
    faux = sim["f_pmp"][2:length(sim["f_pmp"])]
    Paux = sim["Pin"][2:length(sim["f_pmp"])]
else
    faux = []
    Paux = []
end
ω0 = 2*pi*fpmp
dφ = sim["dphi"]
μ_sim = sim["mu_sim"]
μ = collect(μ_sim[1]:μ_sim[2])
pmp_sim_nm = findall(μ .== 0)
pmp_sim = pmp_sim_nm[1]
ind_aux = sim["ind_aux"]
write(logfile,string(ind_aux))
print(string(ind_aux))
# -- Setup the different Parameters --
# ----------------------------------------------------
tR = L*ng/c0
f0 = c0/ng/L
T = 1*tR
θ = ω0/Qc*tR
α= 1/2* (ω0/Q0 + ω0/Qc) * tR

<<<<<<< HEAD
# -- losses --
=======
>>>>>>> 93230dac4268a8762b790cb960096c396711e76c
dω = collect(μ_sim[1]:μ_sim[end])*2*pi/T
Qc_disp=Qc*ones(size(dω))
Qi_disp=Q0*ones(size(dω))

# -- Dispersion --
dβ=dφ/L
dβ=dβ.-dβ[pmp_sim]
dφ_Kerr=θ*Pin*γ*L/α^2

# -- Input Energy --
Ein= zeros(size(μ))
Ein[pmp_sim]=sqrt.(Pin)*length(μ)


for ii=1:length(faux)
    Ein[pmp_sim - ind_aux[ii]] = sqrt.(Paux[ii])*length(μ)
end
Ein_couple=sqrt(θ).*Ein

<<<<<<< HEAD
# -- Noise Background --
Ephoton=ħ*(ω0+dω)
=======
# -- Noise Background -- 
Ephoton=ħ*(ω0.+dω)
>>>>>>> 93230dac4268a8762b790cb960096c396711e76c
phase=2*pi*(rand(1,length(μ)))
array=rand(1,length(μ))
Enoise=array'.*sqrt.(Ephoton/2).*exp.(1im*phase') .*length(μ)


dt=0.1/(sqrt(Pin))
δω_init = δω_init * tR
δω_end = δω_end * tR

t_end = t_end*tR
t_ramp=t_end

# -- Seting up Simulation --
# ----------------------------------------------------

<<<<<<< HEAD
# -- Allocate FFT to go faster --
ifft_plan = plan_ifft(zeros(size(Enoise)))
fft_plan = plan_fft(zeros(size(Enoise)))
=======
# -- Allocate FFT to go faster -- 
ifft_plan = FFTW.plan_ifft(zeros(size(Enoise)))
fft_plan = FFTW.plan_fft(zeros(size(Enoise)))
>>>>>>> 93230dac4268a8762b790cb960096c396711e76c


# -- Sim length --
tol=1e-3
# tol=1e-20
maxiter=6
Nt= round(t_ramp/tR/dt)
t1=0

# -- Detuning ramp --
xx = collect(1:Nt)
Δω_pmp = δω_init.+ xx/Nt * (δω_end - δω_init)

# -- Initiial State --
u0=ifft_plan*Enoise

# -- Output Dict --
S = Dict()
num_probe=1000
S["u_probe"] = 1im*zeros(num_probe, length(u0))
S["Em_probe"] = 1im*zeros(num_probe, length(u0))
S["comb_power"] = zeros(num_probe,)
S["detuning"] = zeros(num_probe,)

# -- Progress Bar --
probe=0
# pb = Progress(Int(Nt),0.5, "Calculating temporal loop:", Int(25))

space = "."^100
logfile =  open(tmp_dir * "log.log" ,"a")
write(logfile,string(Int(round(0))) * "\n")
close(logfile)
probe_pbar = parse(Float64,"0.01")
cnt_pbar = 0
# -- Main Solver --
# ----------------------------------------------------
for it = 1:1:Nt
    global probe_pbar
    global u0
    global probe
    global cnt_pbar
    # print("coucouc\n")
<<<<<<< HEAD
    if it/Nt >=probe_pbar
        cnt_pbar = cnt_pbar + 1
=======
    if it/Nt >= probe_pbar
        cnt_pbar = cnt_pbar + 1 
>>>>>>> 93230dac4268a8762b790cb960096c396711e76c
        # print("\b"^103)
        # print("coucou")
        # print("="^cnt_pbar)
        # print("."^(space-cnt_pbar))
        # print(" ")
        logfile =  open(tmp_dir * "log.log" ,"a")
        write(logfile,string(Int(round(probe_pbar*100))) * "\n")
        close(logfile)
<<<<<<< HEAD
        probe_pbar = probe_pbar + 0.01
    end
=======
        probe_pbar2 = probe_pbar + 0.01
        probe_pbar=probe_pbar2
    end 
>>>>>>> 93230dac4268a8762b790cb960096c396711e76c
    u0=ifft_plan*(fft_plan*(u0) + Ein_couple*dt)
    u1=u0
    cbeta = -α  .+ 1im*Δω_pmp[Int(it)] .+  1im*L*dβ
    halfprop = exp.(cbeta*dt/2)
    cnt = 0
    uhalf = ifft_plan*(halfprop.*(fft_plan*(u0)))
    half1 = ifft_plan*(γ*(fft_plan*( abs.(u0).^2 .*u0) ) )./u0
    for ii = 1:maxiter
<<<<<<< HEAD
        half2 = ifft_plan*(γ*(fft_plan*( abs.(u1).^2.*u1)))./u1
        uv = uhalf .* exp.(1im*L.*(half1 + half2)*dt/2)
        uv2 = ifft_plan*(halfprop.*(fft_plan*(uv))) + Enoise
        if (norm(uv2-u1,2)/norm(u1,2) < tol)
=======
        half2 = ifft_plan*(γ*(fft_plan*( abs.(u1).^2 .*u1)))./u1
        uv = uhalf .* exp.(1im*L.*(half1 + half2)*dt/2)  
        uv2 = ifft_plan*(halfprop.*(fft_plan*(uv)))
        if (LinearAlgebra.norm(uv2-u1,2)/LinearAlgebra.norm(u1,2) < tol)
>>>>>>> 93230dac4268a8762b790cb960096c396711e76c
            u1 = uv2
            break
        else
            u1 = uv2
        end
        cnt += 1
    end

    if (cnt == maxiter)
        logfile =  open(tmp_dir * "log.log" ,"a")
        write(logfile,"Failed to converge to N=" * string(it) * " over Nt=" * string(Nt))
        close(logfile)
    end
    u0=u1

   if (it*num_probe/Nt > probe)
        probe += 1
        S["u_probe"][probe,:]=u0[:,1]
        Em_probe = (fft_plan*(u0))/length(μ)
        Em_probe = Em_probe[:,1]
        S["Em_probe"][probe,:]= Em_probe
        deleteat!(Em_probe,pmp_sim)
        S["comb_power"][probe]=(sum(abs.(Em_probe).^2))/Pin
        S["detuning"][probe] = Δω_pmp[Int(it)]/tR
   end
<<<<<<< HEAD

=======
>>>>>>> 93230dac4268a8762b790cb960096c396711e76c
end
S["Ewg"] = 1im*zeros(size(S["Em_probe"]))
for ii=1:size(S["Em_probe"],1)
    S["Ewg"][ii,:] = Ein/length(μ)-S["Em_probe"][ii,:].*sqrt(θ)
end
S["ω"] = (ω0 .+ dω)

# print("\n")
# print("Congrats, everything has been solved.\n")
# print("Now we are going to save the results in a .h5 file\n")
# print("Everything will be in: ")
logfile =  open(tmp_dir * "log.log" ,"a")
write(logfile,string(Int(round(100))) * "\n")
close(logfile)

SaveResults(tmp_dir, S)
# write(logfile,"Great we have saved everything... Moving back to Python\n")
