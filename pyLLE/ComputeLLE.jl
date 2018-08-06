using HDF5
using Base
using ProgressMeter


function Loadh5Param(dir)
    h5file = dir * "ParamLLEJulia.h5"
    res = Dict()
    sim = Dict()
    sim_name = ["res", "sim"]
    par = [["Qc", "R", "ng", "Qi", "gamma"], ["dphi","Pin", "Tscan", "domega_init", "domega_end", "f_pmp", "mu_sim",  "dispfile", "debug"]] 
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
res, sim = Loadh5Param(tmp_dir)

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
dω_init = sim["domega_init"][1]
t_end = sim["Tscan"][1]
fpmp = sim["f_pmp"][1]
ω0 = 2*pi*fpmp
dφ = sim["dphi"]
Pin = sim["Pin"][1]
μ_sim = sim["mu_sim"]
μ = collect(μ_sim[1]:μ_sim[2])
pmp_sim = find(μ .== 0)[1]

# -- Setup the different Parameters -- 
# ----------------------------------------------------
tR = L*ng/c0
f0 = c0/ng/L 
T = 1*tR
θ = ω0/Qc*tR 
α=1/2* (ω0/Q0 + ω0/Qc) * tR

# -- losses -- 
dω = collect(μ_sim[1]:μ_sim[end])*2*pi/T
Qc_disp=Qc*ones(size(dω))
Qi_disp=Q0*ones(size(dω))
# α_disp=1/2* ((ω0+0*dω)./Qi_disp + (ω0+0*dω)./Qc_disp) * tR
# γ_disp=gamm*(1+0*dω/ω0)

# -- Dispersion -- 
# θ_disp=(ω0+0*dω)./Qc_disp*tR
dβ=dφ/L
dβ=dβ-dβ[pmp_sim]
dφ_Kerr=θ*Pin*γ*L/α^2

# -- Input Energy -- 
Ein= zeros(size(μ))
Ein[pmp_sim]=sqrt.(Pin)*length(μ)
Ein_couple=sqrt(θ).*Ein

# -- Noise Background -- 
Ephoton=ħ*(ω0+dω)
phase=2*pi*(rand(1,length(μ)))
array=rand(1,length(μ))
Enoise=array'.*sqrt.(Ephoton/2/tR).*exp.(1im*phase') .*length(μ)


dt=0.1/(sqrt(Pin)) 
dω_init = dω_init* α
delta_phi = pi^2/8*θ*Pin*γ*L/α^2
dω_end = sim["domega_end"][1]*α + 1*delta_phi


t_end = t_end*tR
t_ramp=t_end




# -- Seting up Simulation -- 
# ----------------------------------------------------

# -- Allocate FFT to go faster -- 
ifft_plan = plan_ifft(zeros(size(Enoise)))
fft_plan = plan_fft(zeros(size(Enoise)))


# -- Sim length -- 
tol=1e-3  
maxiter=6
Nt= round(t_ramp/tR/dt)
t1=0

# -- Detuning ramp -- 
xx = collect(1:Nt)
dω_pmp = dω_init+ xx/Nt * (dω_end- dω_init)

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
pb = Progress(Int(Nt),0.5, "Calculating temporal loop:", Int(25))


# -- Main Solver -- 
# ----------------------------------------------------
for it = 1:1:Nt
    update!(pb, Int(it))
    u0=ifft_plan*(fft_plan*(u0) + Ein_couple*dt)
    u1=u0
    cbeta = -α -1im*dω_pmp[Int(it)] +  1im*L*dβ
    halfprop = exp.(cbeta*dt/2)
    cnt = 0
    uhalf = ifft_plan*(halfprop.*(fft_plan*(u0)))
    half1 = ifft_plan*(γ*(fft_plan*( abs.(u0).^2.*u0) ) )./u0
    for ii = 1:maxiter
        half2 = ifft_plan*(γ*(fft_plan*( abs.(u1).^2.*u1)))./u1
        uv = uhalf .* exp.(1im*L.*(half1 + half2)*dt/2)  
        uv2 = ifft_plan*(halfprop.*(fft_plan*(uv)))
        if (norm(uv2-u1,2)/norm(u1,2) < tol)
            u1 = uv2
            break
        else
            u1 = uv2
        end
        cnt += 1
    end

    if (cnt == maxiter)
        print("Failed to converge to ...")
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
        S["detuning"][probe] = delta_phi
   end
   
end
S["Ewg"] = 1im*zeros(size(S["Em_probe"]))
for ii=1:size(S["Em_probe"],1)
    S["Ewg"][ii,:] = Ein/length(μ_sim)-S["Em_probe"][ii,:].*sqrt(θ)
end
S["ω"] = (ω0 + dω)

SaveResults(tmp_dir, S)
