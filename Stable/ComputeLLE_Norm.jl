using HDF5
using Base
using ProgressMeter
# using FFTW
# using Debug
function Loadh5Param(dir)
    h5file = dir * "ParamLLEJulia.h5"
    sim = Dict()
    par = ["δω_tot", "F2", "δω_ext", "α_init", "α_stop","fact_normE","α_end","Dint","Tscan", "μ_sim", "t_end", "ω", "tR","δω_disp_val"] 
    for it = par
        sim[it] = h5read(h5file,  "sim_norm/" * it)
    end

    return sim

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

print("\n")
print("Starting Julia\n")
tmp_dir = ARGS[1]
sim = Loadh5Param(tmp_dir)
# -- Collect Data ---
# ----------------------------------------------------
c0 = 299792458;
ħ = 1.0545718e-34;


# -- sim parameters --
# ---------------------------------------------------
δω_tot = sim["δω_tot"][1]
δω_ext = sim["δω_ext"][1]
δω_disp = sim["δω_disp_val"][1]
α_init = sim["α_init"][1]
α_end = sim["α_end"][1]
α_stop = sim["α_stop"][1]
Dint = sim["Dint"][:,1]
Tscan = sim["Tscan"][1]
t_end = sim["t_end"][1]
μ_sim = sim["μ_sim"][:,1]
ω = sim["ω"][:,1]
tR = sim["tR"][1]
fact_normE = sim["fact_normE"][1]
F2 = sim["F2"][1]

dβ =  -2 * Dint/(δω_tot)

μ = collect(μ_sim[1]:μ_sim[2])
pmp_sim = find(μ .== 0)[1]

# -- Noise Input -- 
Ephoton=ħ.* ω;
φ = 2*pi*(rand(length(μ)));
array= rand(length(μ));
Enoise=array.*sqrt.(Ephoton).*exp.(1im*φ)*length(μ)* fact_normE

# -- Setup the FFT -- 
ifft_step = plan_ifft(zeros(size(Enoise)))
fft_step = plan_fft(zeros(size(Enoise)))

# Normalize
# ----------------------------------------------------
# dphi_init = 2*dphi_init/δω_tot
# dphi_end = 2*dphi_end/δω_tot

# t_end = t_end*tR*δω_tot/
t_ramp=t_end;
Nt= sim["Tscan"][1];
dt=t_end/Nt; 
t1=0;

# -- Pump -- 
Ein= 1im.*zeros(size(μ));
Ein[pmp_sim]=sqrt.(F2)*length(μ)
# Ein_couple=Ein


# Siμlation Setup
# ----------------------------------------------------
u0 = ifft_step*Enoise; 
tol = 1e-3;  
maxiter = 6;
S = Dict()
num_probe = 1000;
S["u_probe"] = 1im*zeros(num_probe, length(u0))
S["Em_probe"] = 1im*zeros(num_probe, length(u0))
S["detuning"] = zeros(num_probe,)
S["comb_power"] = zeros(num_probe,)



xx = collect(1:Nt);
if typeof(α_stop) == String
    dω_pmp = α_init+ xx/Nt * (α_end- α_init);
else
    dω_pmp = α_init+ xx/Nt * (α_end- α_init);
    ind = indmin(abs.(dω_pmp-α_stop))
    dω_pmp[ind:end] = α_stop
end

# phase=2*pi*(rand(1,length(μ)));
# array=rand(1,length(μ));



print("\n")
print("Starting LLE Loop\n")
probe=0;
pb = Progress(Int(Nt),0.5, "Calculating temporal loop:", Int(25))
for it = 1:1:Nt
    update!(pb, Int(it))
    t1 = (dt+t1) 
    α =  dω_pmp[Int(it)];
    u0 = ifft_step*(fft_step*(u0) + Ein*dt);
    u1 = u0;
    cβ = - (δω_disp +1im*α) + 1im*dβ/2 ;
    halfprop = exp.(cβ.*dt/2)
    uhalf = ifft_step*(halfprop.*(fft_step*(u0)));
    half0 = ifft_step*(1.*(fft_step*( abs.(u0).^2.*u0) ) )
    cnt = 0
    for ii = 1:maxiter
        half1 = ifft_step*(1.*(fft_step*( abs.(u1).^2.*u1) ) ) 
        uv = uhalf .* exp.(1im.*(half0./u0 + half1./u1)*dt/2);  
        uv2 = ifft_step*(halfprop.*(fft_step*(uv)));
        if (norm(uv2-u1,2)/norm(u1,2) < tol)
            u1 = uv2;
            break;
        else
            u1 = uv2;
        end
        cnt += 1
    end

    if (cnt == maxiter)
        print("Failed to converge to ...");
    end
    u0=u1;

    if (it*num_probe/Nt > probe)
        probe += 1;
        S["u_probe"][probe,:]=u0[:,1];
        Em_probe = (fft_step*(u0))/length(μ)
        Em_probe = Em_probe[:,1]
        S["Em_probe"][probe,:]= Em_probe;
        # deleteat!(Em_probe,pmp_sim)
        S["comb_power"][probe]=(sum(abs.(Em_probe).^2))/F2
        S["detuning"][probe] = α;
    end
   
end
S["Ewg"] = 1im*zeros(size(S["Em_probe"]))
for ii=1:size(S["Em_probe"],1)
    S["Ewg"][ii,:] = Ein./length(μ)-(S["Em_probe"][ii,:].*sqrt.(δω_ext/δω_tot.*tR/2))
end
# S["freq"] = (w0 + dw)/(2*pi)

SaveResults(tmp_dir, S)
