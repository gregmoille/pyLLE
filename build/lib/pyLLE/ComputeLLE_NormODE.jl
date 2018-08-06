using HDF5
using Base
using ProgressMeter
using DifferentialEquations
# using FFTW
# using Debug
function Loadh5Param(dir)
    h5file = dir * "ParamLLEJulia.h5"
    sim = Dict()
    par = ["dw_tot", "F2", "dw_ext", "dw_init", "fact_normE","dw_end","Dint","Tscan", "mu_sim", "t_end", "w", "tR"] 
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
        attrs(g)["Description"] = "This gourp. contrain the results"; # an attribute
    end
        # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute

end



print("\n")
print("Starting Julia\n")
tmp_dir = ARGS[1]
sim = Loadh5Param(tmp_dir);
# -- Collect Data ---
# ----------------------------------------------------
c0 = 299792458;
hbar = 1.0545718e-34;


# -- sim parameters --
# ---------------------------------------------------
dw_tot = sim["dw_tot"][1];
dw_ext = sim["dw_ext"][1];
dw_init = sim["dw_init"][1];
dw_end = sim["dw_end"][1];
Dint = sim["Dint"][:,1];
Tscan = sim["Tscan"][1];
t_end = sim["t_end"][1];
mu_sim = sim["mu_sim"][:,1];
w = sim["w"][:,1];
tR = sim["tR"][1];
fact_normE = sim["fact_normE"][1];
F2 = sim["F2"][1];


dbeta =  -2 * Dint/dw_tot;
mu = collect(mu_sim[1]:mu_sim[2]);
pmp_sim = find(mu .== 0)[1];

# -- Noise Input -- 
Ephoton=hbar.*w;
phase=2*pi*(rand(length(mu)));
array=rand(length(mu));
Enoise=array.*sqrt.(Ephoton).*exp.(1im*phase)*length(mu)* fact_normE;

# -- Setup the FFT -- 
ifft_step = plan_ifft(zeros(size(Enoise)));
fft_step = plan_fft(zeros(size(Enoise)));

# Normalize
# ----------------------------------------------------
t_ramp=t_end;
Nt= sim["Tscan"][1];
dt=t_end/Nt; 
t1=0;

# -- Pump -- 
Ein= zeros(size(mu));
Ein[pmp_sim]=sqrt.(F2)*length(mu);
# Ein_couple=Ein





# Simulation Setup
# ----------------------------------------------------
u0=ifft_step*Enoise; 
tol=1e-3;  
maxiter=6;
S = Dict();
num_probe=1000;
S["u_probe"] = 1im*zeros(num_probe, length(u0));
S["Em_probe"] = 1im*zeros(num_probe, length(u0));
S["detuning"] = zeros(num_probe,);
S["comb_power"] = zeros(num_probe,);



xx = collect(1:Nt);


# phase=2*pi*(rand(1,length(mu)));
# array=rand(1,length(mu));



print("\n")
print("Starting LLE Loop\n")

# pb = Progress(Int(floor(t_end)),0.5, "Calculating temporal loop:", Int(25))
prct = 0;
save = 0

function UpdatePB(u, t, integrator)
    str = "";
    test = 100*(t/t_end);
    msg = "Computing LLE ODE:";
    
    bar = "";
 
    # if prct != floor(test)
    prct = floor(test) ;
    for ii = 1:floor(prct/4)
        bar = bar * "|";
    end
    for ii = floor(prct/4)+1:25 
        bar = bar * " ";
    end
    str = msg * bar * string(prct) * "%";
    print(str, "\u1b[1G");
    # end
    

end

function LLEfun(psi, p, t)
    disp = 1im.*(ifft_step*(dbeta.* (fft_step*(psi))))
    NL = 1im*abs.(psi).^2;
    dw_pmp = dw_init+ t/t_end * (dw_end- dw_init);
    losses = 1+1im.*dw_pmp;
    drive = ifft_step*(Ein) ;
    
    # 
    dpsi_t = -losses.* psi + NL.*psi - disp + drive;
    return dpsi_t;
end

tspan = (0.0,t_end)
cb = FunctionCallingCallback(UpdatePB);
prob = ODEProblem(LLEfun,u0,tspan);
sol = solve(prob,Tsit5(),reltol=1e-6,abstol=1e-6,callback=cb);
print("\nSaving data for transfer to Python");

for ii =1:num_probe
    ind = indmin(abs.(sol.t - ii*sol.t[end]/num_probe))
    u0 = sol.u[ind]
    Em_probe = (fft_step*(u0))/length(mu_sim)
    S["u_probe"][ii,:] = u0
    S["Em_probe"][ii,:]= Em_probe;
    S["comb_power"][ii]=(sum(abs.(Em_probe).^2))/F2
    alpha = dw_init+ sol.t[ind]/t_end * (dw_end- dw_init);
    S["detuning"][ii] = alpha;
end





# probe=0;

# for it = 1:1:Nt
#     update!(pb, Int(it))
#     t1=(dt+t1) 
#     alpha=  dw_pmp[Int(it)];
#     u0=ifft_step*(fft_step*(u0) + Ein*dt);
#     u1=u0;
#     cbeta = -1 -1im*alpha -  1im*dbeta;
#     halfprop = exp.(cbeta*dt/2)
#     uhalf = ifft_step*(halfprop.*(fft_step*(u0)));
#     cnt = 0
#     for ii = 1:maxiter
#         half1 = ifft_step*(1.*(fft_step*( abs.(u0).^2.*u0) ) )
#         half2 = ifft_step*(1.*(fft_step*( abs.(u1).^2.*u1) ) ) 
#         uv = uhalf .* exp.(1im.*(half1./u0 + half2./u1)*dt/2);  
#         uv2 = ifft_step*(halfprop.*(fft_step*(uv)));
#         if (norm(uv2-u1,2)/norm(u1,2) < tol)
#             u1 = uv2;
#             break;
#         else
#             u1 = uv2;
#         end
#         cnt += 1
#     end

#     if (cnt == maxiter)
#         print("Failed to converge to ...");
#     end
#     u0=u1;

#    if (it*num_probe/Nt > probe)
#         probe += 1;
#         S["u_probe"][probe,:]=u0[:,1];
#         Em_probe = (fft_step*(u0))/length(mu_sim)
#         Em_probe = Em_probe[:,1]
#         S["Em_probe"][probe,:]= Em_probe;
#         # deleteat!(Em_probe,pmp_sim)
#         S["comb_power"][probe]=(sum(abs.(Em_probe).^2))/F2
#         S["detuning"][probe] = alpha;
#    end
   
# end
S["Ewg"] = 1im*zeros(size(S["Em_probe"]))
for ii=1:size(S["Em_probe"],1)
    S["Ewg"][ii,:] = Ein/length(mu_sim)-S["Em_probe"][ii,:].*sqrt.(dw_ext*tR)
end
# S["freq"] = (w0 + dw)/(2*pi)

SaveResults(tmp_dir, S)
