using HDF5
using Base
using ProgressMeter
# using FFTW
# using Debug
function Loadh5Param(dir)
    h5file = dir * "ParamLLEJulia.h5"
    res = Dict()
    sim = Dict()
    sim_name = ["res", "sim"]
    par = [["Qc", "R", "ng", "Qi", "gamma"], ["dphi","Pin", "Tscan", "dphi_init", "dphi_end", "lbd_pmp", "mu_sim", "mu_fit", "dispfile"]] 
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
tmp_dir = ARGS[1]
print("\n")
print(tmp_dir)
print("\n")
res, sim = Loadh5Param(tmp_dir)
print("Starting Julia\n")
# -- Collect Data ---
# ----------------------------------------------------
c0 = 299792458;
hbar = 6.634e-34/2/pi;
# -- res parameters -- 
ng = res["ng"][1]
R= res["R"][1]
gamm = res["gamma"][1]
L=2*pi*R;

Q0 = res["Qi"][1]
Qc = res["Qc"][1]

# -- sim parameters --
dphi_init = sim["dphi_init"][1]
# dphi_end = sim["dphi_end"]
# lbd = c0./sim["rf"]
t_end = sim["Tscan"][1]
lambda0 = sim["lbd_pmp"][1]
w0 = 2*pi*c0/lambda0
# mu_sim = sim["mu_sim"]
# pump_index = sim["pump_index"]
dphi = sim["dphi"]
Pin = sim["Pin"][1]

mu_sim = sim["mu_sim"]
# print(mu_sim)
mu = collect(mu_sim[1]:mu_sim[2])
# print(mu)
# mu_fit = sim["mu_fit"]   
pmp_sim = find(mu .== 0)[1]
# print(pmp_sim)


j = 1im
# Ns = [mu_sim[1],mu_sim[end]]
tR=L*ng/c0; 
f0=c0/ng/L; 
nlc=-j*L; 
T=1*tR;
theta=w0/Qc*tR;
alpha=1/2* (w0/Q0 + w0/Qc) * tR;

dw = collect(mu_sim[1]:mu_sim[end])*2*pi/T;
Qc_disp=Qc*ones(size(dw));
Qi_disp=Q0*ones(size(dw));
alpha_disp=1/2* ((w0+0*dw)./Qi_disp + (w0+0*dw)./Qc_disp) * tR;
gamma_disp=gamm*(1+0*dw/w0);

theta_disp=(w0+0*dw)./Qc_disp*tR;
dbeta=dphi/L;
dbeta=dbeta-dbeta[pmp_sim];
dphi_Kerr=theta*Pin*gamm*L/alpha^2;

Ein= zeros(size(mu));
Ein[pmp_sim]=sqrt.(Pin)*length(mu)
Ein_couple=sqrt.(theta_disp).*Ein;

tol=1e-3;  
maxiter=6;

dt=0.1/(sqrt(Pin)); 
dphi_init = dphi_init* alpha
delta_phi = pi^2/8*theta*Pin*gamm*L/alpha^2
dphi_end = sim["dphi_end"][1]*alpha + sim["dphi_end"][2]*delta_phi


t_end = t_end*tR
t_ramp=t_end;
Nt= round(t_ramp/tR/dt);
t1=0;

xx = collect(1:Nt);
dphi_pmp = dphi_init+ xx/Nt * (dphi_end- dphi_init);

Ephoton=hbar*(w0+dw);
phase=2*pi*(rand(1,length(mu)));
array=rand(1,length(mu));
Enoise=array'.*sqrt.(Ephoton/2/tR).*exp.(1j*phase') .*length(mu);
ifft_step = plan_ifft(zeros(size(Enoise)))
fft_step = plan_fft(zeros(size(Enoise)))
u0=ifft_step*Enoise; 

S = Dict()
num_probe=1000;
S["u_probe"] = j*zeros(num_probe, length(u0))
S["Em_probe"] = j*zeros(num_probe, length(u0))
S["delta_phi"] = zeros(num_probe,)
S["comb_power"] = zeros(num_probe,)
probe=0;
print(Nt)
pb = Progress(Int(Nt),0.5, "Calculating temporal loop:", Int(25))
for it = 1:1:Nt
    update!(pb, Int(it))
    t1=(dt+t1) 
    delta_phi=  dphi_pmp[Int(it)];

    u0=ifft_step*(fft_step*(u0) + Ein_couple*dt);
    u1=u0;
    cbeta = -alpha_disp -1im*delta_phi +  1im*L*dbeta;
    halfprop = exp.(cbeta*dt/2)
    cnt = 0
    for ii = 1:maxiter
        uhalf = ifft_step*(halfprop.*(fft_step*(u0)));
        half1 = ifft_step*(gamma_disp.*(fft_step*( abs.(u0).^2.*u0) ) )
        half2 = ifft_step*(gamma_disp.*(fft_step*( abs.(u1).^2.*u1) ) ) 
        uv = uhalf .* exp.(j*L.*(half1./u0 + half2./u1)*dt/2);  
        
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
        Em_probe = (fft_step*(u0))/length(mu_sim)
        Em_probe = Em_probe[:,1]
        S["Em_probe"][probe,:]= Em_probe;
        # deleteat!(Em_probe,pmp_sim)
        S["comb_power"][probe]=(sum(abs.(Em_probe).^2))/Pin
        S["delta_phi"][probe] = delta_phi;
   end
   
end
S["Ewg"] = 1im*zeros(size(S["Em_probe"]))
for ii=1:size(S["Em_probe"],1)
    S["Ewg"][ii,:] = Ein/length(mu_sim)-S["Em_probe"][ii,:].*sqrt.(theta_disp)
end
S["freq"] = (w0 + dw)/(2*pi)

SaveResults(tmp_dir, S)
