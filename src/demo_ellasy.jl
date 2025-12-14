using .Ellasy

N = 10
# lad = SingleLadder(N,:RL_con_C, R=1, C=2)
# lad = SingleLadder(N,:RL_con_C, ones(N), 2ones(N), 2ones(N))# 1.0*collect(1:N))
# lad = SingleLadder(N,:LC_con_R, R=1, C=2, L=3)
lad = SingleLadder(N,:RC_con_L, R=1, C=2, L=3)

f_os = buildFirstOrderSystem(lad)

A = f_os[1]
B = f_os[2]
Cout = f_os[3]

using ControlSystems, Plots

sys = ss(A,B,Cout,[0])
using LinearAlgebra
isstable(sys)

ev = eigvals(A)
scatter(real.(ev), imag.(ev))

Tf = 1000.0;
y, t= step(sys, Tf)
plot(step(sys, Tf))

# plot(t[2:end],y[2:end],xaxis=:log)
# plot(t[2:10000],y[2:10000])
