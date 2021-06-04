%matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
import numpy as n
import casadi as ca
x = ca.MX.sym('x',2);
p = ca.MX.sym('p');

z = 1-x[1]**2;
rhs = ca.vertcat(z*x[0] - x[1]+2*ca.tanh(p),x[0])
ode = {'x':x, 'p':p, 'ode':rhs}
F = ca.integrator('F','cvodes', ode, {'tf':1})
u = ca.MX.sym('u',4,1)

x = [0,1]
for k in range(4):
    res = F(x0=x,p=u[k])
    x = res["xf"]
nlp = {'x':u, 'f':ca.dot(u,u), 'g':x}
solver = ca.nlpsol('solver', 'ipopt', nlp)

res = solver(x0=0.2, lbg=0, ubg=0)


x = ca.SX.sym("x", 3)
u_s = ca.SX.sym("u", 3)
rhs = ca.SX.zeros(3,3)

a1 = ca.horzcat(ca.cos(x[2]),-ca.sin(x[2]),0)
a2 = ca.horzcat(ca.sin(x[2]),ca.cos(x[2]),0)
a3 = ca.horzcat(0,0,1)
rhs = ca.vertcat(a1,a2,a3) @ u_s

ode = {'x':x,'p':u_s,'ode':rhs}
u = ca.DM([1,0,0.1])
x_k = ca.DM([0,0,0])
tf = 100
tc = ca.DM(n.linspace(0,tf,101))
mtc = n.linspace(0,tf,101)
model = ca.integrator('vessel', 'cvodes', ode, {'tf':tf, 'grid':mtc} )

out = model(x0=x_k,p=u)['xf'][:,:]

print(x_k.shape)
print(out.shape)
whole = ca.horzcat(x_k,out)
#print(n.array(whole[0,:]))

%matplotlib inline
fig, axs = plt.subplots(1,1)
axs.plot(n.array(whole[1,:])[0,:],
         n.array(whole[0,:])[0,:],
         linewidth=2, label="random diagonal", color='g');
axs.axis('equal')
plt.xlabel('East [m]')
plt.ylabel('North [m]')
fig.show()
