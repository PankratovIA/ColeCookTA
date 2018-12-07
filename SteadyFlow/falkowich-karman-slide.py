#! /usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import *

K = 1.0
gamma = 1.4
x_s = 0.61
fi_x0 = 0.93
assert fi_x0 > K/(gamma+1)
fi_x1 = 2*K/(gamma+1) - fi_x0

def fi(x):
  if x <= x_s:
    return fi_x0*x
  else:
    return fi_x0*x_s + fi_x1*(x-x_s)

a = array([[0.25, -0.5, 1.0],
           [-1.0,  1.0, 0.0],
           [0.25,  0.5, 1.0]])
b = array([0.0, fi_x0, fi(1.0)])
k = linalg.solve(a, b)

def fi_0(x):
  return (k[0]*(x-0.5) + k[1])*(x-0.5) + k[2]

def sch(f0, f1, xDelta, tDelta, eps):
  #n = len(f0)
  #dx = K + (gamma + 1.0)/(2*xDelta)*(f1[3:] - f1[2:-1] + f1[1:-2] - f1[:-3])
  #dxx1 = f1[3:] - 2*f1[2:-1] + f1[1:-2]
  #dxx0 = f1[2:-1] - 2*f1[1:-2] + f1[:-3]
  #b = zeros(n)
  #b[2:-1] = multiply(dxx1, multiply(dxx0, dx)*tDelta/xDelta + (f1[1:-2] - f0[1:-2])*eps) +\
            #multiply(dxx0, multiply(dxx1, dx)*tDelta/xDelta + (f1[2:-1] - f0[2:-1])*eps)
  #a = zeros((n,n))
  #a[0, 0], a[1, 1], a[n-1, n-1] = 1.0, 1.0, 1.0
  #aii_2 = -dxx1
  #aii_1 = f1[3:]*eps - (2.0*eps + 1.0)*f0[2:-1] + (eps + 2.0)*f0[1:-2] - f0[:-3]
  #aii = f1[3:] - (2.0 - eps)*f0[2:-1] + (1.0 - 2.0*eps)*f0[1:-2] + eps*f0[:-3]
  #aii1 = dxx0
  #for i in range(2, n-1):
    #a[i,i-2] = aii_2[i-2]
    #a[i,i-1] = aii_1[i-2]
    #a[i,i] = aii[i-2]
    #a[i,i+1] = aii1[i-2]
  #return linalg.solve(a, b)
  n = len(f0) - 1
  a, b = zeros((n+1,n+1)), zeros(n+1)
  a[0, 0], a[1, 1], a[n, n] = 1.0, 1.0, 1.0
  dxx1 = f1[2] - 2*f1[1] + f1[0]
  for i in range(2, n):
    dxx0 = dxx1
    dxx1 = f1[i+1] - 2*f1[i] + f1[i-1]
    dx = K - (gamma + 1.0)/(2*xDelta)*(f1[i+1] - f1[i] + f1[i-1] - f1[i-2])
    b[i] = dxx1*(dxx0*dx*tDelta/xDelta + (f1[i-1] - f0[i-1])*eps) +\
           dxx0*(dxx1*dx*tDelta/xDelta + (f1[i] - f0[i])*eps)
    a[i,i-2] = -dxx1
    a[i,i-1] = f1[i+1]*eps - (2.0*eps + 1.0)*f0[i] + (eps + 2.0)*f0[i-1] - f0[i-2]
    a[i,i] = f1[i+1] - (2.0 - eps)*f0[i] + (1.0 - 2.0*eps)*f0[i-1] + eps*f0[i-2]
    a[i,i+1] = dxx0
  return linalg.solve(a, b)

xDelta = 1.0/15.0
xArray = arange(0.0, 1.0+xDelta, xDelta)
tDelta = 0.04*xDelta**2
eps = 0.05

f0 = vectorize(fi_0)(xArray)
f1 = f0.copy()
f2 = f0.copy()

from pyx import *

text.set(mode="tex", lfs="10pt")
text.preamble(r"\input cyracc.def \font\tencyr=wncyr10 \def\cyr{\tencyr\cyracc}")

#a b v g d e \"e zh z i {\u i} k l m n o p r s t u f kh c ch sh shch {\cprime} y {\cdprime}  \`e yu ya

uStyle = [graph.style.symbol(graph.style.symbol.circle, size=0.08,\
    symbolattrs=[deco.filled, deco.stroked])]

g = graph.graphxy(width=8, x=graph.axis.linear(title=r"$x$", min=0.0, max=1.0),\
                          y=graph.axis.linear(title=r"$\varphi(x)$"))
g.plot(graph.data.function("y(x)=fi(x)", context=locals()))
g.plot(graph.data.function("y(x)=fi_0(x)", context=locals()))
data = graph.data.list([(xArray[i], f2[i]) for i in range(len(xArray))], x=1, y=2)
g.plot(data, uStyle)
g.finish()

s = [text.halign.boxright, text.valign.middle]
g.stroke(path.line(-1, -1.4, 0, -1.4), [style.linewidth.normal, style.linestyle.solid])
g.text(0.5, -1.5, r"\cyr tochnoe reshenie ")
g.stroke(path.line(-1, -1.9, 0, -1.9), [style.linewidth.normal, style.linestyle.dashed])
g.text(0.5, -2, r"\cyr nachal{\cprime}noe priblizhenie ")
g.stroke(path.circle(-0.5, -2.4, 0.04), [deco.filled])
g.text(0.5, -2.5, r"\cyr reshenie poluchennoe metodom ustanovleniya")

g.writePDFfile("../ch4-1-5-%i" % 0)

for i in range(10000):
  f2 = f1 + sch(f0, f1, xDelta, tDelta, eps)
  f0 = f1
  f1 = f2
  if i % 500 == 0:
    #print i
    g = graph.graphxy(width=8, x=graph.axis.linear(title=r"$x$", min=0.0, max=1.0),\
                              y=graph.axis.linear(title=r"$\varphi(x)$"))
    g.plot(graph.data.function("y(x)=fi(x)", context=locals()))
    g.plot(graph.data.function("y(x)=fi_0(x)", context=locals()))
    data = graph.data.list([(xArray[j], f2[j]) for j in range(len(xArray))], x=1, y=2)
    g.plot(data, uStyle)
    g.finish()

    s = [text.halign.boxright, text.valign.middle]
    g.stroke(path.line(-1, -1.4, 0, -1.4), [style.linewidth.normal, style.linestyle.solid])
    g.text(0.5, -1.5, r"\cyr tochnoe reshenie ")
    g.stroke(path.line(-1, -1.9, 0, -1.9), [style.linewidth.normal, style.linestyle.dashed])
    g.text(0.5, -2, r"\cyr nachal{\cprime}noe priblizhenie ")
    g.stroke(path.circle(-0.5, -2.4, 0.04), [deco.filled])
    g.text(0.5, -2.5, r"\cyr reshenie poluchennoe metodom ustanovleniya")

    g.writePDFfile("../ch4-1-5-%i" % (i/500))
