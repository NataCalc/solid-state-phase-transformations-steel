#!/usr/bin/env python


from pylab import *


def delta(ux0,uxa,uxg,R):
	omega = (uxg-ux0)/max((uxg-uxa),1e-10)
	res = 2.*R*(1-omega)/omega
	return where(res<0.,1e-15,res)


Dx  = 8.232531e-08 
R   = 1.e-2
v   = 6.58331e-11
uxa = 0.0268187
uxg = 0.0756709


ux0 = arange(0.,0.101,0.01)
f	 = (uxg-uxa)*v*delta(ux0,uxa,uxg,R) - Dx*(uxg-ux0)
l	 = -Dx*(uxg-ux0)

plot(ux0,f,'r-')
plot(ux0,l,'k-')
axhline(0.)
show()


