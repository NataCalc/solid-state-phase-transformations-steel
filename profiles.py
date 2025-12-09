#!/usr/bin/env python

import sys
from numpy import *
from scipy.special import erfc
from pylab import *


def theoretical(x,cgi,cg0,D,t,r):
	return cg0 + (cgi-cg0)*erfc(x/2./sqrt(D*t))/erfc(r/2./sqrt(D*t))


def zener(x,cgi,cg0,delta,pos):
	return cg0 + (cgi-cg0)*exp(-(x-pos)/delta)



if __name__ == '__main__':

	ca  = 0.01
	cgi = 0.02
	cg0 = 0.015

	r      = 1e-8 	#m
	deltai = 1e-9	#m
	delta  = 10e-9	#m
	t 	    = 10.	#s
	D 	    = 1e-15	#m^2/s

	xstep = 1e-9	#m
	xmax  = 1e-7	#m
	xg = arange(r+deltai,xmax,xstep)
	#cg = theoretical(xg,cgi,cg0,D,t,r)
	cg = zener(xg,cgi,cg0,delta,r+deltai)

	x  = concatenate(([0.,r,r+deltai],xg))
	c  = concatenate(([ca,ca,cgi],cg))
	
	plot(x,c,'k-')
	xslope = arange(r+deltai,r+deltai+delta,xstep)
	plot(xslope,cgi-(cgi-cg0)/delta*(xslope-r-deltai),'r-')
	xbase = arange(r,xmax,xstep)
	plot(xbase,ones(len(xbase))*cg0,'k:')

	axis([0.,xmax,0.009,0.021])

	show()


	
