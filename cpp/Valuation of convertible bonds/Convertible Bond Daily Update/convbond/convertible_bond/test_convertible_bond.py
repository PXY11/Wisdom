from convertible_bond import convertible_bond_pde as CB
import numpy as np
from time import time
from math import ceil

T = 3.4904
Nt = ceil(252 * 3.4904)
Ns = 500
faceValue = 100.0
conversionPrice = 9.26
S0 = 7.95
Smax = 4 * S0
sigma = 1.3176
finalCouponRate = 0.04
q = 0
p = 0
R = 1
eta = 0
QCTrigger = 15.360488000000002
QPTrigger = 10
theta = 1
r = 0.0345
annuFactor = 252.0
oas = 0.1
noCallTime = np.array([0, 5.0])
noPutTime = np.array([0, 5.0])
noConvertTime = np.array([0, 5.0])
couponTimeRate = np.array([0.02, 0.015, 0.01, 0.008, 0.006])
# t0 = time()

print CB(Nt, Ns, faceValue, conversionPrice, Smax, S0, sigma, T, finalCouponRate, q, p, R, eta,
         QCTrigger, QPTrigger, theta, r, oas, annuFactor, noCallTime, noPutTime, noConvertTime,
         couponTimeRate)

# print time()-t0
