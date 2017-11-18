import math

def BS_Model(S, r, vol, T, X):      #S is current stock price, r is annualized risk-free interest rate, vol is volatility, T is annualized time to expiry, which means T(days) / day count convention (360), X is strike price
    d1 = (math.log(float(S)/X, math.e) + (r + vol ** 2 /2) * T) / vol / math.sqrt(T)
    d2 = d1 - vol * math.sqrt(T)
    C = S * CDF(d1) - X * math.exp(-r*T) * CDF(d2)
    P = X * math.exp(-r*T) * CDF(-d2) - S * CDF(-d1)
    return (C,P)

def CDF(x):
    #'Cumulative distribution function for the standard normal distribution'
    return (1.0 + math.erf(x / math.sqrt(2.0))) / 2.0



c,p = BS_Model(70, 0.12, 0.2, 90.0/360, 80.0)

print ("Call price ", c)
print ("Put price ", p)

