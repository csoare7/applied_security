import Crypto.Util.number as number
import math

b = 2**32
N = 667
n_N = int( math.ceil( math.log( N, b ) ) )
rho = b ** n_N
rho_2 = pow( rho , 2, N )
omega = ( -number.inverse ( N, rho ) ) % b
print b, N, n_N , rho , rho_2 , omega

