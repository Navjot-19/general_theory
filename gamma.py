import sympy as sp

t, r, theta, phi = sp.symbols('t r theta phi')
coords = [t, r, theta, phi]

# metric tensor 
g = sp.Matrix([
    [-(1 - 2/r), 0, 0, 0],
    [0, 1/(1 - 2/r), 0, 0],
    [0, 0, r**2, 0],
    [0, 0, 0, r**2 * sp.sin(theta)**2]
])

# Calculate the inverse metric tensor
g_inv = g.inv()

# Function to calculate Christoffel symbols
def christoffel_symbols(metric, coords):
    n = len(coords)
    Gamma = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
    
    for l in range(n):
        for m in range(n):
            for n_ in range(n):
                
                Gamma[l][m][n_] = sp.Rational(1, 2) * sum(
                    g_inv[l, k] * (
                        sp.diff(metric[k, m], coords[n_]) +
                        sp.diff(metric[k, n_], coords[m]) -
                        sp.diff(metric[m, n_], coords[k])
                    ) for k in range(n)
                )
                Gamma[l][m][n_] = sp.simplify(Gamma[l][m][n_])
    
    return Gamma

# Christoffel symbols
Gamma = christoffel_symbols(g, coords)

# Display the Christoffel symbols
for l in range(len(coords)):
    for m in range(len(coords)):
        for n_ in range(len(coords)):
            if Gamma[l][m][n_] != 0:
                print(f"Γ^{coords[l]}_{coords[m]}{coords[n_]} = {Gamma[l][m][n_]}")
    


