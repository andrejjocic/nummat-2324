module Vaje05

export potencna

"""
    v, lambda = potencna(A, x0)

Izračunaj lastni vektor `v` za lastno vrednost `lambda` za matriko `A` s potenčno metodo
z začetnim približkom za lastni vektor `x0`.
"""
function potencna(A, x0, maxit=100, tol=1e-10)
  x = x0
  ls, index = findmax(abs, x)
  ls = x[index]
  x = x / ls # normiramo, da je maksimalni element enak 1
  for i = 1:maxit
    x = A * x
    ln, index = findmax(abs, x)
    ln = x[index]
    x = x / ln # x normiramo, da je maksimalni element enak 1
    if abs(ln - ls) < tol
      println("Potenčna metoda se je končala po $i korakih.")
      return x, ln
    end
    ls = ln
  end
  throw("Potenčna metoda ne konvergira v $maxit korakih.")
end

end # module Vaje05
