module GaussLegendre

export integrate_GLQ1

"""
    v = integrate_GLQ1(f, a, b, n)

Approximate the integral of `f` from `a` to `b` by dividing [a, b] into `n` sub-intervals
and approximating each with the Gauss-Legendre quadrature rule on 2 points.
The function `f` is evaluated a total of `2n` times.
"""
function integrate_GLQ1(f, a, b, n)
    x = range(a, b, length=n+1)
    h2 = (x[2] - x[1]) / 2 # (b - a) / 2
    node_offset = h2 / sqrt(3)

    v = 0.0
    for i in 1:n
        mean = (x[i] + x[i+1]) / 2
        v += f(mean - node_offset) + f(mean + node_offset)
    end
    return h2 * v
end

end # module GaussLegendre
