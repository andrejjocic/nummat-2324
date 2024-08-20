# Gauss-Legendrove kvadrature

Andrej Jočić

## Naloga

Izpelji Gauss-Legendreovo integracijsko pravilo na dveh točkah

$$\int_{0}^{h} f(x) \, dx = A f(x_1) + B f(x_2) + R_f$$

vključno s formulo za napako $R_f$. Izpelji sestavljeno pravilo za $\int_{a}^{b} f(x) \, dx$ in napiši program, ki to pravilo uporabi za približno računanje integrala. Oceni, koliko izračunov funkcijske vrednosti je potrebnih, za izračun približka za

$$\int_{0}^{5} \frac{\sin(x)}{x} \, dx$$

na 10 decimalk natančno.

## Uporaba kode

Funkcija `v = integrate_GLQ1(f, a, b, n)` izračuna približek integrala funkcije `f` na intervalu `[a, b]`. Interval razdeli na `n` podintervalov in na vsakem uporabi Gauss-Legendreovo kvadraturo z dvema točkama.

## Zagon testov

V načinu `pkg` napišemo `activate GaussLegendre` in nato `test`.

## Ustvarjanje poročila

Poročilo `scripts/porocilo.ipynb` je napisano v obliki Jupyter notebooka, ki ga enostavno poženemo da se izvedejo vsi izračuni.