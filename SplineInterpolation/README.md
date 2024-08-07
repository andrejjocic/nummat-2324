# Naravni zlepek

Andrej Jočić

## Opis naloge

Danih je $n$ interpolacijskih točk $(x_i,f_i)$ , $i=1,2,...,n$. **Naravni interpolacijski kubični zlepek** $S$ je funkcija, ki izpolnjuje naslednje pogoje:

1. $S(x_i)=f_i, i=1,2,...,n.$

2. $S$ je polinom stopnje 3 ali manj na vsakem podintervalu $[x_i,x_{i+1}]$, $i=1,2,...,n−1$.

3. $S$ je dvakrat zvezno odvedljiva funkcija na interpolacijskem intervalu $[x_1,x_n]$

4. $S^{′′}(x_1)=S^{′′}(x_n)=0$.

Zlepek $S$ določimo tako, da postavimo

$S(x)=S_i(x)=a_i+b_i(x−x_i)+c_i(x−x_i)^2+d_i(x−x_i)^3, x∈[x_i,x_{i+1}]$,

nato pa izpolnimo zahtevane pogoje.

## Uporaba kode

Funkcija `S = interpolate(x, y)` vrne objekt tipa `CubicSpline`, ki interpolira točke (`x[i]`, `y[i]`) z naravnim zlepkom. Točke morejo biti urejene po `x` v naraščajočem vrstnem redu.

Funkcija `evaluate(cs::CubicSpline, x)` vrne vrednost zlepka v točki `x`.

Funkcija `plot_spline(cs::CubicSpline, resolution = 50)` nariše graf zlepka na definicijskem območju.
Odseki so izmenično pobarvani z rdečo in modro barvo.
Parameter `resolution` določa koliko točk bo evalviranih za risanje vsakega odseka.

## Zagon testov

V načinu `pkg` napišemo `activate SplineInterpolation` in nato `test`.

## Ustvarjanje poročila

Poročilo `scripts/report.ipynb` je napisano v obliki Jupyter notebooka, ki ga enostavno poženemo da se izvedejo vsi izračuni in prikažejo grafi.