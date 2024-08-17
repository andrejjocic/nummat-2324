# Obhod Lune

Andrej Jočić

## Opis naloge

Sondo Appolo pošljite iz Zemljine orbite na tir z vrnitvijo brez potiska (free-return trajectory), ki obkroži Luno in se vrne nazaj v Zemljino orbito. Rešujte sistem diferencialnih enačb, ki ga dobimo v koordinatnem sistemu, v katerem Zemlja in Luna mirujeta (omejen krožni problem treh teles). Naloge ni potrebno reševati na 10 decimalk.

## Uporaba kode

Funkcija `sol = solve_CR3BP(problem, init_pos, init_vel, duration, integrator_params)` vrne objekt tipa `IVP_solution`, ki vsebuje tabelo časov `sol.t` in tabelo stanj `sol.y` ob teh časih. Stanje `sol.y[i]` je vektor `[x, y, z, vx, vy, vz]` položaja in hitrosti v trenutku `sol.t[i]`.
Problem je definiran z objektom tipa `CR3BP`, ki ga ustvarimo z `problem = CR3BP(m1, m2)`, kjer sta `m1` in `m2` masi primarnega in sekundarnega telesa (privzeti vrednosti sta za Zemljo in Luno).

Začetni položaj sonde je podan z 3D vektorjem `init_pos`. Podan naj bo relativno na masno središče sistema, in sicer v brezdimenzijskih koordinatah kjer 1 pomeni razdaljo med primarnim in sekundarnim telesom. Najlažje ga nastavimo relativno na eno od teles, npr. `init_pos = problem.pos2 + [0.1, 0, 0]`.
Začetna hitrost sonde je podana z 3D vektorjem `init_vel`. Pri magnitudi hitrosti pazite, da funkcija uporablja  brezdimenzijske časovne enote, kjer $2\pi$ ustreza eni orbitalni periodi sekundarnega telesa.

Čas trajanja simulacije je podan z `duration`, ki je v prej omenjenih brezdimenzijskih časovnih enotah. Metodo reševanja sistema diferencialnih enačb nastavimo z `integrator_params`, ki je bodisi `ParamsRK4(n)` za metodo Runge-Kutta 4. reda (z `n` fiksnimi koraki) ali `ParamsDOPRI5(ε, h0, η)` za adaptivno Dormand-Prince metodo reda 5(4). Pri slednji je `ε` toleranca napake, `h0` velikost začetnega koraka in `η` faktor, ki "za vsak slučaj" zmanjša na novo izračunan korak.

Funckija `L4, L5 = stable_Lagrange_points(problem)` vrne položaje stabilnih Lagrangeevih točk L4 in L5 v brezdimenzijskih koordinatah.

## Zagon testov

V načinu `pkg` napišemo `activate CRThreeBodyProblem` in nato `test`.

## Ustvarjanje poročila

Poročilo `scripts/report.ipynb` je napisano v obliki Jupyter notebooka, ki ga enostavno poženemo da se izvedejo vsi izračuni in prikažejo grafi.