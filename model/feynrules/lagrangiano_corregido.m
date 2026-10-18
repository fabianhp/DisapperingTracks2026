(* Derivada covariante del quintuplete, CORREGIDA.
   En la version original el termino del[vec[m,i],n] estaba dentro de la suma
   sobre jj (5 valores) y aaa (3 valores), de modo que la derivada quedaba
   multiplicada por 15 y el termino cinetico por 225 (lo muestra la salida de
   la celda "Lminimal" del notebook: -225/4 (d_nu V0_mu)^2 ...).
   Aqui la suma actua solo sobre el termino de gauge. *)

DevQ[m_, n_, i_] := del[Vquint[m, i], n] -
   I gw Sum[gen5[[aa, i, jj]] Vquint[m, jj] Wi[n, aa], {jj, 1, 5}, {aa, 1, 3}];

Lminimal := Block[{mu, nu, iii, feynmangaugerules},
   feynmangaugerules = If[Not[FeynmanGauge], {G0 | GP | GPbar -> 0}, {}];
   ExpandIndices[
     -1/4 Sum[HC[DevQ[mu, nu, ii] - DevQ[nu, mu, ii]] (DevQ[mu, nu, ii] - DevQ[nu, mu, ii]), {ii, 1, 5}]
     + 1/2 mv mv anti[Vquint[mu, iii]] Vquint[mu, iii],
     FlavorExpand -> {SU2Q, SU2W}] /. feynmangaugerules];

Lnonminimal := Block[{mu, nu, feynmangaugerules},
   feynmangaugerules = If[Not[FeynmanGauge], {G0 | GP | GPbar -> 0}, {}];
   ExpandIndices[
     kappa*I gw Sum[HC[Vquint[mu, i]] FS[Wi, mu, nu, a] gen5[[a, i, j]] Vquint[nu, j], {i, 1, 5}, {j, 1, 5}, {a, 1, 3}]
     + xiinv*Sum[HC[DevQ[mu, mu, ii]] DevQ[nu, nu, ii], {ii, 1, 5}],
     FlavorExpand -> {SU2Q, SU2W}] /. feynmangaugerules];

(* Comprobaciones a hacer tras cargar el modelo:
   1. Table[gen5[[a]].gen5[[b]] - gen5[[b]].gen5[[a]], ...] debe dar I Eps[a,b,c] gen5[[c]].
   2. La salida de Lminimal debe empezar con -1/4 (d_nu V0_mu)^2, no con -225/4.
   3. En el UFO, el vertice (a, V-, V+) debe tener acoplamiento de modulo ee
      (carga 1) y el (a, V--, V++) de modulo 2 ee, no 3 Sqrt[5]/2 ee y 3 Sqrt[5] ee.
   4. Con esta normalizacion, la estructura de Yang-Mills (la del W) corresponde
      a kappa = -1/2 en este Lnonminimal; kappa = 1 no es el caso "tipo W".
   5. El coeficiente de LVVpi se dedujo con los generadores reescalados y sin
      V_ud: hay que rederivarlo con los generadores corregidos. *)
