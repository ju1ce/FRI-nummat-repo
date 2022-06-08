## Razprsene matrike:

Razprsena matrika je implementirana kot dve matrike, V in I, tipa Matrix. Inicializiramo jo s tem, da podamo matriki V in I.

Implementacija getindex je preprosta, samo sledimo enačbi da pridobimo pravilno vrednost iz matrike V.

Pri implementaciji setindex! ločimo 3 primere:

- Vhod je 0. V tem primeru najdemo element na mestu i,j in ga izbrišemo, tako da premaknemo vse desne elemente v vrstici eno mesto levo
- Na mestu i,j je že element. V tem primeru samo zamenjamo element v matriki V.
- Mesto i,j je prazno. V tem primeru vstavimo element na pripadajoče mesto, naslednje elemente v vrstici pa pomaknemo eno mesto desno. Če je vrstica polna, dodamo stolpec celotni matriki.

```
A = RazprsenaMatrika([1 2; 1 0; 1 2], [1 2; 2 0 ; 2 3]) —>

[1 2 0]
[0 1 0]
[0 1 2]

```

Definiramo tudi fukcijo, ki pretvori matriko v razprseno matriko tako, da inicializira prazno matriko in vstavlja vrednosti.

SOR iteracija je implementirana po enačbah. Grafi števila iteracij v odvisnosti od parametra omega so v notebooku hw2.ipynb.

Prvi primer je poenostavljen primer, najden na Wikipedii:

```
mat = RazprsenaMatrika([4 -1 -6 0; -5 -4 10 8; 9 4 -2 0; 1 -7 5 0], [1 2 3 0; 1 2 3 4; 2 3 4 0; 1 3 4 0] )

4×4 RazprsenaMatrika{Int64}:
  4  -1  -6   0
 -5  -4  10   8
  0   9   4  -2
  1   0  -7   5
  
b = [2, 21, -12, -6]

x0 = [0.0, 0, 0, 0]

```

Ta primer najhitreje konvergira pri omegi 0.55.

Drug primer je reševanje vložitve naključnega grafa s 100 vozlišči v 2d, kjer je 5 vozlišč fiksiranih. V tem primeru najhitreje konvergira pri omegi 1, torej ko metoda deluje enako kot gauss-seidel iteracija.
