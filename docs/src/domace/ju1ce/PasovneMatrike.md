 
## Domača naloga 1: pasovne matrike

### Podatkovni tipi

Podatkovi tip PasovnaMatrika hrani vektor vrednosti na diagonali, vektor vektorjev vrednosti pod diagonalo in vektor vektorjev vrednosti nad diagonalo. ZgornjePasovnaMatrika in SpodnjePasovnaMatrika sta enaka, le da hranita vektor vektorjev vrednosti le nad/pod diagonalo.

PasovnaMatrika([1,1,1],[[2,2],[3]],[4,4])

=

[1 4 0]
[2 1 4]
[3 2 1]

SpodnjePasovnaMatrika([1,1,1],[[2,2],[3]])

=

[1 0 0]
[2 1 0]
[3 2 1]

### Množenje s vektorjem

Pri množenju s vektorjem se diagonala in vsak vektor ob diagonali zmnoži po vrednostih, nato pa se te vektorje sešteje

### Deljenje s leve

Računanje sistema A*x=b. Računa se preko gaussove eliminacije, nato pa s obratnim vstavljanjem. Pri ZgornjePasovnaMatrika se korak gausove eliminacije preskoči, pri SpodnjePasovnaMatrika pa je poenostavljeno obratno vstavljanje.

## lu razcep

Za lu razcep uporabimo gaussovo eliminacijo, pri kateri pa l vrednosti shranjujemo v SpodnjePasovnaMatrika L. Končno matriko pretvorimo v ZgornjePasovnaMatrika U.
