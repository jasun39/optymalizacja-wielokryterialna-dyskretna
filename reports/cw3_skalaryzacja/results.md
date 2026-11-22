### Wariant (a) – front wypukły, U = [-1, 1]² [a_convex]

Tablica wypłat (wiersz i = F w minimum kryterium i):

| min | F1 | F2 |
|---|---|---|
| F1 | 0.0000 | 1.9991 |
| F2 | 1.9991 | 0.0000 |

ideal = (0.000, 0.000), nadir = (1.999, 1.999), r = (-0.200, -0.200)

| metoda | zadań | rozwiązanych | niezdominowanych | unikalnych niezd. | czas [s] | najbliższy idealnemu: x | F(x) | d₂ |
|---|---|---|---|---|---|---|---|---|
| S1 | 19 | 19 | 19 | 19 | 0.01 | (0.500, 0.500) | (0.500, 0.500) | 0.3537 |
| S2 | 21 | 21 | 21 | 21 | 0.03 | (0.500, 0.500) | (0.500, 0.500) | 0.3537 |
| S3 (p=2) | 19 | 19 | 19 | 19 | 0.03 | (0.500, 0.500) | (0.500, 0.500) | 0.3537 |
| S3 (p=∞) | 19 | 19 | 19 | 19 | 0.02 | (0.500, 0.500) | (0.500, 0.500) | 0.3537 |

Zgodność z rozwiązaniem analitycznym (P = odcinek x1 = x2 ∈ [0, 1], front √F1 + √F2 = √2), punkty niezdominowane:

- S1: max |√F1 + √F2 − √2| = 1.15e-11, max |x1 − x2| = 1.76e-06
- S2: max |√F1 + √F2 − √2| = 5.70e-13, max |x1 − x2| = 6.85e-07
- S3 (p=2): max |√F1 + √F2 − √2| = 9.04e-13, max |x1 − x2| = 1.13e-06
- S3 (p=∞): max |√F1 + √F2 − √2| = 5.59e-10, max |x1 − x2| = 1.32e-06

### Wariant (a) – front niewypukły (Fonseca–Fleming), U = [-1, 1]² [a_fonseca]

Tablica wypłat (wiersz i = F w minimum kryterium i):

| min | F1 | F2 |
|---|---|---|
| F1 | 0.0000 | 0.9817 |
| F2 | 0.9817 | 0.0000 |

ideal = (0.000, 0.000), nadir = (0.982, 0.982), r = (-0.098, -0.098)

| metoda | zadań | rozwiązanych | niezdominowanych | unikalnych niezd. | czas [s] | najbliższy idealnemu: x | F(x) | d₂ |
|---|---|---|---|---|---|---|---|---|
| S1 | 19 | 19 | 19 | 19 | 0.02 | (0.677, 0.677) | (0.002, 0.978) | 0.9966 |
| S2 | 21 | 21 | 21 | 21 | 0.06 | (0.006, 0.006) | (0.626, 0.638) | 0.9107 |
| S3 (p=2) | 19 | 19 | 19 | 19 | 0.04 | (-0.000, 0.000) | (0.632, 0.632) | 0.9107 |
| S3 (p=∞) | 19 | 19 | 19 | 19 | 0.04 | (-0.000, 0.000) | (0.632, 0.632) | 0.9107 |

Rozpiętość F1 punktów niezdominowanych (pokrycie frontu):

- S1: 17 różnych wartości F1, największa luka ΔF1 = 0.977
- S2: 21 różnych wartości F1, największa luka ΔF1 = 0.128
- S3 (p=2): 18 różnych wartości F1, największa luka ΔF1 = 0.482
- S3 (p=∞): 19 różnych wartości F1, największa luka ΔF1 = 0.083

### Wariant (a) – wielomiany 3. stopnia, U = koło [a_cubic_ball]

Tablica wypłat (wiersz i = F w minimum kryterium i):

| min | F1 | F2 |
|---|---|---|
| F1 | 0.0000 | 0.9990 |
| F2 | 0.6801 | -1.3365 |

ideal = (0.000, -1.337), nadir = (0.680, 0.999), r = (-0.068, -1.570)

| metoda | zadań | rozwiązanych | niezdominowanych | unikalnych niezd. | czas [s] | najbliższy idealnemu: x | F(x) | d₂ |
|---|---|---|---|---|---|---|---|---|
| S1 | 19 | 19 | 19 | 19 | 0.02 | (0.451, -0.535) | (0.248, -0.387) | 0.5465 |
| S2 | 21 | 21 | 21 | 21 | 0.04 | (0.454, -0.542) | (0.254, -0.402) | 0.5468 |
| S3 (p=2) | 19 | 19 | 19 | 19 | 0.05 | (0.470, -0.577) | (0.283, -0.489) | 0.5520 |
| S3 (p=∞) | 19 | 19 | 19 | 16 | 0.04 | (0.444, -0.520) | (0.237, -0.351) | 0.5468 |

### Wariant (b) – jedno kryterium liniowe, U = [-1, 1]³ [b_mixed]

Tablica wypłat (wiersz i = F w minimum kryterium i):

| min | F1 | F2 | F3 |
|---|---|---|---|
| F1 | 0.0000 | 1.9992 | -0.0003 |
| F2 | 1.9992 | 0.0000 | -0.0003 |
| F3 | 3.0000 | 5.0000 | -3.0000 |

ideal = (0.000, 0.000, -3.000), nadir = (3.000, 5.000, -0.000), r = (-0.300, -0.500, -3.300)

| metoda | zadań | rozwiązanych | niezdominowanych | unikalnych niezd. | czas [s] | najbliższy idealnemu: x | F(x) | d₂ |
|---|---|---|---|---|---|---|---|---|
| S1 | 36 | 36 | 36 | 36 | 0.02 | (-0.100, -0.500, -0.900) | (1.070, 1.470, -1.500) | 0.6810 |
| S2 | 121 | 80 | 80 | 41 | 0.32 | (-0.112, -0.500, -0.888) | (1.051, 1.500, -1.500) | 0.6802 |
| S3 (p=2) | 36 | 36 | 36 | 36 | 0.07 | (-0.193, -0.501, -0.809) | (0.942, 1.710, -1.502) | 0.6818 |
| S3 (p=∞) | 36 | 36 | 36 | 36 | 0.05 | (-0.050, -0.485, -0.920) | (1.084, 1.345, -1.455) | 0.6842 |

### Wariant (b) – jedno kryterium liniowe, U = kula [b_mixed_ball]

Tablica wypłat (wiersz i = F w minimum kryterium i):

| min | F1 | F2 | F3 |
|---|---|---|---|
| F1 | 0.0000 | 1.9992 | -0.0003 |
| F2 | 1.0000 | 0.1716 | -0.0005 |
| F3 | 1.0000 | 2.9984 | -1.7321 |

ideal = (0.000, 0.172, -1.732), nadir = (1.000, 2.998, -0.000), r = (-0.100, -0.111, -1.905)

| metoda | zadań | rozwiązanych | niezdominowanych | unikalnych niezd. | czas [s] | najbliższy idealnemu: x | F(x) | d₂ |
|---|---|---|---|---|---|---|---|---|
| S1 | 36 | 36 | 36 | 36 | 0.05 | (0.071, -0.214, -0.500) | (0.301, 1.158, -0.643) | 0.7797 |
| S2 | 121 | 88 | 88 | 68 | 0.38 | (-0.014, -0.289, -0.563) | (0.401, 1.302, -0.866) | 0.7555 |
| S3 (p=2) | 36 | 36 | 36 | 36 | 0.10 | (0.023, -0.274, -0.571) | (0.401, 1.213, -0.821) | 0.7572 |
| S3 (p=∞) | 36 | 36 | 36 | 36 | 0.07 | (-0.029, -0.318, -0.607) | (0.470, 1.315, -0.954) | 0.7659 |
