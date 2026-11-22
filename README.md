# optymalizacja-wielokryterialna-dyskretna

## Aplikacja do eksperymentów

Aplikacja webowa w Streamlit:

```bash
conda env create -f environment.yml   # albo: pip install -r requirements.txt
conda activate owd_env
streamlit run app.py
```

Struktura projektu:

```
app.py                      # punkt wejścia Streamlit (menu stron)
run_experiments.py          # seria eksperymentów z pliku TOML (ćw. 2)
owd/                        # logika – bez Streamlita
├── algorithms/             # algorytmy wyznaczania punktów niezdominowanych (ćw. 1)
│   ├── without_filtration.py
│   ├── with_filtration.py
│   └── ideal_point.py
├── data_generator.py       # generowanie danych
├── experiment.py           # uruchamianie algorytmów, kontrola poprawności, CSV
└── scalarization.py        # metody skalaryzacji S1, S2, S3 (ćw. 3)
app_pages/                  # strony aplikacji
├── single_experiment.py    # „Pojedynczy eksperyment”
├── experiment_results.py   # „Wyniki eksperymentów”
└── scalarization_page.py   # „Skalaryzacja”
experiments/                # konfiguracje serii eksperymentów (TOML)
results/                    # wyniki run_experiments.py (poza repozytorium)
reports/                    # sprawozdania – katalog na ćwiczenie
└── cw3_skalaryzacja/       # sprawozdanie.md, generate_figures.py, fig/
```

- wybór algorytmu (`ALGORITHMS` w `owd/experiment.py` – nowy algorytm wystarczy tam dopisać jako `(moduł, nazwa_funkcji)`; funkcja przyjmuje listę punktów, minimalizuje wszystkie kryteria, zwraca `(P, liczba_porownan)` i porównuje punkty przez funkcję `dominates(a, b)` z poziomu modułu – aplikacja ją opakowuje, żeby zaznaczyć punkty porównywane),
- generowanie danych: rozkład, N, k, seed (🎲 losuje), zaokrąglanie,
- tabela „Kryteria”: dla każdego kryterium osobno nazwa (domyślnie `f1`, `f2`, ...), kierunek `min`/`max`, średnia i odchylenie; przyciski „Wszystkie min/max”; kryteria maksymalizowane są negowane przed wywołaniem algorytmu,
- „Uruchom wszystkie algorytmy” liczy wszystkie algorytmy na tych samych danych – wybór algorytmu w panelu przełącza wtedy wyświetlany wynik,
- automatyczna kontrola poprawności: każdy wynik jest sprawdzany z definicją dominacji Pareto,
- interaktywny wykres (Plotly) 2D lub 3D z warstwami: zbiór punktów, zdominowane, porównywane, niezdominowane, front Pareto, punkt idealny; po najechaniu widać id i wartości punktu,
- tabela danych ze statusem i liczbą porównań, filtrem i sortowaniem; zaznaczone wiersze są wyróżnione na wykresie,
- pobieranie CSV: wynik (każdy punkt ze statusem i liczbą porównań; metadane w liniach `#` – `pd.read_csv(path, comment="#")`) oraz historia wszystkich uruchomień,
- czas w wynikach mierzony jest w osobnym przebiegu bez zliczania porównań poszczególnych punktów, żeby to śledzenie nie zawyżało czasu.

## Strona „Skalaryzacja” (ćw. 3)

Wyznacza dyskretną aproksymację zbiorów P(U, F) i F_P(U) dla problemów ciągłych trzema metodami skalaryzacji. Obliczenia są w `owd/scalarization.py`, strona w `app_pages/scalarization_page.py`.

- warianty: (a) dwa kryteria F = (F1, F2), U ⊂ R² – **obie funkcje muszą być nieliniowe**; (b) trzy kryteria G = (G1, G2, G3), U ⊂ R³ – **co najmniej jedno nieliniowe**; inaczej przycisk „Oblicz” jest zablokowany. Funkcja jest liniowa (afiniczna) wtedy i tylko wtedy, gdy jej hesjan jest macierzą zerową (`sympy.hessian`),
- U: kwadrat / sześcian [-1, 1]ⁿ albo koło / kula ‖x‖ ≤ 1,
- funkcje wpisuje się jako formuły w `x1, x2, x3` (np. `x1**2 + (x2 - 1)**2`, `exp(-x1)`); domyślne funkcje i gotowe zestawy (także z niewypukłym frontem) wybiera się z listy, „↺” przywraca domyślne; każda funkcja jest wyświetlana jako wzór (LaTeX) z oznaczeniem liniowa / nieliniowa,
- metody: S1 – funkcja liniowa Σλᵢ Fᵢ, S2 – ε-ograniczenia (min Fⱼ przy Fᵢ ≤ aᵢ, aᵢ między punktem idealnym a nadirem z tablicy wypłat), S3 – odległość od punktu r dominującego F_P(U) w metryce ‖·‖_{λ,p} (p = 1, 2, ∞ lub dowolne p ≥ 1; dla p = ∞ zadanie jest przekształcone do postaci gładkiej ze zmienną pomocniczą),
- parametry λ (S1, S3) i z (S2) zmieniają się ze stałym krokiem ρ; każde zadanie rozwiązuje `scipy.optimize.minimize` metodą SLSQP (odpowiednik `fmincon`; ograniczenia jako `Bounds` / `NonlinearConstraint`, gradienty z `sympy` – jakobian `Matrix.jacobian`, w S3 przybliżany przez scipy) z kilku punktów startowych (hipersześcian łaciński `scipy.stats.qmc` w kwadracie, `scipy.stats.uniform_direction` w kole) oraz z rozwiązania dla poprzedniego parametru; metryki i odległości liczy `numpy.linalg.norm`,
- otrzymane punkty przechodzą przez nasz algorytm filtracji (`owd/algorithms/with_filtration.py`) – punkty zdominowane (np. słabo niezdominowane z S2 / S3) są oznaczone,
- rozwiązanie najbliższe punktowi idealnemu: wybierane spośród niezdominowanych, odległość liczona po normalizacji każdego kryterium zakresem (nadir − punkt idealny), metryka p = 1, 2 lub ∞; zaznaczone na obu wykresach (z linią od punktu idealnego) i w tabeli, a w „Porównaniu metod” zestawione dla wszystkich metod,
- tabela wyników: numer zadania, parametry, x, F(x), wartość funkcji skalaryzującej, odległość od punktu idealnego, „jak nr” (ten sam punkt co we wcześniejszym zadaniu) i status; filtr (np. tylko unikalne niezdominowane), zaznaczone wiersze są wyróżnione na wykresach,
- wykresy w przestrzeni ocen i decyzji: w 2D spline (`scipy.interpolate.PchipInterpolator` w przestrzeni ocen, `make_interp_spline` w przestrzeni decyzji), w 3D powierzchnia z triangulacji siatki parametrów (`scipy.spatial.Delaunay`); zakładka „Porównanie metod” nakłada wyniki wszystkich metod; tabela wyników i CSV.

## Seria eksperymentów (skrypt, ćw. 2)

Eksperymenty na wielu losowych zbiorach uruchamia `run_experiments.py`. Parametry są w pliku TOML – nowy eksperyment to kopia `experiments/eksperyment_B.toml` ze zmienionymi wartościami (rozmiary N i liczba powtórzeń K, listy liczby kryteriów, algorytmy, rozkład, kierunki, seed).

```bash
python run_experiments.py experiments/eksperyment_B.toml --quick   # szybki test: 2 powtórzenia na rozmiar
python run_experiments.py experiments/eksperyment_B.toml           # pełny eksperyment
python run_experiments.py --resume results/<katalog>               # dokończenie przerwanego (Ctrl+C)
python run_experiments.py --summarize results/<katalog>            # ponowne przeliczenie podsumowania
```

- każdy zbiór danych ma seed wyliczony z `base_seed`, N, k i numeru powtórzenia – ten sam plik konfiguracji zawsze daje te same dane, a wszystkie algorytmy liczą na identycznych zbiorach,
- wyniki trafiają do `results/<nazwa>_<data>/`: `raw.csv` (każdy zbiór × algorytm, zapisywany na bieżąco), `summary.csv` (średnia, odchylenie, min, max dla N, k, algorytmu), `comparisons.csv` i `time_ms.csv` (tabele średnich do sprawozdania), `report.html` (wykresy),
- `--workers N` liczy kilka zbiorów równolegle – liczba porównań się nie zmienia, ale czasy są wtedy mniej wiarygodne.

### Strona „Wyniki eksperymentów”

Przegląda wyniki `run_experiments.py`: wybór katalogu z `results/` (albo wczytanie pliku `raw.csv`), filtry algorytmów, N i k oraz opcja pominięcia niepoprawnych wyników. Zakładki:

- **Wykresy** – liczba porównań, czas lub |P| w zależności od k albo N (średnia ± odchylenie albo wykres pudełkowy, skala logarytmiczna),
- **Tabele** – wybrana statystyka (średnia, odchylenie, mediana, min, max) dla N, k i algorytmu, także jako iloraz względem wybranego algorytmu; pobieranie CSV,
- **Liczba niezdominowanych** – |P| (albo |P|/N) na zbiór danych, liczone z poprawnych wyników,
- **Poprawność** – odsetek wyników zgodnych z definicją dominacji Pareto,
- **Dane surowe** – przefiltrowany `raw.csv`.

Strona czyta `raw.csv` na bieżąco – przycisk „Odśwież” pokazuje postęp trwającego eksperymentu.

## Sprawozdania

Każde ćwiczenie ma katalog w `reports/` ze sprawozdaniem w Markdownie, rysunkami i skryptem, który je generuje.

```bash
python reports/cw3_skalaryzacja/generate_figures.py          # rysunki fig/ i tabele results.md (wymaga matplotlib)
cd reports/cw3_skalaryzacja && pandoc sprawozdanie.md -o sprawozdanie.pdf --pdf-engine=xelatex -V geometry:margin=2cm
```
