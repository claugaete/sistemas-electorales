# Análisis de sistemas electorales para la Cámara de Diputadas y Diputados de Chile

**Informe disponible en el *notebook* [`informe.ipynb`](informe.ipynb).**

A lo largo del presente repositorio se utilizan los datos de las elecciones
parlamentarias de Chile de 2021 para simular la distribución de escaños en la
Cámara de Diputadas y Diputados bajo distintos sistemas electorales. Para cada
sistema, se presentan ventajas y desventajas, justificadas mediante un análisis
numérico de los resultados (e.g. cálculo del [índice de
Gallagher](https://en.wikipedia.org/wiki/Gallagher_index) y el [índice de
Loosemore-Hanby](https://en.wikipedia.org/wiki/Loosemore%E2%80%93Hanby_index)
para ver disproporcionalidad, y de la [cantidad efectiva de
partidos](https://en.wikipedia.org/wiki/Effective_number_of_parties) para ver
fragmentación).

Se proponen los siguientes sistemas (en cursiva los que aún no se han
analizado):

- Proporcional por distrito
    - D'Hondt estándar (actualmente usado)
    - *D'Hondt con umbral nacional (actualmente en discusión)*
    - *D'Hondt sin pactos*
- *Proporcional nacional*
    - *D'Hondt nacional*
    - *D'Hondt nacional con umbral electoral*
    - *D'Hondt nacional sin pactos*
- *Bipartidistas*
    - *Binomial*
    - *Uninominal (first-past-the-post)*
- *Mixtos*
    - *Mixed-Member Proportional / Additional Member System*
    - *D'Hondt con escaños niveladores*
    - *Binomial con escaños niveladores*
- *Sistema biproporcional*