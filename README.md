# Análisis de sistemas electorales para la Cámara de Diputadas y Diputados de Chile

**Informe disponible en el *notebook* [`informe.ipynb`](informe.ipynb).**

A lo largo del presente repositorio se utilizan los datos de las elecciones
parlamentarias de Chile de 2021 para simular la distribución de escaños en la
Cámara de Diputadas y Diputados bajo distintos sistemas electorales. Para cada
sistema, se presentan ventajas y desventajas, justificadas mediante un análisis
numérico de los resultados (p.ej. cálculo del [índice de
Gallagher](https://en.wikipedia.org/wiki/Gallagher_index) y el [índice de
Loosemore-Hanby](https://en.wikipedia.org/wiki/Loosemore%E2%80%93Hanby_index)
para ver desproporcionalidad, y de la [cantidad efectiva de
partidos](https://en.wikipedia.org/wiki/Effective_number_of_parties) para ver
fragmentación).

Se analizan los siguientes sistemas:

- Proporcional por distrito
    - D'Hondt estándar (actualmente usado)
    - D'Hondt con umbral nacional (actualmente en discusión)
    - D'Hondt sin pactos
- Proporcional nacional
    - D'Hondt nacional con umbral bajo
    - D'Hondt nacional con umbral alto
- Bipartidistas
    - Binominal
    - Uninominal (first-past-the-post)
- Sistema biproporcional
- Mixtos
    - Mixed-Member Proportional / Additional-Member System
    - D'Hondt con escaños niveladores

Finalmente, en base a los resultados obtenidos en los análisis anteriores, se
propone un *framework* para un **nuevo sistema electoral**, que permita un
compromiso entre las características que se buscan en el contexto político
chileno. Este sistema es analizado un poco más en profundidad mediante el
ajuste de sus parámetros.