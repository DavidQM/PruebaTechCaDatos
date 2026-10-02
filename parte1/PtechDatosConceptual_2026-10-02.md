# Prueba técnica · Parte 1

## Conceptual, estándares y arquitectura

| | |
|---|---|
| **Autor** | Alexis David Quintero Montoya |
| **Fecha** | 2 de octubre de 2026 |
| **Documento** | Prueba técnica · Profesional en Ciencia de Datos SIATA · Calidad de datos |

---

### Contenido

- [Pregunta 1.1](#pregunta-11): Calidad del dato hidrometeorológico
- [Pregunta 1.2](#pregunta-12): Niveles del dato y diccionario de banderas
- [Pregunta 1.3](#pregunta-13): Indicadores de calidad
- [Pregunta 1.4](#pregunta-14): Pseudocódigo para completitud

---

## Pregunta 1.1

> Explique con sus propias palabras qué significa que un dato hidrometeorológico *tenga calidad*. Su respuesta debe:

### (a)

Para definir ello hay que partir de la definición de **Calidad de datos**. Entiéndase por Calidad de Datos como una propiedad del registro y de la serie: si está completo, si es válido, si es coherente con otras variables y si llegó a tiempo. Se evalúa mirando los datos, es decir cumple con unas dimensiones mínimas.

Un dato hidrometeorológico tiene calidad cuando ha pasado por un proceso de análisis basado en dimensiones y que cumple con unas políticas de uso dispuestas por personal experto. Con ello, se listarán unas reglas para etiquetar los registros y/o las series del registro hidrometeorológico. En el contexto SIATA, mínimamente deberá cumplir con las dimensiones completitud, consistencia, exactitud y unicidad (ISO/IEC 25012).

### (b)

Anteriormente se expuso una definición corta de calidad de dato, por otra parte, **la calidad de un instrumento** está relacionada a la capacidad del sensor de medir bien: exactitud, calibración, trazabilidad metrológica y emplazamiento. Lo cual depende no solo de los materiales del sensor, sino de las condiciones externas que lo rodean, el lugar donde fue aprovisionada y la periodicidad con la que se le realizan los mantenimientos.

Por otra parte, **la calidad del proceso** es qué tan confiable es la forma de producir el dato. En el contexto SIATA está relacionada con el grupo como Mantenimiento, Telemetría, Calidad del Dato y Sistemas. Todos ellos interactúan en la producción y disposición de los datos de las estaciones hidrometeorológicas.

### (c)

Un buen ejemplo de ello son las pruebas de rango (Deltas) solo descartan lo imposible, no lo equivocado. Un valor puede estar dentro de todos los límites físicos y estar mal.

**Ejemplo:**

- Un pluviómetro con un obstáculo cerca, o desnivelado, mide menos de lo que llovió sin dar nunca un valor absurdo.

---

## Pregunta 1.2

> Complete la siguiente tabla y justifique en un párrafo la diferencia entre cada nivel.

| Aspecto | Dato crudo | Dato validado con metadatos | Dato para consumo del usuario |
|---|---|---|---|
| **Definición** | Lectura tal como llega del sensor, sin modificar. | Registro que pasó la cadena de pruebas, con banderas y metadatos. Equivale a los estados en revisión o definitivo, según la revisión humana. | Producto derivado de datos validados, agregado (10 minutos, horario, diario) y documentado para un uso concreto. |
| **Unidades y resolución** | Las del sensor y su resolución nativa, sin conversiones (por ejemplo 0,254 mm por pulso en una estación pluvio). | Unidades normalizadas y declaradas en los metadatos, con resolución e incertidumbre documentadas. | Unidades explícitas en cada campo y redondeo acorde con la resolución del instrumento, sin decimales de más. |
| **Banderas de calidad** | Dependerá de las dimensiones de calidad que se tengan en las políticas. | Banderas propias con diccionario y versión del algoritmo (ok, faltante, fuera de rango, sospechoso, incoherente), no excluyentes entre sí. | Resumen por valor agregado: porcentaje de datos válidos que lo respaldan y estado del dato. |
| **Metadatos mínimos** | Identificador de estación, fecha y hora, valor y fecha y hora de arribo. | Los anteriores más coordenadas, altura del sensor, unidades, resolución, intervalo de muestreo, zona horaria, versión del algoritmo de validación, diccionario de banderas e historial de mantenimiento y calibración. | Los de la columna anterior que el usuario necesita para interpretar el valor, más licencia de uso, fecha de generación, versión del producto, contacto y forma de citar (ISO 19115). |
| **Tratamiento de faltantes** | No se rellenan. Las mediciones (-999) se conservan como llegaron. | El analista pasa a nulo con bandera de faltante. Se imputa solo en huecos cortos y siempre marcado como imputado. | Si no se cumple el mínimo de datos definido en la política de uso para esa serie, el agregado no se calcula. |
| **Usuario típico** | Técnicos de operación y mantenimiento, verificadores. | Grupo de geociencias y sus analistas, científicos de datos, modeladores. | Ciudadanía, entidades que toman decisiones, investigadores y desarrolladores que consumen por API. |
| **Riesgo de mal uso** | Decidir con datos que nadie ha verificado. | Usar la serie sin leer las banderas o sin saber que cambió la versión de validación. | Leer un agregado como si fuera una medición directa, ignorar la cobertura o comparar estaciones con metodologías distintas. |
| **Qué se busca** | Conservar la evidencia para poder volver al original y reprocesar. | Confiabilidad y trazabilidad. | Claridad y reutilización. |

### Diccionario de banderas

A continuación se propone un diccionario de banderas para este ejercicio:

| Código | Bandera | Cuándo se asigna | Qué pasa con el dato |
|:---:|---|---|---|
| `OK` | Aprobado | Pasó todas las pruebas aplicables | Entra a los productos |
| `F` | Faltante | No hay lectura en la ventana esperada, o llegó el -999 | No entra a agregados; cuenta contra la cobertura |
| `X` | Formato inválido | Fecha ilegible o valor no numérico | Se excluye; se conserva el crudo |
| `D` | Duplicado | Más de un registro para la misma estación, variable y ventana | Se conserva uno; el resto va a cuarentena |
| `R` | Fuera de rango | Excede los límites físicos de la variable | Se excluye, con causa registrada |
| `S` | Sospechoso | Pasa el rango, pero falla salto o persistencia | Se conserva marcado, para revisión humana |
| `I` | Incoherente | Contradice otra variable o un sensor redundante | Se conserva marcado, para revisión humana |
| `E` | Inconsistencia espacial | No concuerda con estaciones vecinas | Revisión humana, sin rechazo automático |
| `M` | Imputado | Valor estimado en un hueco corto | Se publica solo si el producto lo permite, marcado |
| `C` | Corregido | Modificado tras la verificación | Se conserva el original junto al corregido |
| `K` | Mantenimiento | Lectura durante una intervención declarada | Se excluye de productos y de la cobertura |

---

## Pregunta 1.3

Propongo tres indicadores de dimensiones distintas: completitud, validez y consistencia. Los tres se expresan como porcentaje o coeficiente, se calculan por estación y variable, y se pueden comparar entre estaciones y entre variables.

### Indicador A · Cobertura de datos válidos

| Campo | Descripción |
|---|---|
| **Nombre y dimensión** | Cobertura de datos válidos. Dimensión: completitud. |
| **Objetivo y pregunta** | Responde a la pregunta de qué proporción de las lecturas que la estación debía entregar en el periodo llegó y se puede usar. |
| **Fórmula** | `C = (N_válidos / N_esperados) × 100`<br><br>**N_esperados** es la duración efectiva del periodo en segundos dividida por el intervalo de muestreo de la estación y variable, descontando mantenimientos declarados.<br><br>**N_válidos** es el número de ventanas de muestreo distintas con valor no nulo, distinto del centinela (-999) y con una bandera aceptada. |
| **Unidad y periodicidad** | Porcentaje. Se calcula por hora (operación), por día y por mes. |
| **Fuente y agregación** | Tabla de registros (`codigo`, `fecha_hora`, `valor`, `calidad`) y tabla de metadatos con el intervalo de cada serie. Se calcula por estación y variable, y luego se resume por red. |
| **Umbrales** | 🟢 **Aceptable:** 95 % o más.<br>🟡 **Alerta:** de 67 % a menos de 95 %.<br>🔴 **Crítico:** menos de 67 %.<br><br>El 67 % equivale a 20 de 30 días, el mínimo que exige el IDEAM para calcular un total mensual de precipitación. Por debajo de eso el agregado oficial no existe. El 95 % es una propuesta mía para una red de alerta y habría que calibrarlo con la historia de la red. En los archivos de prueba la cobertura va de 98,8 % a 100 %. |
| **Responsable y acción** | Responsable de operación de la red, con el custodio de la variable. Ante alerta, revisar comunicación y alimentación del equipo. Ante crítico durante dos días seguidos, programar visita de mantenimiento y recuperar los datos de la memoria local del datalogger. |
| **Limitaciones y sesgos** | Depende de que el intervalo de muestreo esté bien registrado. No distingue un hueco de un día de mil minutos sueltos, así que conviene reportar también la duración del hueco más largo. Qué cuenta como válido depende del diccionario de calidad: en el nivel 803, por ejemplo, la cobertura es de 99,9 % contando solo presencia y de 87,9 % si se exige calidad igual a 1. |

### Indicador B · Validez por rango y centinelas

| Campo | Descripción |
|---|---|
| **Nombre y dimensión** | Tasa de registros válidos por rango. Dimensión: validez. |
| **Objetivo y pregunta** | Responde a la pregunta de qué proporción de los registros recibidos cumple formato, tipo y rango físico permitido para la variable. |
| **Fórmula** | `V = (N_en_rango / N_recibidos) × 100`<br><br>`N_recibidos` son los registros con marca de tiempo válida. `N_en_rango` son los que además tienen valor numérico, no son centinela y quedan dentro de `[mínimo, máximo]` de la variable.<br><br>Los límites se tabulan por variable en el catálogo de reglas y se revisan con percentiles históricos. |
| **Unidad y periodicidad** | Porcentaje. Cálculo diario y seguimiento semanal de la tendencia. |
| **Fuente y agregación** | Tabla de crudos o validados. Se calcula por estación y variable, y se resume por red. |
| **Umbrales** | 🟢 **Aceptable:** 99 % o más.<br>🟡 **Alerta:** de 95 % a menos de 99 %.<br>🔴 **Crítico:** menos de 95 %.<br><br>Son propuestas. Como los rangos solo descartan lo imposible, el rechazo esperado es muy bajo y un 1 % ya merece mirar el sensor. En la estación 35 hay 7 registros con -999 en cada canal. En el piranómetro 6004 hay 98 minutos por encima de 1 361 W/m², un valor cercano a la irradiancia solar fuera de la atmósfera (\*), que serían candidatos a revisión. |
| **Responsable y acción** | Custodio de la variable y técnico de mantenimiento. Revisar sensor y cableado, corregir el rango si estaba mal definido y marcar los registros como rechazados con su causa. |
| **Limitaciones y sesgos** | Un rango amplio deja pasar errores plausibles (ver 1.1 c). Un rango demasiado estrecho rechaza extremos reales, justo en los eventos que más interesan, y sesga el indicador a favor de condiciones normales. Los límites dependen de la estación, por ejemplo de su altitud. |

### Indicador C · Concordancia con la estación de referencia

| Campo | Descripción |
|---|---|
| **Nombre y dimensión** | Concordancia con referencia. Dimensión: consistencia (con elementos de exactitud relativa). |
| **Objetivo y pregunta** | Responde a la pregunta de si la serie se comporta como lo hacen sus referencias: estaciones vecinas de entorno similar, un sensor redundante, o una estación aguas arriba o abajo en la misma corriente. |
| **Fórmula** | `r = corr(serie_estación, serie_referencia)`<br><br>Correlación de Pearson sobre series diarias en ventana móvil de 30 días. Para nivel se aplica el rezago que maximiza la correlación, que representa el desfase de la escorrentía. Con sensores redundantes se reporta además el cociente de acumulados diarios. |
| **Unidad y periodicidad** | Coeficiente adimensional entre -1 y 1. Se recalcula cada día con la ventana de 30 días. |
| **Fuente y agregación** | Series ya validadas por las pruebas de formato y rango, agregadas a diario, y una tabla de pares estación-referencia mantenida en los metadatos. Se calcula por par y se resume por variable. |
| **Umbrales** | 🟢 **Aceptable:** r de 0,7 o más.<br>🟡 **Alerta:** de 0,5 a menos de 0,7.<br>🔴 **Crítico:** menos de 0,5.<br><br>El manual de validación hidrológica del IDEAM usa 0,7 como criterio para comparar estaciones cuando no hay balance en la misma corriente. El de meteorología usa r mayor que 0,5 como condición mínima para corregir o generar datos con una ecuación de regresión. Tomé esos dos valores como bordes. |
| **Responsable y acción** | Analista de datos. Contrastar con la hoja de inspección y el historial de mantenimiento, revisar la cota cero o el emplazamiento y dejar el dato en estado en revisión hasta aclarar. |
| **Limitaciones y sesgos** | La lluvia convectiva es muy local, así que dos pluviómetros cercanos pueden diferir sin que haya error. La correlación no detecta un sesgo constante: un sensor que mide 20 % de menos puede correlacionar perfectamente, por eso se acompaña del cociente o del sesgo medio. Depende de elegir bien las referencias. |

---

## Pregunta 1.4

### Pseudocódigo para completitud

```text
ENTRADAS
  registros : tabla (codigo, fecha_hora, valor, calidad)
  tramos    : tabla (codigo, variable, intervalo_s, fecha_inicio, fecha_fin)
              una fila por tramo de operación con el mismo intervalo de muestreo
  periodo   : (inicio, fin)          # fin es exclusivo
  calidad_ok: lista de banderas aceptadas
  centinela : -999

FUNCION cobertura(registros, tramos, periodo, calidad_ok, centinela)

  PARA CADA mes M DENTRO DE periodo:
    PARA CADA (codigo, variable) EN tramos:

      esperados ← 0
      validos   ← 0

      PARA CADA tramo T DE ESA SERIE:
        desde ← MAX(inicio_de(M), T.fecha_inicio)
        hasta ← MIN(fin_de(M), T.fecha_fin)
        SI hasta <= desde: CONTINUAR

        # ventanas que la estación debía entregar en el tramo
        esperados ← esperados + ENTERO( segundos(hasta - desde) / T.intervalo_s )
        esperados ← esperados - ventanas_en_mantenimiento(T, desde, hasta)

        # registros del tramo
        d ← registros DONDE codigo = T.codigo
                       Y fecha_hora >= desde Y fecha_hora < hasta

        # un registro es válido si tiene valor útil y bandera aceptada
        d ← d DONDE valor NO ES NULO
                Y valor <> centinela
                Y calidad EN calidad_ok

        # cada registro cae en su ventana; los duplicados cuentan una vez
        ventanas ← CONJUNTO( PISO(fecha_hora / T.intervalo_s) PARA CADA registro EN d )
        validos  ← validos + TAMAÑO(ventanas)
      FIN PARA

      SI esperados > 0:
        cobertura ← 100 * validos / esperados
      SINO:
        cobertura ← NULO

      GUARDAR (codigo, variable, M, esperados, validos, cobertura)
    FIN PARA
  FIN PARA

  # resumen por red: se suman esperados y validos, NUNCA se promedian porcentajes
  RETORNAR tabla de resultados
FIN FUNCION
```

---
