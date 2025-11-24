# Modelamiento Estadístico

Este repositorio contiene un proyecto de modelamiento estadístico aplicado sobre datos reales / simulados con el objetivo de analizar fenómenos, estimar relaciones e interpretar resultados usando técnicas estadísticas avanzadas.

## Objetivos

- Construir modelos estadísticos para explicar / predecir variables de interés.  
- Realizar inferencia sobre los parámetros de los modelos y validar supuestos estadísticos (normalidad, heterocedasticidad, multicolinealidad, etc.).  
- Evaluar la capacidad predictiva de los modelos y hacer comparaciones entre diferentes aproximaciones.

## Modelos y Técnicas Utilizadas

En este proyecto se han empleado varios modelos estadísticos, entre ellos:

- **Regresión lineal simple y múltiple**: para capturar relaciones lineales entre variables dependientes e independientes.  
- **Modelos de regresión con variables categóricas o transformadas**: para acomodar datos no lineales o variables dummy.  
- **Análisis de series de tiempo** (si aplica): para modelar dinámicas temporales entre variables.  
- **Modelos de datos de panel** (si se usan paneles): para controlar la variabilidad no observada entre unidades.  
- **Modelos de regresión con errores robustos**: para corregir heterocedasticidad o dependencia serial cuando sea necesario.  
- **Validación y diagnóstico**: pruebas de normalidad, heterocedasticidad, multicolinealidad, autocorrelación; así como métricas de ajuste.

## Flujo del Proyecto

1. **Carga y preprocesamiento de datos**: se limpian los datos, se crean variables derivadas y se preparan conjuntos de entrenamiento / prueba (si aplica).  
2. **Análisis exploratorio**: estadísticas descriptivas, gráficos, correlaciones, matrices de dispersión.  
3. **Estimación de modelos**: se ajustan los modelos descritos arriba según el problema y los datos.  
4. **Diagnóstico de modelos**: evaluación de supuestos, tests estadísticos y análisis de residuos.  
5. **Validación**: estimación de errores de predicción, validación cruzada o partición de datos, comparación entre modelos.  
6. **Interpretación y conclusiones**: análisis de los parámetros estimados, su significancia, implicaciones prácticas y teóricas.  

## Uso

1. Clona este repositorio:
   ```bash
   git clone https://github.com/YokoMolina/Modelamiento-Estadistico.git
