#  Simulador de Pandemia COVID-19 - Guía de Uso

##  Descripción

Simulador epidemiológico interactivo que modela la propagación de COVID-19 Evolucionado en Europa usando ecuaciones diferenciales (Modelo SIRD) y algoritmos de Divide y Vencerás para análisis de datos.

---

##  Instalación y Ejecución

### Requisitos Previos

```bash
pip install geopandas matplotlib numpy tkinter
```

**Nota:** `tkinter` viene preinstalado con Python en la mayoría de sistemas.

### Ejecutar el Simulador

```bash
python simuladorPandemia.py
```

---

##  Cómo Usar el Simulador

### Paso 1: Seleccionar País Inicial 🌍

1. En el panel izquierdo, usa el **dropdown** para elegir el país donde iniciará el brote
2. Puedes elegir cualquiera de los 35 países europeos disponibles
3. Ejemplo: Selecciona "Italy", "Germany" o "Poland"

### Paso 2: Iniciar Simulación 

1. Presiona el botón **"▶️ Iniciar Simulación"**
2. Espera unos segundos mientras se calculan los 300 días de propagación
3. La simulación comenzará automáticamente en modo Play

**Mensaje inicial:** Verás una ventana emergente con información del país seleccionado y una nota importante: *"La recuperación inicia alrededor del día 100"*

### Paso 3: Visualizar la Propagación 

#### Controles de Visualización (Panel Derecho - Superior)

Haz clic en los botones para cambiar el tipo de mapa:

- **🔴 Infectados**: Muestra casos activos (mapa rojo)
- **⚫ Muertes**: Muestra fallecidos acumulados (mapa gris)
- **🟢 Recuperados**: Muestra recuperaciones (mapa verde)
- **🔵 Susceptibles**: Muestra población vulnerable (mapa azul)
- **⏸️ Pausar / ▶️ Play**: Controla la animación

#### Navegación Temporal

- **Slider horizontal**: Arrastra para moverte entre los días 0-300
- **Contador**: Muestra "Día: X / 300"
- Puedes usar el slider mientras la simulación está pausada o en reproducción

### Paso 4: Analizar Resultados 

#### Botones de Análisis (Panel Izquierdo - Centro)

Usa estos botones para obtener estadísticas usando algoritmos de Divide y Vencerás:

1. **🔴 Mayor Infectados**: Encuentra el país con el pico más alto de infectados
2. **⚫ Mayor Muertes**: Encuentra el país con más fallecidos
3. **🟢 Mayor Recuperados**: Encuentra el país con más recuperaciones
4. **📅 Día Pico**: Identifica el día con más infectados en toda Europa
5. **🔽 Ordenar por Infectados**: Ordena la tabla de mayor a menor infectados
6. **🔽 Ordenar por Muertes**: Ordena la tabla de mayor a menor muertes
7. **🔽 Ordenar por Recuperados**: Ordena la tabla de mayor a menor recuperados

Cada botón mostrará una ventana emergente con los resultados del análisis.

### Paso 5: Revisar la Tabla de Datos �

La tabla (panel izquierdo - inferior) muestra en tiempo real:

- **País**: Nombre del país
- **Sistema**: 🟢 Desarrollado | 🟡 Normal | 🔴 Precario
- **Infectados**: Casos activos en el día actual
- **Muertes**: Fallecidos acumulados
- **Recuperados**: Personas recuperadas

**Fila TOTAL**: Al final de la tabla, suma de todos los países

---

##  Eventos Importantes Durante la Simulación

###  Día ~100: Fase de Recuperación

Cuando llegues al día 95-110, verás:

- **Título del mapa cambia a verde**: "🏥 ¡FASE DE RECUPERACIÓN INICIADA! 🏥"
- **Mensaje en consola**:
  ```
  ********************************************************
  🏥 DÍA 100: ¡FASE DE RECUPERACIÓN INICIADA!
  Los sistemas de salud comienzan a recuperar pacientes
  ********************************************************
  ```
- **Efecto visible**: Las curvas de infectados comienzan a descender

###  Día 110-140: Inicio de Vacunación

Los países comienzan a vacunar según su sistema:
- 🟢 Desarrollados: Día 110
- 🟡 Normales: Día 120
- 🔴 Precarios: Día 140

---

¡Disfruta explorando la dinámica de propagación de pandemias! 🦠📊
