#Graficos ggplot2

library(ggplot2)
library(pheatmap)
library(ggrepel)
data(iris)  # Cargamos la base de datos iris
head(iris)  # Vista previa de los datos

#1. Gráfico de Dispersión (Scatter Plot)

ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Gráfico de Dispersión de Sepal.Length vs Sepal.Width",
       x = "Longitud del Sépalo",
       y = "Ancho del Sépalo")

#2. Histograma

ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
  geom_histogram(binwidth = 0.3, alpha = 0.7, position = "identity") +
  theme_minimal() +
  labs(title = "Distribución de la Longitud del Sépalo",
       x = "Longitud del Sépalo",
       y = "Frecuencia")
#3. Boxplot (Diagrama de Cajas)

ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribución de Sepal.Length por Especie",
       x = "Especie",
       y = "Longitud del Sépalo")

#4. Gráfico de Barras

ggplot(iris, aes(x = Species, fill = Species)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Cantidad de Observaciones por Especie",
       x = "Especie",
       y = "Conteo")

#5. Gráfico de líneas

ggplot(iris, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Relación entre Sepal.Length y Petal.Length",
       x = "Longitud del Sépalo",
       y = "Longitud del Pétalo")

#6 Gráfico de Densidad

ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Densidad de la Longitud del Sépalo por Especie",
       x = "Longitud del Sépalo",
       y = "Densidad")

#7 Gráfico de Facetas

ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
geom_point() +
  facet_wrap(~ Species) +
  theme_minimal() +
  labs(title = "Gráfico de Dispersión por Especie",
       x = "Longitud del Sépalo",
       y = "Ancho del Sépalo")

#8 Gráfico de Violín

ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species)) +
  geom_violin(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribución de Sepal.Length por Especie",
       x = "Especie",
       y = "Longitud del Sépalo")

