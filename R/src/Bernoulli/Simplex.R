library(archetypes)
library(dplyr)

source("create_hex_data.R")

radius=10
z<-  Z # toy[c(7,23,247,298,359),]
k=4

#women
#z=toy[c(60,64,85,152),]
#k=4


params <- t(z)  #z
coefs <- t(S) # t(ae) #alphas

order = NULL
labels_cex = 1
labels = NULL
show_labels = TRUE
points_col = "#00000044"
points_pch = 19
points_cex = 1
projection = simplex_projection
show_points = TRUE
show_circle = TRUE
circle_col = "lightgray"
show_edges = TRUE
edges_col = "lightgray"
show_direction = TRUE
direction_length = 1
directions_col = points_col

object=c()
object$alphas = coefs

proj_z <- projection(params, r = radius - 1)

print(dim(proj_z))
print(dim(coefs))
proj_h <- coefs %*% proj_z

# Extract the relevant feature names for each sample
col_names <- col_names[col_names != "drug_name"]

hex_data <- create_hex_data(X, proj_h, col_names, num_bins = 30, selected_feature = NULL)




proj_labels <- proj_z
t <- cbind(x = acos(proj_z[, "x"] / (radius-1)), y = asin(proj_z[, "y"] / (radius-1)))
proj_labels <- cbind(x = radius * cos(t[, "x"]), y = radius * sin(t[, "y"]))

proj_circle <- list(center = cbind(x = 0, y = 0), radius = radius - 1)
proj_edges <- proj_z[as.integer(combn(1:k, 2)), ]

proj_directions <- vector("list", length = nrow(object$alphas))

for ( j in 1:nrow(object$alphas)) {
  s <- proj_h[j, , drop = FALSE]
  d <- matrix(NA_real_, ncol = 2, nrow = ncol(object$alphas))
  for ( i in 1:ncol(object$alphas) ) {
    e <- proj_z[i, , drop = FALSE]

    v <- e - s
    m <- sqrt(sum(v^2))
    v <- v / m

    px <- s[1] + v[1] * direction_length * object$alphas[j, i]
    py <- s[2] + v[2] * direction_length * object$alphas[j, i]

    d[i, ] <- c(px, py)
  }
  proj_directions[[j]] <- list(s = s, e = d)
}

if ( is.null(order) )
  order <- 1:k

if ( is.null(labels) )
  labels <- sprintf("A%s", order)

if ( length(points_col) == 1 )
  points_col <- rep(points_col, nrow(coefs))

if ( length(points_cex) == 1 )
  points_cex <- rep(points_cex, nrow(coefs))

if ( length(directions_col) == 1)
  directions_col <- rep(directions_col, nrow(coefs))


library(ggplot2)
library(grid)


############################# TRYING SOMETHING ELSE



# Load the necessary package
library(ggplot2)
library(gridExtra)

# Make sure that proj_z is a data frame with a column named x
proj_z <- as.data.frame(proj_z)

# Create the base plot
p <- ggplot(proj_z, aes(x = x, y = y)) +
  geom_point(color = "transparent") +
  coord_fixed(xlim = c(-radius, radius), ylim = c(-radius, radius)) +
  theme_void()

print(p)

if (show_circle) {
  # Extract the x and y coordinates of the center of the circle
  x0 <- proj_circle$center[1]
  y0 <- proj_circle$center[2]

  # Create a sequence of points to draw a circle
  theta <- seq(0, 2 * pi, length.out = 100)
  x <- x0 + (radius - 1) * cos(theta)
  y <- y0 + (radius - 1) * sin(theta)

  # Add the circle to the plot
  p <- p +
    geom_path(data = data.frame(x = x, y = y), aes(x = x, y = y), color = circle_col)
}

if (show_edges) {
  p <- p +
    geom_path(data = proj_edges, aes(x = x, y = y), color = edges_col)
}

if (show_labels) {
  p <- p +
    geom_text(data = proj_labels, aes(x = x, y = y, label = labels), size = 5*labels_cex)
}

print(p)
# Reshape the data for faceting
hex_data_long <- hex_data %>%
  pivot_longer(cols = starts_with("feature_"),
               names_to = "feature",
               names_pattern = "feature_(.*)_count",
               values_to = "feature_count")


# Select a specific feature to visualize
selected_feature <- col_names[1]

# Create a binary indicator for the presence of the selected feature
hex_data <- hex_data %>%
  mutate(feature_present = sapply(features, function(f) any(f == selected_feature)))

# Create the hexbin plot with counts and feature colors
p <- p + geom_hex(data = hex_data, aes(x = x, y = y, fill = feature_present), stat = "identity") +
  scale_size_continuous(range = c(1, 10)) +
  scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "blue"), name = paste("Presence of", selected_feature)) +
  theme_minimal() +
  labs(title = paste("Hexbin Plot with Counts and Presence of", selected_feature),
       size = "Count",
       fill = paste("Presence of", selected_feature)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())




print(p)


