# R (version 4.4.3)


library(jsonlite)
library(sf)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(data.table)
library(pbapply)
library(future)
plan("multisession", workers = 16)


# SET PATHS
wdir <- "" # set your working directory HERE
data_dir <- file.path(wdir, 'json') # this folder contains the geojson files for the 132 samples, as classified by QuPath classifier
clin_path <- "" # set path to clinical annotation file HERE

out_dir <- file.path(wdir, "polygon_analysis")
if(!dir.exists(out_dir)) dir.create(out_dir)

samples <- sapply(dir(data_dir, pattern = ".geojson", full.names = FALSE), function(x) gsub(".geojson", "", x))

num_threads <- 16

# test different buffer sizes to find the one that best captures the tumor morphology without over-smoothing
# In Caire et al., a buffer of 10 pixels was used
buffer_test <- c(10, 50)

metric_list_all <- list()

for(bt in buffer_test){
    metrics_list <- pblapply(
        samples,
        function(s) {

            all.cells <- jsonlite::fromJSON(file.path(data_dir, paste0(s, ".geojson")))$features 
            all.cell.type <- data.frame(cellID = all.cells$id, all.cells$properties$classification)
            all.cell.type <- all.cell.type[1:dim(all.cell.type)[1],]
            json.obj <- st_read(file.path(data_dir,paste0(s, ".geojson")))
            json.obj <- st_set_crs(json.obj, NA)
            json.obj <- json.obj[json.obj$objectType!="annotation","geometry"]
            json.obj <- cbind(json.obj, "cellID" = all.cell.type$cellID)
            json.obj <- cbind(json.obj, "Type" = all.cell.type$name)
            json.obj <- json.obj[json.obj$Type == "Tumor",]
            
            bbox <- st_bbox(json.obj)
            w = bbox[3]-bbox[1]
            h = bbox[4]-bbox[2]

            json.obj_simple <- st_make_valid(json.obj[st_is_valid(json.obj),])
            json.obj_simple <- st_union(json.obj_simple)
            json.obj_simple <- st_simplify(json.obj_simple, dTolerance = 5)
            json.obj_simple <- st_buffer(json.obj_simple, dist = bt)
            json.obj_simple <- st_buffer(json.obj_simple, dist = (bt/2))
            json.obj_simple <- st_make_valid(json.obj_simple)

            png(file.path(out_dir, paste0(s, "_polygons_simple_bt", bt, ".png")), width = round(w/10), height = round(h/10))
            p <- ggplot() + 
            geom_sf(data = json.obj_simple, fill = "red") +
            theme_void() +
            theme(
                legend.position = "none",
                panel.background = element_rect(fill = "black", color = NA),
                plot.background = element_rect(fill = "black", color = NA)
            )
            print(p)
            dev.off()

            polys <- st_cast(json.obj_simple, "POLYGON")
            areas <- st_area(polys)
            perims <- st_length(st_cast(polys, "MULTILINESTRING"))
            coords <- st_coordinates(st_cast(polys[i, ], "LINESTRING"))[,1:2]
            })

            metrics <- data.frame(
                sample = s,
                n_cells = nrow(json.obj),
                n_polygons = length(polys),
                polygons_on_area = length(polys)/sum(as.numeric(areas)),
                polygons_on_cells = length(polys)/nrow(json.obj),
                area = sum(areas),
                perimeter = weighted.mean(as.numeric(perims), w = as.numeric(areas), na.rm = TRUE),
                perimeter_on_area = weighted.mean(as.numeric(perims/areas), w = as.numeric(areas), na.rm = TRUE)
            )

        }, 
        cl = num_threads
    )
    metrics_df <- do.call(rbind, metrics_list)
    metric_list_all[[bt]] <- metrics_df
}

metrics_df <- metric_list_all[[10]]
metrics_df$buffer <- 10

metrics_df1 <- metric_list_all[[50]]
metrics_df1$buffer <- 50

metrics_df <- rbind(metrics_df, metrics_df1)

write.csv(metrics_df, file.path(out_dir, "polygon_metrics.csv"), row.names = FALSE)
metrics_df <- fread(file.path(out_dir, "polygon_metrics.csv"))



# load clinical annotation for primaries
clin_df <- fread(clin_path)

clin_df$Dist_rel <- ifelse(clin_df$Dist_rel == 0, "BBG", "BBB")
metrics_df <- left_join(metrics_df, clin_df %>% select(ID_isto, Dist_rel), by = c("sample" = "ID_isto"))


write.csv(metrics_df %>% filter(buffer == 10) %>%
    select(sample, polygons_on_area, perimeter_on_area, Dist_rel),
    file.path(out_dir, "polygon_metrics.csv"), row.names = FALSE)


metrics_plot <- filter(metrics_df, !is.na(Dist_rel)) %>%
    gather(key = "Metric", value = "value", -sample, -Dist_rel, -buffer)

pdf(file.path(out_dir, "polygon_metrics_boxplots.pdf"), useDingbats = FALSE)
for(bt in buffer_test) {
    metrics_plot_bt <- filter(metrics_plot, buffer == bt)
    p <- ggplot(metrics_plot_bt, aes(x = Dist_rel, y = value, fill = Dist_rel)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
        facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
        theme_minimal() +
        theme(
            legend.position = "none",
            strip.text = element_text(size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        xlab("") +
        ylab("Value") +
        ggtitle(paste("Buffer used: ", bt))

        print(p)
}
dev.off()








make_scale_bar_r <- function(x_vals, y_vals, microns_per_pixel = 0.12028) {
  # Adds a scale bar to a ggplot.
    #Parameters:
      #x_vals: vector of x coordinates in pixels
      #y_vals: vector of y coordinates in pixels
      #microns_per_pixel: conversion factor from pixels to microns. Default is conversion factor for commercial CosMx instrument
  # Example usage:
    # scale_bar = make_scale-bar_r(x_vals = cell_meta$CenterX_global_px, y_vals = cell_meta$CenterY_global_px)
    # ggplot() + scale_bar$bg + scale_bar$rect + scale_bar$label
  
  # Calculate x-axis range
  x_range <- range(x_vals, na.rm = TRUE)
  x_length <- diff(x_range)
  x_length_um <- x_length * microns_per_pixel
  # Target scale length ~1/4 of the x-axis
  target <- x_length_um / 4
  # Compute order of magnitude
  order <- 10^floor(log10(target))
  mantissa <- target / order

  # Round mantissa to nearest 1, 2, or 5
  nice_mantissa <- if (mantissa < 1.5) {
    1
  } else if (mantissa < 3.5) {
    2
  } else if (mantissa < 7.5) {
    5
  } else {
    10
  }

  # Final scale length in pixels
  scale_length_um <- nice_mantissa * order
  scale_length_px <- scale_length_um / microns_per_pixel
  # Format label
  scale_label <- if (scale_length_um >= 1000) {
    paste0(scale_length_um / 1000, " mm")
  } else {
    paste0(scale_length_um, " µm")
  }

  # Set coordinates for the scale bar 
  x_start <- x_range[2] - scale_length_px * 1.1
  x_end <- x_range[2] - scale_length_px * 0.1
  y_pos <- min(y_vals, na.rm = TRUE) + scale_length_px * 0.1

  # Generate scale bar background, scale bar and annotation to return
  list(
    bg = annotation_custom(
    grob = rectGrob(gp = gpar(fill = "white", alpha = 0.8, col = NA)),
    xmin = x_start- scale_length_px*0.05, xmax = x_end+ scale_length_px*0.05,
    ymin = y_pos- scale_length_px*0.05, ymax = y_pos + scale_length_px*0.3
  ),
    rect = annotation_custom(
      grob = rectGrob(gp = gpar(fill = "black")),
      xmin = x_start, xmax = x_end,
      ymin = y_pos, ymax = y_pos + scale_length_px * 0.05
    ),
    label = annotation_custom(
      grob = textGrob(scale_label, gp = gpar(col = "black"), just = "center", vjust = 0),
      xmin = (x_start + x_end)/2, xmax = (x_start + x_end)/2,
      ymin = y_pos + scale_length_px * 0.1, ymax = y_pos + scale_length_px * 0.1
    )
  )
}

plot_samples <- c("11-02527-YAP", "16-12905-YAP")


for(s in plot_samples){

    bt <- 10

    all.cells <- jsonlite::fromJSON(file.path(data_dir, paste0(s, ".geojson")))$features 
    all.cell.type <- data.frame(cellID = all.cells$id, all.cells$properties$classification)
    all.cell.type <- all.cell.type[1:dim(all.cell.type)[1],]
    json.obj <- st_read(file.path(data_dir,paste0(s, ".geojson")))
    json.obj <- st_set_crs(json.obj, NA)
    json.obj <- json.obj[json.obj$objectType!="annotation","geometry"]
    json.obj <- cbind(json.obj, "cellID" = all.cell.type$cellID)
    json.obj <- cbind(json.obj, "Type" = all.cell.type$name)
    json.obj <- json.obj[json.obj$Type == "Tumor",]
                
    bbox <- st_bbox(json.obj)
    w = bbox[3]-bbox[1]
    h = bbox[4]-bbox[2]

    json.obj_simple <- st_make_valid(json.obj[st_is_valid(json.obj),])
    json.obj_simple <- st_union(json.obj_simple)
    json.obj_simple <- st_simplify(json.obj_simple, dTolerance = 5)
    json.obj_simple <- st_buffer(json.obj_simple, dist = bt)
    json.obj_simple <- st_buffer(json.obj_simple, dist = (bt/2))
    json.obj_simple <- st_make_valid(json.obj_simple)

    pdf(file.path(out_dir, paste0(s, "_polygons_simple_bt", bt, ".pdf")), width = w/4000, height = h/4000, useDingbats = FALSE)
    p <- ggplot() + 
    geom_sf(data = json.obj_simple, fill = "red", color = "red", linewidth = 0) +
    theme_void() +
    theme(
        legend.position = "none",
        panel.background = element_rect(fill = "black", color = NA),
        plot.background = element_rect(fill = "black", color = NA)
    )
    scale_bar = make_scale_bar_r(
        x_vals = st_coordinates(json.obj)[,1],
        y_vals = st_coordinates(json.obj)[,2],
        microns_per_pixel = 0.2210
    )
    print(p + scale_bar$bg + scale_bar$rect + scale_bar$label)
    dev.off()
    
    png(file.path(out_dir, paste0(s, "_polygons_simple_bt", bt, ".png")), width = w/10, height = h/10)
    p <- ggplot() + 
      geom_sf(data = json.obj_simple, fill = "red", color = "red", linewidth = 0) +
      theme_void() +
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = "black", color = NA),
        plot.background = element_rect(fill = "black", color = NA)
      )
    print(p)
    dev.off()
}
