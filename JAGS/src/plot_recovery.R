#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This file contains functions to visualize parameter recovery results
#
# FUNCTIONS IN THIS FILE:
# -----------------------
# 1. plot_recovery_faceted()
#    Creates a 2x2 faceted scatterplot of true vs estimated parameters,
#    split by sample size (I) and number of trials (J).
#
# 2. save_recovery_plot()
#    Saves a recovery plot to file in multiple formats (png, pdf, eps).
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# plot_recovery_faceted: Create recovery scatterplots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
plot_recovery_faceted <- function(data, true_col, est_col, main_title, 
                                   xlab, ylab, I_vals, J_vals,
                                   by_condition = TRUE,
                                   col_cond1 = "goldenrod", 
                                   col_cond2 = "steelblue",
                                   use_transparency = TRUE) {
    
    # Create 2x2 layout for I x J facets
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(2, 0, 2, 0))
    
    # Set colors based on whether transparency is supported
    if (use_transparency) {
        pt_col1 <- adjustcolor(col_cond1, 0.35)
        pt_col2 <- adjustcolor(col_cond2, 0.35)
        pt_col_single <- adjustcolor(col_cond1, 0.4)
    } else {
        pt_col1 <- col_cond1
        pt_col2 <- col_cond2
        pt_col_single <- col_cond1
    }
    
    for (i_val in I_vals) {
        for (j_val in J_vals) {
            # Subset data for this facet
            idx <- data$I == i_val & data$J == j_val
            sub_data <- data[idx, ]
            
            # Set up plot
            plot(sub_data[[true_col]], sub_data[[est_col]], 
                 type = "n",
                 xlab = xlab, ylab = ylab,
                 main = paste0("I = ", i_val, ", J = ", j_val))
            
            # Add identity line
            abline(0, 1, lty = 2, col = "gray40")
            
            # Add points
            if (by_condition && "condition" %in% names(sub_data)) {
                idx1 <- sub_data$condition == 0
                idx2 <- sub_data$condition == 1
                points(sub_data[[true_col]][idx1], sub_data[[est_col]][idx1], 
                       col = pt_col1, pch = 16, cex = 0.8)
                points(sub_data[[true_col]][idx2], sub_data[[est_col]][idx2], 
                       col = pt_col2, pch = 16, cex = 0.8)
            } else {
                points(sub_data[[true_col]], sub_data[[est_col]], 
                       col = pt_col_single, pch = 16, cex = 1)
            }
        }
    }
    
    # Add overall title
    mtext(main_title, outer = TRUE, cex = 1.2, line = 0.5)
    
    # Add legend at bottom
    if (by_condition) {
        mtext("Condition 1 = goldenrod, Condition 2 = steelblue", 
              outer = TRUE, side = 1, cex = 0.9, line = 0.5)
    }
    
    # Reset layout
    par(mfrow = c(1, 1))
    
    invisible(NULL)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# save_recovery_plot: Save recovery plot to eps, png, and pdf formats
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
save_recovery_plot <- function(filename_base, output_dir, data, true_col, est_col,
                                main_title, xlab, ylab, I_vals, J_vals,
                                by_condition = TRUE, width = 10, height = 8) {
    
    saved_files <- character(3)
    
    for (i in 1:3) {
        ext <- c("png", "pdf", "eps")[i]
        filepath <- file.path(output_dir, paste0(filename_base, ".", ext))
        
        # EPS does not support transparency
        use_trans <- (ext != "eps")
        
        if (ext == "png") {
            png(filepath, width = width, height = height, units = "in", res = 300)
        } else if (ext == "pdf") {
            pdf(filepath, width = width, height = height)
        } else {
            setEPS()
            postscript(filepath, width = width, height = height)
        }
        
        plot_recovery_faceted(data, true_col, est_col, main_title, 
                               xlab, ylab, I_vals, J_vals, by_condition,
                               use_transparency = use_trans)
        dev.off()
        
        saved_files[i] <- filepath
    }
    
    invisible(saved_files)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# plot_recovery_population: Create population parameter recovery plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
plot_recovery_population <- function(data, params, param_colors, main_title,
                                      I_vals, J_vals, use_transparency = TRUE) {
    
    # Subset to requested parameters
    data_sub <- data[data$param %in% params, ]
    
    # Create 2x2 layout
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(2, 0, 2, 0))
    
    for (i_val in I_vals) {
        for (j_val in J_vals) {
            # Subset data for this facet
            idx <- data_sub$I == i_val & data_sub$J == j_val
            sub_data <- data_sub[idx, ]
            
            # Set up empty plot
            plot(sub_data$true_value, sub_data$estimate, type = "n",
                 xlab = "True Value", ylab = "Estimated Value",
                 main = paste0("I = ", i_val, ", J = ", j_val))
            
            # Add identity line
            abline(0, 1, lty = 2, col = "gray50")
            
            # Add points for each parameter
            for (p_idx in seq_along(params)) {
                p_name <- params[p_idx]
                p_data <- sub_data[sub_data$param == p_name, ]
                
                if (use_transparency) {
                    pt_col <- adjustcolor(param_colors[p_idx], 0.7)
                } else {
                    pt_col <- param_colors[p_idx]
                }
                
                points(p_data$true_value, p_data$estimate, 
                       col = pt_col, pch = 16, cex = 1.5)
            }
        }
    }
    
    # Add overall title
    mtext(main_title, outer = TRUE, cex = 1.2, line = 0.5)
    
    # Add legend at bottom
    param_labels <- gsub("_", " ", params)  # Make labels more readable
    mtext(paste(param_labels, collapse = " | "), 
          outer = TRUE, side = 1, cex = 0.8, line = 0.5)
    
    # Reset layout
    par(mfrow = c(1, 1))
    
    invisible(NULL)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# save_recovery_population: Save population recovery plot to eps, png, and pdf
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
save_recovery_population <- function(filename_base, output_dir, data, 
                                      params, param_colors, main_title,
                                      I_vals, J_vals, width = 10, height = 8) {
    
    saved_files <- character(3)
    
    for (i in 1:3) {
        ext <- c("png", "pdf", "eps")[i]
        filepath <- file.path(output_dir, paste0(filename_base, ".", ext))
        
        # EPS does not support transparency
        use_trans <- (ext != "eps")
        
        if (ext == "png") {
            png(filepath, width = width, height = height, units = "in", res = 300)
        } else if (ext == "pdf") {
            pdf(filepath, width = width, height = height)
        } else {
            setEPS()
            postscript(filepath, width = width, height = height)
        }
        
        plot_recovery_population(data, params, param_colors, main_title,
                                  I_vals, J_vals, use_transparency = use_trans)
        dev.off()
        
        saved_files[i] <- filepath
    }
    
    invisible(saved_files)
}
