## Function for plotting vulnerability donuts
## Reads table_mpa and makes 6 donuts of ecological vulnerability and 6 donuts of social vulnerability

  df <- table_mpa
  
  ## Transform tibble to data frame
  df <- as.data.frame (df)
  
  df [df == "Low"] <- "Low (percentile < 20)"
  df [df == "Moderate"] <- "Moderate (percentile 20 - 40)"
  df [df == "High"] <- "High (percentile 40 - 60)"
  df [df == "Very high"] <- "Very high (percentile 60 - 80)"
  df [df == "Extreme"] <- "Extreme (percentile > 80)"
  
  for(i in 1:nrow(df)) # each row is a scenario. makes two plots (one ecolocial, one social)
  {
    
    ## Ecological donut
    eco.AC   <- df [i,"EcologicalAdaptiveCapacity"]
    eco.S    <- df [i,"EcologicalSensitivity"]
    eco.E    <- df [i,"EcologicalExposure"]
    
    eco.Q_AC <- df [i,"EcologicalAdaptiveCapacityQ"]
    eco.Q_S  <- df [i,"EcologicalSensitivityQ"]
    eco.Q_E  <- df [i,"EcologicalExposureQ"]
    
    Index    <- toupper (df [i,"EcologicalVulnerability"])
    Index_short <- substr(Index, start=1, stop= unlist(gregexpr(pattern=" \\(", Index)) -1 )
    
    if (Index == "LOW (PERCENTILE < 20)") col.Index <- "olivedrab4"
    if (Index == "MODERATE (PERCENTILE 20 - 40)") col.Index <- "#0073C2FF"
    if (Index == "HIGH (PERCENTILE 40 - 60)") col.Index <- "#EFC000FF"
    if (Index == "VERY HIGH (PERCENTILE 60 - 80)") col.Index <- "darkorange2"
    if (Index == "EXTREME (PERCENTILE > 80)") col.Index <- "#CD534CFF"
    
    df.plot <- data.frame ("Dimension"=c("AC", "S", "E"), "Value"=c(eco.AC, eco.S, eco.E), "Category"=c(eco.Q_AC, eco.Q_S, eco.Q_E))
    
    ## Color for the categories
    df.plot$Color <- rep(NA, 3)
    df.plot$Color [which (df.plot$Category == "Low (percentile < 20)")] <- "olivedrab4" 
    df.plot$Color [which (df.plot$Category == "Moderate (percentile 20 - 40)")] <- "#0073C2FF" #blue
    df.plot$Color [which (df.plot$Category == "High (percentile 40 - 60)")] <- "#EFC000FF" #yellow
    df.plot$Color [which (df.plot$Category == "Very high (percentile 60 - 80)")] <- "darkorange2"
    df.plot$Color [which (df.plot$Category == "Extreme (percentile > 80)")] <- "#CD534CFF" #red
    
    .plot.title <- paste ("Ecological", paste(mpa, df[i,"Scenario"], sep=" "), sep=" - ")
    
    ##############################################################################
    library(dplyr)
    #breaks       = c( 0, 0.1, 0.2, 0.3, 0.4,0.6)
    breaks       = c( "Low (percentile < 20)", "Moderate (percentile 20 - 40)", "High (percentile 40 - 60)", "Very high (percentile 60 - 80)", "Extreme (percentile > 80)")
    df.plot$bins <- match (df.plot$Category, breaks)
    Colors <- c("olivedrab4", "#0073C2FF", "#EFC000FF", "darkorange2", "#CD534CFF")
    df.plot$bins <- factor(df.plot$Category, levels=breaks)
    df.plot$Category_short <- substr(df.plot$Category, start=1, stop= unlist(gregexpr(pattern=" \\(", df.plot$Category)) )
    ##############################################################################
    
    ## The plot
    library(ggplot2)
    
    plot.name <- paste (paste("figures/", .plot.title, sep=""), ".png", sep="")
    png (here::here(plot.name), width=2400, height=2400, unit="px", res=350)
    
    
    donut <- ggplot (data=df.plot, aes(x=2, y=Value, fill=bins)) + ## values of the donut
      geom_bar(stat = "identity", color = "white") +
      coord_polar(theta = "y", start = 0) +
      geom_text(aes(label = Dimension), position=position_stack(vjust=0.5), color = "white", fontface="bold", size=8) +
      geom_text(data=df.plot, aes(x=0.6, y=0.0, label="Vulnerability"), inherit.aes = FALSE, color="black",   fontface="bold", size=6) +
      geom_text(data=df.plot, aes(x=0.1, y=0.0, label=Index_short)          , inherit.aes = FALSE, color=col.Index, fontface="bold", size=10) +
      scale_fill_manual("Categories", values = Colors, drop = FALSE) + 
      theme_void() + ## clean the background
      xlim(0.1, 2.5) +
      ggtitle(.plot.title) +
      theme (plot.title = element_text(hjust=0.5, vjust=-4, size = 15, face="bold"), 
             legend.position="right",
             legend.text = element_text(size=12),
             legend.title = element_text(size=16))
    print (donut)
    
    
    dev.off() 
    
    ## Social donut
    social.AC   <- df [i,"SocialAdaptiveCapacity"]
    social.S    <- df [i,"SocialSensitivity"]
    social.E    <- df [i,"SocialExposure"]
    
    social.Q_AC <- df [i,"SocialAdaptiveCapacityQ"]
    social.Q_S  <- df [i,"SocialSensitivityQ"]
    social.Q_E  <- df [i,"SocialExposureQ"]
    
    Index    <- toupper (df [i,"SocialVulnerability"])
    Index_short <- substr(Index, start=1, stop= unlist(gregexpr(pattern=" \\(", Index)) -1 )
    
    if (Index == "LOW (PERCENTILE < 20)") col.Index <- "olivedrab4"
    if (Index == "MODERATE (PERCENTILE 20 - 40)") col.Index <- "#0073C2FF"
    if (Index == "HIGH (PERCENTILE 40 - 60)") col.Index <- "#EFC000FF"
    if (Index == "VERY HIGH (PERCENTILE 60 - 80)") col.Index <- "darkorange2"
    if (Index == "EXTREME (PERCENTILE > 80)") col.Index <- "#CD534CFF"
    
    df.plot <- data.frame ("Dimension"=c("AC", "S", "E"), "Value"=c(social.AC, social.S, social.E), "Category"=c(social.Q_AC, social.Q_S, social.Q_E))
    
    ## Color for the categories
    df.plot$Color <- rep(NA, 3)
    df.plot$Color [which (df.plot$Category == "Low (percentile < 20)")] <- "olivedrab4" 
    df.plot$Color [which (df.plot$Category == "Moderate (percentile 20 - 40)")] <- "#0073C2FF" #blue
    df.plot$Color [which (df.plot$Category == "High (percentile 40 - 60)")] <- "#EFC000FF" #yellow
    df.plot$Color [which (df.plot$Category == "Very high (percentile 60 - 80)")] <- "darkorange2"
    df.plot$Color [which (df.plot$Category == "Extreme (percentile > 80)")] <- "#CD534CFF" #red
    
    .plot.title <- paste ("Social", paste(mpa, df[i,"Scenario"], sep=" "), sep=" - ")
    
    ##############################################################################
    #breaks       = c( 0, 0.1, 0.2, 0.3, 0.4,0.6)
    breaks       = c( "Low (percentile < 20)", "Moderate (percentile 20 - 40)", "High (percentile 40 - 60)", "Very high (percentile 60 - 80)", "Extreme (percentile > 80)")
    df.plot$bins <- match (df.plot$Category, breaks)
    Colors <- c("olivedrab4", "#0073C2FF", "#EFC000FF", "darkorange2", "#CD534CFF")
    df.plot$bins <- factor(df.plot$Category, levels=breaks)
    df.plot$Category_short <- substr(df.plot$Category, start=1, stop= unlist(gregexpr(pattern=" \\(", df.plot$Category)) )
    ##############################################################################
    
    ## The plot
    plot.name <- paste (paste("figures/", .plot.title, sep=""), ".png", sep="")
    png (here::here(plot.name), width=2400, height=2400, unit="px", res=350)
    
    
    donut <- ggplot (data=df.plot, aes(x=2, y=Value, fill=bins)) + ## values of the donut
      geom_bar(stat = "identity", color = "white") +
      coord_polar(theta = "y", start = 0) +
      geom_text(aes(label = Dimension), position=position_stack(vjust=0.5), color = "white", fontface="bold", size=8) +
      geom_text(data=df.plot, aes(x=0.6, y=0.0, label="Vulnerability"), inherit.aes = FALSE, color="black",   fontface="bold", size=6) +
      geom_text(data=df.plot, aes(x=0.1, y=0.0, label=Index_short)          , inherit.aes = FALSE, color=col.Index, fontface="bold", size=10) +
      scale_fill_manual("Categories", values = Colors, drop = FALSE) + 
      theme_void() + ## clean the background
      xlim(0.1, 2.5) +
      ggtitle(.plot.title) +
      theme (plot.title = element_text(hjust=0.5, vjust=-4, size = 15, face="bold"), 
             legend.position="right",
             legend.text = element_text(size=12),
             legend.title = element_text(size=16))
    print (donut)
    
    
    dev.off() 
    
  }
