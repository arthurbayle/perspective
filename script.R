# =============================================================================
# SCRIPT FOR "Satellite remote sensing of alpine vegetation dynamics
#             : challenges and perspectives"
# Purpose: Compute the Simpson Diversity Index of alpine vegetation communities
#          across a sequence of spatial resolutions, for Figure 3 panel b.
# =============================================================================

# Set working directory
setwd("...")

# Load required libraries
library(terra)
library(vectormetrics)
library(abdiv)
library(sf)

# -----------------------------------------------------------------------------
# Load vegetation map of Niwot Ridge
# -----------------------------------------------------------------------------

# Load the vegetation shapefile
NR = vect("...")

# -----------------------------------------------------------------------------
# Define spatial resolutions and pixel selection criterion
# -----------------------------------------------------------------------------

# Sequence of pixel sizes, from 500 m to 5 m, by 5 m decrements
SEQ = seq(from=500, to=5, by=-5)

# Reclassification matrix used to select pixels based on vegetation cover
# Values between 0 and 0.9 are set to NA
# Values between 0.9 and 2 are set to 1
m <- c(0, 0.9, NA,
       0.9, 2, 1)

rclmat <- matrix(m, ncol=3, byrow=TRUE)

# -----------------------------------------------------------------------------
# Compute Simpson Diversity Index for each spatial resolution
# -----------------------------------------------------------------------------

# Create an empty list to store Simpson Diversity Index values
LIST = list()

# Loop over each pixel size
for(i in 1:length(SEQ)){
  
  # ---------------------------------------------------------------------------
  # Prepare grid for current spatial resolution
  # ---------------------------------------------------------------------------
  
  # Current pixel size
  SIZE = SEQ[i]
  
  # Create an empty raster using the current spatial resolution
  RAST = rast(NR, res=c(SIZE, SIZE))
  
  # Rasterize the vegetation polygons and calculate cover proportion per pixel
  RAST = rasterize(NR, RAST, cover=T)
  
  # Apply pixel selection criterion
  RAST = classify(RAST, rclmat)
  
  # Convert selected raster pixels back to polygons
  V = as.polygons(RAST, aggregate=F)
  
  # Add an empty column to store the Simpson Diversity Index
  V$Simpson = NA
  
  # ---------------------------------------------------------------------------
  # Compute Simpson Diversity Index for each selected pixel
  # ---------------------------------------------------------------------------
  
  for(y in 1:length(V)){
    
    # Print progress
    print(paste0("Res: ", SEQ[i],
                 " - Computing Simpson Diversity Index ...",
                 round((y/length(V))*100, 2), " %"))
    
    # Extract current pixel
    TMP = V[y]
    
    # Crop vegetation map to the current pixel
    NRtmp = crop(NR, TMP)
    
    # Optional plot for visual inspection
    # plot(NRtmp, col=rainbow(length(unique(NRtmp$VEG_ID))))
    
    # Compute area of each vegetation class within the current pixel
    AREA = vm_c_ca(st_as_sf(NRtmp), "VEG_ID")
    
    # Format class-area output as a data frame
    DF = data.frame(CLASS = AREA$class, VALUE = AREA$value)
    
    # Convert vegetation-class areas to percentages
    DF$PROP = (DF$VALUE/sum(DF$VALUE))*100
    
    # Compute Simpson Diversity Index for the current pixel
    V$Simpson[y] = simpson(DF$PROP)
    
  }
  
  # Store Simpson Diversity Index values for the current spatial resolution
  LIST[[i]] = V$Simpson
  
  # Name list element using the corresponding pixel size
  names(LIST)[[i]] = SEQ[i]
  
}

# -----------------------------------------------------------------------------
# Plot and save results
# -----------------------------------------------------------------------------

# Plot Simpson Diversity Index distributions across spatial resolutions
boxplot(LIST, outline=F)

# Save the full R workspace
save.image("...")
