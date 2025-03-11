// Complete Batch Export Script for Images and Annotations

// Get the current project and its image list
def project = getProject()
if (project == null) {
    println "No project open."
    return
}
def entries = project.getImageList()

// Define export directories for images and annotation CSV files
def exportDir = buildFilePath(PROJECT_BASE_DIR, 'exportedImages')
mkdirs(exportDir)
def annotationDir = buildFilePath(PROJECT_BASE_DIR, 'exportedAnnotations')
mkdirs(annotationDir)

// Loop through each project entry and process each image
for (entry in entries) {
    def imageData = entry.readImageData()
    def server = imageData.getServer()
    def name = GeneralTools.getNameWithoutExtension(entry.getImageName())
    
    // Export full image as a TIFF
    def outFile = buildFilePath(exportDir, name + '.tif')
    writeImage(server, outFile)
    println "Exported image: " + outFile

    // Retrieve annotation objects directly from the image's hierarchy
    def hierarchy = imageData.getHierarchy()
    def annotations = hierarchy.getAnnotationObjects()
    
    // Prepare a CSV file to export annotation details for this image
    def annotationFile = buildFilePath(annotationDir, name + '_annotations.csv')
    def fw = new FileWriter(annotationFile)
    fw.write("RegionName,ROIType,BBox_X,BBox_Y,BBox_Width,BBox_Height\n")
    
    // Loop through each annotation and extract key properties
    annotations.each { annot ->
        def regionName = annot.getName() != null ? annot.getName() : ""
        def roi = annot.getROI()
        def roiType = roi.getClass().getSimpleName()
        // Use individual methods to get the bounding box coordinates
        def x = roi.getBoundsX()
        def y = roi.getBoundsY()
        def w = roi.getBoundsWidth()
        def h = roi.getBoundsHeight()
        fw.write("${regionName},${roiType},${x},${y},${w},${h}\n")
    }
    fw.close()
    println "Exported annotations: " + annotationFile
}
