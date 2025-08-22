import qupath.lib.roi.ROIs
import qupath.lib.objects.PathAnnotationObject
import qupath.lib.regions.ImagePlane

// Parámetros
double spotDiameter = 55.0     // en micrómetros
double spotRadius = spotDiameter / 2.0
double spacing = 100.0  // sin superposición

// Obtener información de la imagen
def imageData = getCurrentImageData()
def server = imageData.getServer()
def cal = server.getPixelCalibration()

// Obtener tamaño en µm
double imageWidth = server.getWidth() * cal.pixelWidth
double imageHeight = server.getHeight() * cal.pixelHeight

// Crear plano de imagen 2D (z=0, t=0)
def plane = ImagePlane.getDefaultPlane()

// Crear anotaciones circulares
def annotations = []

for (double y = spotRadius; y < imageHeight; y += spacing) {
    for (double x = spotRadius; x < imageWidth; x += spacing) {

        // Convertir a píxeles
        double px = x / cal.pixelWidth
        double py = y / cal.pixelHeight
        double pw = spotDiameter / cal.pixelWidth
        double ph = spotDiameter / cal.pixelHeight

        def roi = ROIs.createEllipseROI(
            px - pw / 2.0, py - ph / 2.0, pw, ph, plane
        )
        def annotation = new PathAnnotationObject(roi)
        annotations << annotation
    }
}

// Agregar las anotaciones al proyecto
addObjects(annotations)
print "Círculos creados: " + annotations.size()
