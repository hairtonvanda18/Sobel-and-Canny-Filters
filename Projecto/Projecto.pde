import java.util.Arrays;
PImage img,img2;
float rotx = PI / 4;
float roty = PI / 4;
double[][] gradMat, thetaMat;
double[][] strength;
//Kernel para vertical
float[][] Gx = {
    {-1, 0,1},
    {-2,0,2},
    {-1,0,1}
};

//Kernel para horizontal                
float[][] Gy = {
    {1,2,1},
    {0,0,0},
    {-1,-2,-1}
};
//Kernel para vertical
float[][] Gx1 = {
    {1, 0,-1},
    {1,0,-1},
    {1,0,-1}
};

//Kernel para horizontal                
float[][] Gy1 = {
    {1,1,1},
    {0,0,0},
    {-1,-1,-1}
};
double[][] kernel = gaussianTerms(3, 1);
void setup() {
    size(640, 360,P3D);
    
    img = loadImage("sample.png");
    img2 = loadImage("prewitt.png");
    gradMat = new double[img.height][img.width];
    thetaMat = new double[img.height][img.width];
    strength = new double[img.height][img.width];
    textureMode(NORMAL);
    fill(255);
    stroke(color(44, 48, 32));
    //testFilterAccuracy(img2,prewitt(img));
}

void draw() {
    background(255);
    noStroke();
    translate(width / 2.0, height / 2.0, -100);
    rotateX(rotx);
    rotateY(roty);
    scale(90);
    TexturedCube();
    
    //image(canny(img,0.2,0.4),0,0);
}
PImage canny(PImage img,float min,float max) {
    return doubleThresholding(nonMax(sobel(gaussianBlur(img))), min, max);
}
PImage gaussianBlur(PImage img) {
  img.loadPixels();
    // Create an opaque image of the same size as the original
    PImage edgeImg = createImage(img.width, img.height, RGB);

    // Loop through every pixel in the image
    for (int y = 1; y < img.height - 1; y++) { // Skip top and bottom edges
        for (int x = 1; x < img.width - 1; x++) { // Skip left and right edges
            float sum = 0; // Kernel sum for this pixel
            for (int ky = -1; ky <= 1; ky++) {
                for (int kx = -1; kx <= 1; kx++) {
                    // Calculate the adjacent pixel for this kernel point
                    int pos = (y + ky) * img.width + (x + kx);
                    // Image is grayscale, red/green/blue are identical
                    float val = red(img.pixels[pos]);
                    // Multiply adjacent pixels based on the kernel values
                    sum += kernel[ky + 1][kx + 1] * val;
                }
            }
            // For this pixel in the new image, set the gray value
            // based on the sum from the kernel
            edgeImg.pixels[y * img.width + x] = color(sum);
        }
    }
    // State that there are changes to edgeImg.pixels[]
    edgeImg.updatePixels();
    return edgeImg;
}
double gaussian(int x, double sigma) {
    double c = 2.0 * sigma * sigma;
    return Math.exp(-x * x / c) / Math.sqrt(c * Math.PI);
}
double[][] gaussianTerms(int kernelSize, double sigma) {
    double[][] terms = new double[kernelSize][kernelSize];
    for (int i = 0; i < kernelSize; ++i) {
        for (int j = 0; j < kernelSize; ++j) {
            terms[i][j] = gaussian(j - kernelSize / 2, sigma);
        }
    }
    return terms;
}
PImage sobel(PImage img) {
    img.loadPixels();
    //Criação de uma imagem com as mesmas dimensões da original
    PImage imgFiltrada = createImage(img.width, img.height, RGB);
    for (int y = 1; y < img.height - 1; y++) { // Imagem original altura 
        for (int x = 1; x < img.width - 1; x++) { // Imagem original largura 
            float sumy = 0;
            float sumx = 0;
            for (int ky = -1; ky <= 1; ky++) {
                for (int kx = -1; kx <= 1; kx++) {
                    // Calcule o pixel adjacente para este ponto do kernel
                    int pos = (y + ky) * img.width + (x + kx);
                    // pegar o valor de 0-255 , podia ser blue ou green, mas a imagem é em tons de cinza, então é irrelevante
                    float val = red(img.pixels[pos]);
                    // Multiplique os pixels adjacentes com base nos valores do kernel para X e Y
                    sumy += Gy[ky + 1][kx + 1] * val;
                    sumx += Gx[ky + 1][kx + 1] * val;
                }
            }
            // Para este pixel na nova imagem, defina o valor de cinza
            // produto vectorial dos dois kernels
            gradMat[y][x] = Math.min(255, Math.sqrt(Math.pow(sumy, 2) + Math.pow(sumx, 2)));
            thetaMat[y][x] = Math.atan2(sumy, sumx);
            strength[y][x] = 0;
            //Inserção do novo pixel

            imgFiltrada.pixels[y * img.width + x] = color(Math.round(gradMat[y][x]));
            //imgFiltrada.pixels[y*img.width + x] = color(sumy);
        }
    }


    // actualizar a imagem filtrada
    imgFiltrada.updatePixels();
    return imgFiltrada;
}
PImage prewitt(PImage img) {
    img.loadPixels();
    //Criação de uma imagem com as mesmas dimensões da original
    PImage imgFiltrada = createImage(img.width, img.height, RGB);


    for (int y = 1; y < img.height - 1; y++) { // Imagem original altura 
        for (int x = 1; x < img.width - 1; x++) { // Imagem original largura 
            float sumy = 0;
            float sumx = 0;
            for (int ky = -1; ky <= 1; ky++) {
                for (int kx = -1; kx <= 1; kx++) {
                    // Calcule o pixel adjacente para este ponto do kernel
                    int pos = (y + ky) * img.width + (x + kx);
                    // pegar o valor de 0-255 , podia ser blue ou green, mas a imagem é em tons de cinza, então é irrelevante
                    float val = red(img.pixels[pos]);
                    // Multiplique os pixels adjacentes com base nos valores do kernel para X e Y
                    sumy += Gy1[ky + 1][kx + 1] * val;
                    sumx += Gx1[ky + 1][kx + 1] * val;
                }
            }
            // Para este pixel na nova imagem, defina o valor de cinza
            // produto vectorial dos dois kernels
            gradMat[y][x] = Math.min(255, Math.sqrt(Math.pow(sumy, 2) + Math.pow(sumx, 2)));
            thetaMat[y][x] = Math.atan2(sumy, sumx);
            //Inserção do novo pixel

            imgFiltrada.pixels[y * img.width + x] = color(Math.round(gradMat[y][x]));
            //imgFiltrada.pixels[y*img.width + x] = color(sumy);
        }
    }


    // actualizar a imagem filtrada
    imgFiltrada.updatePixels();
    return imgFiltrada;
}
PImage nonMax(PImage img) {
    img.loadPixels();
    int count = 0;

    for (int i = 1; i < img.height - 1; i++) {
        for (int j = 1; j < img.width - 1; j++) {
            double maxValue = gradMat[i][j];
            if (thetaMat[i][j] == 0) {
                if (j > 0) {
                    if (maxValue < gradMat[i][j - 1]) {
                        maxValue = gradMat[i][j - 1];
                    }
                } else if (j < img.width - 1) {
                    if (maxValue < gradMat[i][j + 1]) {
                        maxValue = gradMat[i][j + 1];
                    }
                }
            } else if (thetaMat[i][j] == 45) {
                if (j < img.width - 1 && i > 0) {
                    if (maxValue < gradMat[i - 1][j + 1]) {
                        maxValue = gradMat[i - 1][j + 1];
                    }
                } else if (j > 0 && i < img.height - 1) {
                    if (maxValue < gradMat[i + 1][j - 1]) {
                        maxValue = gradMat[i + 1][j - 1];
                    }
                }
            } else if (thetaMat[i][j] == 90) {
                if (j > 0) {
                    if (maxValue < gradMat[i - 1][j]) {
                        maxValue = gradMat[i - 1][j];
                    }
                } else if (j < img.height - 1) {
                    if (maxValue < gradMat[i + 1][j]) {
                        maxValue = gradMat[i + 1][j];
                    }
                }
            } else if (thetaMat[i][j] == 135) {
                if (j > 0 && i > 0) {
                    if (maxValue < gradMat[i - 1][j - 1]) {
                        maxValue = gradMat[i - 1][j - 1];
                    }
                } else if (j < img.width - 1 && i < img.height - 1) {
                    if (maxValue < gradMat[i + 1][j + 1]) {
                        maxValue = gradMat[i + 1][j + 1];
                    }
                }
            }
            if (maxValue != gradMat[i][j]) {
                img.pixels[count] = color(0);
                gradMat[i][j] = 0;
            }
            count++;
        }
    }
    img.updatePixels();
    return img;
}
color[][] getPixelMatrix(PImage img) {
    img.loadPixels();
    color[][] pixelMatrix = new color[img.height][img.width];

    int line = 0;
    int col = 0;

    for (int i = 0; i < img.width * img.height; i++) {
        if (col == img.width) {
            col = 0;
            line++;
        }
        pixelMatrix[line][col] = img.pixels[i];
        col++;
    }
    return pixelMatrix;
}
PImage doubleThresholding(PImage img, float minThresh, float maxThresh) {
    img.loadPixels();
    int count = 0;
    color[][] pixels = getPixelMatrix(img);
    for (int i = 0; i < img.height; i++) {
        for (int j = 0; j < img.width; j++) {
            double intensity = gradMat[i][j] / 255.0;
            if (intensity < minThresh) {
                img.pixels[count] = color(0);
                strength[i][j] = 0;
            } else if (intensity > maxThresh) {
                strength[i][j] = 2;
                img.pixels[count] = color(255);
            } else {
                strength[i][j] = 1;
            }
            count++;
        }
    }
    img.updatePixels();
    //Applying Hysteresis
    img.loadPixels();
    pixels = getPixelMatrix(img);
    count = 0;
    for (int i = 0; i < img.height; i++) {
        for (int j = 0; j < img.width; j++) {

            if (strength[i][j] != 1) {
                count++;
                continue;
            }
            boolean connected = false;
            //Checking for connections to strong neighbours
            if (j > 0)
                connected = connected || (strength[i][j - 1] == 2);
            if (j < img.width - 1)
                connected = connected || (strength[i][j + 1] == 2);
            if (j < img.width - 1 && i > 0)
                connected = connected || (strength[i - 1][j + 1] == 2);
            if (j > 0 && i < img.height - 1)
                connected = connected || (strength[i + 1][j - 1] == 2);
            if (i > 0)
                connected = connected || (strength[i - 1][j] == 2);
            if (i < img.height - 1)
                connected = connected || (strength[i + 1][j] == 2);
            if (j > 0 && i > 0)
                connected = connected || (strength[i - 1][j - 1] == 2);
            if (j < img.width - 1 && i < img.height - 1)
                connected = connected || (strength[i + 1][j + 1] == 2);
            if (!connected) {
                img.pixels[count] = color(0);
            } else
                img.pixels[count] = color(255);
            count++;
        }
    }
    img.updatePixels();
    return img;
}
void testFilterAccuracy(PImage refImg, PImage filterImg){
    refImg.filter(GRAY);
    filterImg.filter(GRAY);
    refImg.loadPixels();
    filterImg.loadPixels();

    float refSize = refImg.width * refImg.height;
    float filterSize = filterImg.width * filterImg.height;
    assert refSize == filterSize;

    float error = 0;
    for(int i = 0; i < refSize; i++){
        color refPix = refImg.pixels[i];
        color filterPix = filterImg.pixels[i];
        error += abs(red(filterPix) - red(refPix));
    }

    print("Error médio absoluto por pixel: ");
    print(error / refSize);
}
void TexturedCube() {
    beginShape(QUADS);
    texture(img);


    // +Z "front" face
    vertex(-1, -1, 1, 0, 0);
    vertex(1, -1, 1, 1, 0);
    vertex(1, 1, 1, 1, 1);
    vertex(-1, 1, 1, 0, 1);
    endShape();
    beginShape(QUADS);
    texture(canny(img,0.2,0.4));
    // -Z "back" face
    vertex(1, -1, -1, 0, 0);
    vertex(-1, -1, -1, 1, 0);
    vertex(-1, 1, -1, 1, 1);
    vertex(1, 1, -1, 0, 1);
    endShape();
    beginShape(QUADS);
    texture(canny(img,0.1,0.3));
    // +Y "bottom" face
    vertex(-1, 1, 1, 0, 0);
    vertex(1, 1, 1, 1, 0);
    vertex(1, 1, -1, 1, 1);
    vertex(-1, 1, -1, 0, 1);
    endShape();
    beginShape(QUADS);
    texture(prewitt(img));
    // -Y "top" face
    vertex(-1, -1, -1, 0, 0);
    vertex(1, -1, -1, 1, 0);
    vertex(1, -1, 1, 1, 1);
    vertex(-1, -1, 1, 0, 1);
    endShape();
    beginShape(QUADS);
    texture(gaussianBlur(img));
    // +X "right" face
    vertex(1, -1, 1, 0, 0);
    vertex(1, -1, -1, 1, 0);
    vertex(1, 1, -1, 1, 1);
    vertex(1, 1, 1, 0, 1);
    endShape();
    beginShape(QUADS);
    texture(sobel(img));
    // -X "left" face
    vertex(-1, -1, -1, 0, 0);
    vertex(-1, -1, 1, 1, 0);
    vertex(-1, 1, 1, 1, 1);
    vertex(-1, 1, -1, 0, 1);

    endShape();
}

void mouseDragged() {
    float rate = 0.01;
    rotx += (pmouseY - mouseY) * rate;
    roty += (mouseX - pmouseX) * rate;
}
