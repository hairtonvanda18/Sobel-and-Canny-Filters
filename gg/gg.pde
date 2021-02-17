PImage img;
float rotx = PI/4;
float roty = PI/4;
double gradMat,thetaMat;
//Kernel para vertical
float[][] Gx = {{ -1, 0, 1}, 
                { -2,  0, 2}, 
                { -1, 0, 1}};
                    
//Kernel para horizontal                
float[][] Gy = {{ 1, 2, 1}, 
                { 0,  0, 0}, 
                { -1, -2, -1}};
double[][] kernel = gaussianTerms(3,1);
void setup() {
  size(640, 360, P3D);
  img = loadImage("CC.PNG");
  textureMode(NORMAL);
  fill(255);
  stroke(color(44,48,32));
}

void draw() {
  background(255);
  noStroke();
  translate(width/2.0, height/2.0, -100);
  rotateX(rotx);
  rotateY(roty);
  scale(90);
  TexturedCube();
}
PImage canny(PImage img){
 return sobel(gaussianBlur(img)); 
}
PImage gaussianBlur(PImage img){
  img.loadPixels();

  // Create an opaque image of the same size as the original
  PImage edgeImg = createImage(img.width, img.height, RGB);

  // Loop through every pixel in the image
  for (int y = 1; y < img.height-1; y++) {   // Skip top and bottom edges
    for (int x = 1; x < img.width-1; x++) {  // Skip left and right edges
      float sum = 0; // Kernel sum for this pixel
      for (int ky = -1; ky <= 1; ky++) {
        for (int kx = -1; kx <= 1; kx++) {
          // Calculate the adjacent pixel for this kernel point
          int pos = (y + ky)*img.width + (x + kx);
          // Image is grayscale, red/green/blue are identical
          float val = red(img.pixels[pos]);
          // Multiply adjacent pixels based on the kernel values
          sum += kernel[ky+1][kx+1] * val;
        }
      }
      // For this pixel in the new image, set the gray value
      // based on the sum from the kernel
      edgeImg.pixels[y*img.width + x] = color(sum);
    }
  }
  // State that there are changes to edgeImg.pixels[]
  edgeImg.updatePixels();
  return edgeImg;
}
double gaussian(int x,double sigma){
  double c = 2.0 * sigma*sigma;
  return Math.exp(-x *x / c)/ Math.sqrt(c*Math.PI);
}
double[][] gaussianTerms(int kernelSize,double sigma){
  double[][] terms = new double[kernelSize][kernelSize];
  for(int i = 0;i < kernelSize;++i){
    for(int j=0;j<kernelSize;++j){
      terms[i][j] = gaussian(j-kernelSize/2,sigma);
    }
  }
  return terms;
}

PImage sobel(PImage img){
  img.loadPixels();
  //Criação de uma imagem com as mesmas dimensões da original
  PImage imgFiltrada = createImage(img.width, img.height, RGB);
  
  
   for (int y = 1; y < img.height-1; y++) {// Imagem original altura 
    for (int x = 1; x < img.width-1; x++) { // Imagem original largura 
      float sumy = 0; 
      float sumx=0;
      for (int ky = -1; ky <= 1; ky++) {
        for (int kx = -1; kx <= 1; kx++) {
          // Calcule o pixel adjacente para este ponto do kernel
          int pos = (y + ky)*img.width + (x + kx);
          // pegar o valor de 0-255 , podia ser blue ou green, mas a imagem é em tons de cinza, então é irrelevante
          float val = red(img.pixels[pos]);
          // Multiplique os pixels adjacentes com base nos valores do kernel para X e Y
          sumy += Gy[ky+1][kx+1] * val;
          sumx += Gx[ky+1][kx+1] * val;
        }
      }
      // Para este pixel na nova imagem, defina o valor de cinza
      // produto vectorial dos dois kernels
      gradMat=Math.min(255,Math.sqrt(Math.pow(sumy,2)+Math.pow(sumx,2)));
      thetaMat=Math.atan2(sumy,sumx);
      //Inserção do novo pixel
      
      imgFiltrada.pixels[y*img.width + x] = color(Math.round(gradMat));
      //imgFiltrada.pixels[y*img.width + x] = color(sumy);
    }
  }
  
  
// actualizar a imagem filtrada
imgFiltrada.updatePixels();
return imgFiltrada;
}

void TexturedCube() {
  beginShape(QUADS);
  texture(img);

  
  // +Z "front" face
  vertex(-1, -1,  1, 0, 0);
  vertex( 1, -1,  1, 1, 0);
  vertex( 1,  1,  1, 1, 1);
  vertex(-1,  1,  1, 0, 1);
  endShape();
  beginShape(QUADS);
  texture(canny(img));
  // -Z "back" face
  vertex( 1, -1, -1, 0, 0);
  vertex(-1, -1, -1, 1, 0);
  vertex(-1,  1, -1, 1, 1);
  vertex( 1,  1, -1, 0, 1);
  endShape();
  beginShape(QUADS);
  texture(img);
  // +Y "bottom" face
  vertex(-1,  1,  1, 0, 0);
  vertex( 1,  1,  1, 1, 0);
  vertex( 1,  1, -1, 1, 1);
  vertex(-1,  1, -1, 0, 1);
  endShape();
  beginShape(QUADS);
  texture(img);
  // -Y "top" face
  vertex(-1, -1, -1, 0, 0);
  vertex( 1, -1, -1, 1, 0);
  vertex( 1, -1,  1, 1, 1);
  vertex(-1, -1,  1, 0, 1);
  endShape();
  beginShape(QUADS);
  texture(gaussianBlur(img));
  // +X "right" face
  vertex( 1, -1,  1, 0, 0);
  vertex( 1, -1, -1, 1, 0);
  vertex( 1,  1, -1, 1, 1);
  vertex( 1,  1,  1, 0, 1);
  endShape();
  beginShape(QUADS);
  texture(sobel(img));
  // -X "left" face
  vertex(-1, -1, -1, 0, 0);
  vertex(-1, -1,  1, 1, 0);
  vertex(-1,  1,  1, 1, 1);
  vertex(-1,  1, -1, 0, 1);

  endShape();
}

void mouseDragged() {
  float rate = 0.01;
  rotx += (pmouseY-mouseY) * rate;
  roty += (mouseX-pmouseX) * rate;
}
