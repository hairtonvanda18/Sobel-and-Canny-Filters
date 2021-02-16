PImage img;
float rotx = PI/4;
float roty = PI/4;
//Kernel para vertical
float[][] Gx = {{ -1, 0, 1}, 
                { -2,  0, 2}, 
                { -1, 0, 1}};
                    
//Kernel para horizontal                
float[][] Gy = {{ 1, 2, 1}, 
                { 0,  0, 0}, 
                { -1, -2, -1}};
void setup() {
  size(640, 360, P3D);
  img = loadImage("5.jpg");
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
  TexturedCube(Sobel(img));
}
PImage Sobel(PImage img){
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
      double r=Math.min(255,Math.sqrt(Math.pow(sumy,2)+Math.pow(sumx,2)));
      //Inserção do novo pixel
      
      imgFiltrada.pixels[y*img.width + x] = color(Math.round(r));
      //imgFiltrada.pixels[y*img.width + x] = color(sumy);
    }
  }
  
  
// actualizar a imagem filtrada
imgFiltrada.updatePixels();
return imgFiltrada;
}

void TexturedCube(PImage tex) {
  beginShape(QUADS);
  texture(tex);

  // Given one texture and six faces, we can easily set up the uv coordinates
  // such that four of the faces tile "perfectly" along either u or v, but the other
  // two faces cannot be so aligned.  This code tiles "along" u, "around" the X/Z faces
  // and fudges the Y faces - the Y faces are arbitrarily aligned such that a
  // rotation along the X axis will put the "top" of either texture at the "top"
  // of the screen, but is not otherwised aligned with the X/Z faces. (This
  // just affects what type of symmetry is required if you need seamless
  // tiling all the way around the cube)
  
  // +Z "front" face
  vertex(-1, -1,  1, 0, 0);
  vertex( 1, -1,  1, 1, 0);
  vertex( 1,  1,  1, 1, 1);
  vertex(-1,  1,  1, 0, 1);

  // -Z "back" face
  vertex( 1, -1, -1, 0, 0);
  vertex(-1, -1, -1, 1, 0);
  vertex(-1,  1, -1, 1, 1);
  vertex( 1,  1, -1, 0, 1);

  // +Y "bottom" face
  vertex(-1,  1,  1, 0, 0);
  vertex( 1,  1,  1, 1, 0);
  vertex( 1,  1, -1, 1, 1);
  vertex(-1,  1, -1, 0, 1);

  // -Y "top" face
  vertex(-1, -1, -1, 0, 0);
  vertex( 1, -1, -1, 1, 0);
  vertex( 1, -1,  1, 1, 1);
  vertex(-1, -1,  1, 0, 1);

  // +X "right" face
  vertex( 1, -1,  1, 0, 0);
  vertex( 1, -1, -1, 1, 0);
  vertex( 1,  1, -1, 1, 1);
  vertex( 1,  1,  1, 0, 1);

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
