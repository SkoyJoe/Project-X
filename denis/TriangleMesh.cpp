// ========================================================================= //
// Authors: Roman Getto, Matthias Bein                                       //
// mailto:roman.getto@gris.informatik.tu-darmstadt.de                        //
//                                                                           //
// GRIS - Graphisch Interaktive Systeme                                      //
// Technische Universität Darmstadt                                          //
// Fraunhoferstrasse 5                                                       //
// D-64283 Darmstadt, Germany                                                //
//                                                                           //
// Content: Simple class for reading and rendering triangle meshes            //
// ========================================================================== //

#include <iostream>
#include <fstream>
#include <float.h>
#include <GL/glut.h>
#include "TriangleMesh.h"

void TriangleMesh::calculateNormals() {
  normals.resize(vertices.size());
  // TODO: calculate normals for each vertex
  std::vector<Normal> triangleNorms;
  triangleNorms.resize(triangles.size());
  //Norms of triangles ( (x-y) x (x-z) )
  //
 std::cout<<"================="<<std::endl;
 std::cout<<"LOADING NORMALS"<<std::endl;
 std::cout<<"================="<<std::endl;
/*  for (int i=0;i<triangles.size();i++){
  Vec3f a = vertices[triangles[i].x] - vertices[triangles[i].y];
  Vec3f b = vertices[triangles[i].x] - vertices[triangles[i].z];
  
  // a x b
  Normal temp;
  temp.x = a.y*b.z - a.z*b.y;
  temp.y = a.z*b.x - a.x*b.z;
  temp.z = a.x*b.y - a.y*b.x;
  
  triangleNorms[i]=temp;
}
  //finds 3 triangles that has vertex[i] and sums it

  for(int i=0; i<vertices.size();i++){
    Normal temp(0,0,0);
    int flag=0;
     for(int j=0; j<triangles.size();j++){
        if (flag>=3) break;
        if(triangles[j].x ==i || triangles[j].y ==i || triangles[j].z ==i){
        temp+=triangleNorms[j];
        flag++;
       
}
        
    

}
     
    if(flag!=3)std::cout<<"Something is wrong with normals or the object is not completely in 3D like cloth "<<flag<<" Point"<<i<<std::endl;
    normals[i]=temp;
    
}
*/
for (int i=0;i<triangles.size();i++){
  Vec3f a = vertices[triangles[i].x] - vertices[triangles[i].y];
  Vec3f b = vertices[triangles[i].x] - vertices[triangles[i].z];
  
  // a x b
  Normal temp;
  temp.x = a.y*b.z - a.z*b.y;
  temp.y = a.z*b.x - a.x*b.z;
  temp.z = a.x*b.y - a.y*b.x;
  //max 3 normals for vertex 
  normals[triangles[i].x]+=temp;
  normals[triangles[i].y]+=temp;
  normals[triangles[i].z]+=temp;
}
  // normalize normals
  for (Normals::iterator nit=normals.begin(); nit!=normals.end(); ++nit) {
     //the normalize() function returns a boolean which can be used if you want to check for erroneous normals
	 (*nit).normalize();
   }
 std::cout<<"================="<<std::endl;
 std::cout<<"  LOADING ENDED  "<<std::endl;
 std::cout<<"================="<<std::endl;
}

// ===============================
// === CONSTRUCTOR, DESTRUCTOR ===
// ===============================

TriangleMesh::TriangleMesh() {
  clear();
}

TriangleMesh::~TriangleMesh() {
  clear();
}

void TriangleMesh::clear() {
  // clear mesh data
  vertices.clear();
  triangles.clear();
  normals.clear();
}

// ================
// === RAW DATA ===
// ================

vector<Vec3f>& TriangleMesh::getPoints() {
  return vertices;
}
vector<Vec3i>& TriangleMesh::getTriangles() {
	return triangles;
}

vector<Vec3f>& TriangleMesh::getNormals() {
  return normals;
}

void TriangleMesh::flipNormals() {
  for (Normals::iterator it = normals.begin(); it != normals.end(); ++it) {
    (*it) *= -1.0;
  }
}

// =================
// === LOAD MESH ===
// =================

void TriangleMesh::loadLSA(const char* filename) {  
  std::ifstream in(filename);
  if (!in.is_open()) {
    cout << "loadLSA: can not open " << filename << endl;
    return;
  }
  char s[256];
  in >> s;
  // first word: LSA
  if (!(s[0] == 'L' && s[1] == 'S' && s[2] == 'A')) return;
  // get number of vertices nv, faces nf, edges ne and baseline distance
  std::cout<<"================"<<std::endl;
  std::cout<<"       LSA      "<<std::endl;
  std::cout<<"================"<<std::endl;
  int nv, nf, ne;
  float baseline;
  in >> nv;
  in >> nf;
  in >> ne;
  in >> baseline;
  if (nv <= 0 || nf <= 0) return;
  // clear any existing mesh
  clear();
  // read vertices  
  vertices.resize(nv);
  // TODO: read alpha, beta, gamma for each vertex and calculate verticex coordinates
  std::cout<<"================"<<std::endl;
  std::cout<<"LOADING VERTICES"<<std::endl;
  std::cout<<"================"<<std::endl;
  for(int i=0;i<nv;i++){
   float x,y,z,alpha,beta,gamma;
   	in>>alpha;
   	in>>beta;
   	in>>gamma; 
   	z=-baseline/(tan(alpha*M_PI/180)+tan(beta*M_PI/180));
   	x= z*tan(beta*M_PI/180) ;
   	//y Berechnung
   	y=sqrt(x*x+z*z)*tan(gamma*M_PI/180);
   	
   	vertices[i]=Vertex(x,y,z); 
 std::cout<< vertices[i]<<std::endl;
} 

  // read triangles
  triangles.resize(nf);
  // TODO: read all triangles from the file
std::cout<<"================="<<std::endl;
 std::cout<<"LOADING TRIANGLES"<<std::endl;
 std::cout<<"================="<<std::endl;
  triangles.resize(nf);
  // TODO: read triangles from the file
  for(int j=0;j<nf;j++){
  int dump, first ,second, third;
  in>>dump;
  in>>first;
  in>>second;
  in>>third;
 
 triangles[j]=Triangle(first,second,third);
 std::cout<< triangles[j]<<std::endl;



}
 std::cout<<"================="<<std::endl;
 std::cout<<" LOADING ENDED   "<<std::endl;
 std::cout<<"================="<<std::endl;
  // calculate normals
	calculateNormals();
}

void TriangleMesh::loadOFF(const char* filename) {
  std::ifstream in(filename);
  if (!in.is_open()) {
    cout << "loadOFF: can not find " << filename << endl;
    return;
  }
  char s[256];
  in >> s;
  // first word: OFF
  if (!(s[0] == 'O' && s[1] == 'F' && s[2] == 'F')) return;
  // get number of vertices nv, faces nf, edges ne and (baseline distance)
    std::cout<<"================"<<std::endl;
  std::cout<<"       OFF     "<<std::endl;
  std::cout<<"================"<<std::endl;
  int nv, nf, ne;
  in >> nv;
  in >> nf;
  in >> ne;
  if (nv <= 0 || nf <= 0) return;
  // clear any existing mesh
  clear();
  // read vertices 
  std::cout<<"================"<<std::endl;
  std::cout<<"LOADING VERTICES"<<std::endl;
  std::cout<<"================"<<std::endl;
  vertices.resize(nv);
  // TODO: read all vertices from the file
  for(int i=0;i<nv;i++){
	float x,y,z;
        in>>x;
        in>>y;
        in>>z;
 
  vertices[i]=Vertex(x,y,z);
 std::cout<< vertices[i]<<std::endl;
	}
  // read triangles
 std::cout<<"================="<<std::endl;
 std::cout<<"LOADING TRIANGLES"<<std::endl;
 std::cout<<"================="<<std::endl;
  triangles.resize(nf);
  // TODO: read triangles from the file
  for(int j=0;j<nf;j++){
  int dump, first ,second, third;
  in>>dump;
  in>>first;
  in>>second;
  in>>third;

 triangles[j]=Triangle(first,second,third);
 std::cout<< triangles[j]<<std::endl;



}
 std::cout<<"================="<<std::endl;
 std::cout<<" LOADING ENDED   "<<std::endl;
 std::cout<<"================="<<std::endl;

  // calculate normals
	calculateNormals();
}

// ==============
// === RENDER ===
// ==============

void TriangleMesh::draw() {
  if (triangles.size() == 0) return;
  // TODO: draw triangles with immediate mode
glBegin(GL_TRIANGLES);
  for (int i = 0; i < triangles.size(); ++i)
  {
    glVertex3d(vertices[triangles[i].x].x, vertices[triangles[i].x].y, vertices[triangles[i].x].z);
    glVertex3d(vertices[triangles[i].y].x, vertices[triangles[i].y].y, vertices[triangles[i].y].z);
    glVertex3d(vertices[triangles[i].z].x, vertices[triangles[i].z].y, vertices[triangles[i].z].z);
  }
  glEnd();
}