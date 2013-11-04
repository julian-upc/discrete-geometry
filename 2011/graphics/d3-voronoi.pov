#version 3.5;

#include "colors.inc"
#include "textures.inc"

global_settings {
  assumed_gamma 1.0
  max_trace_level 50
}

// ----------------------------------------


sky_sphere {
  pigment {
    gradient y
    color_map {
      [ (1-cos(radians(-30)))/2 color CornflowerBlue ]
      [ (1-cos(radians(160)))/2 color MidnightBlue ]
    }
    scale 2
    translate -1
  }
}



// first, the camera position
camera {
  //  orthographic
  location <-2.5,-2.5,5>
  sky <0,0,1>
  look_at <2,1,0.5>
}

// now, some light
light_source {
  <-20,-20,20>
  color rgb <1,1,1>
}

light_source {
  <0,0,20>
  color rgb <1,1,1>
}

// the spheres

#macro colored_sphere(RA, R, G, B, T)
  sphere {
    <0,0,0>, RA
    pigment { color rgbt<R,G,B,T> }
  }
#end
    
#declare blue_sphere   = colored_sphere(0.70710678, 0,0,1,0.7)
#declare lightblue_sphere   = colored_sphere(0.70710678, 0,0,1,0.9)

#macro fcc (the_sphere)
  union {
    object { the_sphere translate<0,0,0> } 
    object { the_sphere translate<0,2,0> } 
    object { the_sphere translate<1,-1,0> } 
    object { the_sphere translate<1,1,0> } 
    object { the_sphere translate<1,3,0> } 
    object { the_sphere translate<2,0,0> } 
    object { the_sphere translate<2,2,0> } 
    object { the_sphere translate<3,-1,0> } 
    object { the_sphere translate<3,1,0> } 
    object { the_sphere translate<3,3,0> } 
    object { the_sphere translate<4,0,0> } 
    object { the_sphere translate<4,2,0> } 
    
    object { the_sphere translate<0,1,1> } 
    object { the_sphere translate<1,0,1> } 
    object { the_sphere translate<1,2,1> } 
    object { the_sphere translate<2,-1,1> } 
    object { the_sphere translate<2,1,1> } 
    object { the_sphere translate<2,3,1> } 
    object { the_sphere translate<3,0,1> } 
    object { the_sphere translate<3,2,1> } 
    object { the_sphere translate<4,1,1> } 
    
    object { the_sphere translate<0,0,2> } 
    object { the_sphere translate<0,2,2> } 
    object { the_sphere translate<1,-1,2> } 
    object { the_sphere translate<1,1,2> } 
    object { the_sphere translate<1,3,2> } 
    object { the_sphere translate<2,0,2> } 
    object { the_sphere translate<2,2,2> } 
    object { the_sphere translate<3,-1,2> } 
    object { the_sphere translate<3,1,2> } 
    object { the_sphere translate<3,3,2> } 
    object { the_sphere translate<4,0,2> } 
    object { the_sphere translate<4,2,2> } 
  }
#end

#macro ptA() <0,0,1> #end
#macro ptB1() < 0.5, 0.5, 0.5> #end
#macro ptB2() <-0.5, 0.5, 0.5> #end
#macro ptB3() <-0.5,-0.5, 0.5> #end
#macro ptB4() < 0.5,-0.5, 0.5> #end
#macro ptC1() <1,0,0> #end
#macro ptC2() <0,1,0> #end
#macro ptC3() <-1,0,0> #end
#macro ptC4() <0,-1,0> #end
#macro ptD1() < 0.5, 0.5,-0.5> #end
#macro ptD2() <-0.5, 0.5,-0.5> #end
#macro ptD3() <-0.5,-0.5,-0.5> #end
#macro ptD4() < 0.5,-0.5,-0.5> #end
#macro ptE() <0,0,-1> #end

#macro voronoi(R,G,B,T)
  union {
    polygon {4, ptA(), ptB1(), ptC2(), ptB2() pigment{color rgbt<R,G,B,T>}}
    polygon {4, ptA(), ptB2(), ptC3(), ptB3() pigment{color rgbt<R,G,B,T>}}
    polygon {4, ptA(), ptB3(), ptC4(), ptB4() pigment{color rgbt<R,G,B,T>}}
    polygon {4, ptA(), ptB4(), ptC1(), ptB1() pigment{color rgbt<R,G,B,T>}}

    polygon {4, ptC1(), ptD1(), ptC2(), ptB1() pigment{color rgbt<R,G,B,T>}}
    polygon {4, ptC2(), ptD2(), ptC3(), ptB2() pigment{color rgbt<R,G,B,T>}}
    polygon {4, ptC3(), ptD3(), ptC4(), ptB3() pigment{color rgbt<R,G,B,T>}}
    polygon {4, ptC4(), ptD4(), ptC1(), ptB4() pigment{color rgbt<R,G,B,T>}}

    polygon {4, ptE(), ptD1(), ptC2(), ptD2() pigment{color rgbt<R,G,B,T>}}
    polygon {4, ptE(), ptD2(), ptC3(), ptD3() pigment{color rgbt<R,G,B,T>}}
    polygon {4, ptE(), ptD3(), ptC4(), ptD4() pigment{color rgbt<R,G,B,T>}}
    polygon {4, ptE(), ptD4(), ptC1(), ptD1() pigment{color rgbt<R,G,B,T>}}
}
#end
  
  
fcc(lightblue_sphere)
object {
  voronoi(1,0,0,0.5) translate<1,2,1>
}