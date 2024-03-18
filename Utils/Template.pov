#version 3.7;

#include "colors.inc"

global_settings {
    assumed_gamma 1.2
    ambient_light <1, 1, 1> // <red, green, blue}
}

#default {
finish {
    ambient 0.1
    diffuse 0.9
  }
}

camera { 
        location  <0,15,0>
        look_at 0
        up <0,1,0>
		right x
        angle {{ Angle }}
}

light_source{
    <0,20,5> White
}

#macro Vector (x_offset, y_offset, phi, theta, r, g, b)
    cone {
        <0, -1, 0>, 1.0, <0, 1, 0>, 0.0
        texture {

            pigment {
                color rgb<r, g, b>
            }
            
            finish {
                specular 0.1
            }
        }
    scale <0.5, 0.5, 0.5>
    rotate <0, 0, -degrees(theta)>
    rotate <0, -degrees(phi), 0>
    translate <x_offset, 0, y_offset>
    }
#end

