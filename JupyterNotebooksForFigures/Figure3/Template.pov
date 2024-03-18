#include "colors.inc"


global_settings {
    assumed_gamma 1.2
    ambient_light <1, 1, 1>
}

#default {
finish {
    ambient 0.1
    diffuse 0.9
  }
}

#declare L = {{ L }};

camera { 
    location  <L, L, -L>
    look_at <0, 0.55*L, 0>
    right <0.5,0,0>
}

light_source{
    <L, L, -L> White
}


#declare sphere_radius = 0.1;

{{ Sphere Sweep }}