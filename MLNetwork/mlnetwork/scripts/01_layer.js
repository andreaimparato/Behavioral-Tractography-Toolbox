

 /* src/01_layer.js - start*/ 

 function layer() {
    for (obj of zzzz) {
  
      let obj_level = obj.name;
      let obj_z = parseFloat(obj.position_x);
      let obj_floor_current_color = obj.floor_current_color;
      let obj_geometry_parameters_width = parseFloat(obj.geometry_parameters_width);
      let obj_last_layer_scale = "1";
      let obj_name = obj.name;
      let obj_position_x = parseFloat(obj.position_x);
      let obj_position_y = parseFloat(obj.position_y);
      let obj_position_z = parseFloat(obj.position_z);
      let obj_rotation_x = parseFloat(obj.rotation_x);
      let obj_rotation_y = parseFloat(obj.rotation_y);
      let obj_rotation_z = parseFloat(obj.rotation_z);
  
      if (!obj.hide) {
        push(),
          squareColor = color(255, 255, 255, .50),
          noStroke(),
          fill('rgba(255,255,255, 0.5)'),
          translate(0, 0, (obj_z * layerZcoef) - 5);
  
  /*
          for (var x = 0; x < width; x += width / 10) {
            for (var y = 0; y < height; y += height / 10) {
              stroke(0);
              strokeWeight(1);
              line(x, 0, x, height);
              line(0, y, width, y);
            }
          };
  */
  
          specularMaterial(0),
          plane(1e3, 1e3),
          translate(-450, -420, 10),
          squareColor = color(255, 255, 255, .85),
          fill('rgba(255,255,255, 0.85)'),
          textSize(60);
        if (obj_level > -1) text('level' + obj_level, -30, -10);
        pop();
      }
    }
    if (document.querySelectorAll('.menu-layer.listalayer button').length == 0) layerList();
  }
  
  
  
  function parseobjectLoaded() {
      sfereJson = objectLoaded[0],
      zzzz = objectLoaded[1],
      connection = objectLoaded[2];
      console.log('connection2 : ', connection)
      //oggettoSfere = objectLoaded[0];
  }
  
  
  function base() {
    push(), fill(80), rotateY(0), box(500, 3, 500), pop();
  }
  
  
  
  function textRotationY(quanto) {
    textYrotation += quanto;
  }
  function textRotationX(quanto) {
    textXrotation += quanto;
  }
  
  