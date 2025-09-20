
 /* src/03_edge_draw.js - start*/ 

 function linee(connection) {


    for (let obj of connection) {
  
      let obj_from = obj.src;
      let obj_to = obj.trg;
  
      let layerFrom = nodeListPosition[obj_from].layer;
      let layerTo = nodeListPosition[obj_to].layer;
      let isHide_to = nodeListPosition[obj_to].hide?nodeListPosition[obj_to].hide:false;
      let isHide_from = nodeListPosition[obj_from].hide?nodeListPosition[obj_from].hide:false;
  
      if (!isHide_to && !isHide_from && !zzzz[layerFrom].hide && !zzzz[layerTo].hide) {
        if (
          nodeListPosition[obj_to] &&
          nodeListPosition[obj_from] &&
          nodeListPosition[obj_from].x &&
          nodeListPosition[obj_to].x
        ) {
  
          let xFrom = nodeListPosition[obj_from].x,
            yFrom = nodeListPosition[obj_from].y,
            zFrom = nodeListPosition[obj_from].z * layerZcoef,
  
            xTo = nodeListPosition[obj_to].x,
            yTo = nodeListPosition[obj_to].y,
            zTo = nodeListPosition[obj_to].z * layerZcoef;
  
          let colore = obj.color,
            stroke_weight = obj.size * lineeCoef;
  
          push();
          smooth();
          stroke(colore);
          smooth();
          strokeWeight(stroke_weight);
          smooth();
          specularMaterial(90);
          line(xFrom, yFrom, zFrom, xTo, yTo, zTo);
          let sphereToSize = sfereJson.filter(({ name }) => name === obj_to.split('+')[0]).filter(({ layer }) => layer === obj_to.split('+')[1])[0].scale_x;
          let sphereFromSize = sfereJson.filter(({ name }) => name === obj_from.split('+')[0]).filter(({ layer }) => layer === obj_from.split('+')[1])[0].scale_x;
          let sphereA = { position: createVector(xFrom, yFrom, zFrom), radius: ((12 + 5 * sphereFromSize) * sfereCoef) };
          let sphereB = { position: createVector(xTo, yTo, zTo), radius: ((12 + 5 * sphereToSize) * sfereCoef) };
          push()
          //console.log("freccia: ", obj.arrow)
          if (obj.arrow === "true") {
            if (xFrom > xTo) {
              let buildArrow = orientCylinder(sphereA, sphereB, stroke_weight, colore, stroke_weight, obj_from, obj_to)
            } else {
              let buildArrow = orientCylinder_invert(sphereB, sphereA, stroke_weight, colore, stroke_weight, obj_from, obj_to)
              //let buildArrow = orientCylinder(sphereA, sphereB, stroke_weight, colore, stroke_weight, obj.from, obj.to)
            }
          }
          pop()
          pop();
        }
      }
    }
    return true;
  }
  
  
  
  
  
  
  function orientCylinder(sphere1, sphere2, size, colore, lineSize, froml, tol) {
    push();
    let distance = sphere1.position.dist(sphere2.position);
    let direction = p5.Vector.sub(sphere2.position, sphere1.position).normalize();
    let middlePoint1 = p5.Vector.add(sphere1.position, p5.Vector.mult(direction, sphere1.radius));
    let middlePoint2 = p5.Vector.sub(sphere2.position, p5.Vector.mult(direction, sphere2.radius));
    let middle = p5.Vector.add(middlePoint1, middlePoint2).div(2);
    let axis = createVector(0, 1, 0).cross(direction);
    let angle = createVector(0, 1, 0).angleBetween(direction);
    translate(middle);
    rotate(angle, axis);
    noStroke();
    fill(colore)
    //cylinder(size, distance - sphere1.radius - sphere2.radius - 20);
    push();
    translate(0, (distance - sphere1.radius - sphere2.radius) / 2, 0);
    specularMaterial(60);
    cone(3 + lineSize, 16 + lineSize);
    fill(0, 0, 0)
    let fontSizeLarge = (12 + textSizeCoeff)
    textSize(fontSizeLarge);
    let testoLabel = froml + ' > ' + tol;
    // text(testoLabel, -30, -30);
    pop();
    pop();
    return true;
  }
  
  
  
  
  
  function orientCylinder_invert(sphere1, sphere2, size, colore, lineSize, froml, tol) {
    push();
    let distance = sphere1.position.dist(sphere2.position);
    let direction = p5.Vector.sub(sphere2.position, sphere1.position).normalize();
    let middlePoint1 = p5.Vector.add(sphere1.position, p5.Vector.mult(direction, sphere1.radius));
    let middlePoint2 = p5.Vector.sub(sphere2.position, p5.Vector.mult(direction, sphere2.radius));
    let middle = p5.Vector.add(middlePoint1, middlePoint2).div(2);
    let axis = createVector(0, 1, 0).cross(direction);
    let angle = createVector(0, 1, 0).angleBetween(direction);
    translate(middle);
    rotate(angle, axis);
    noStroke();
    fill(colore)
    //cylinder(size, distance - sphere1.radius - sphere2.radius - 20);
    push();
    translate(0, -(distance - sphere1.radius - sphere2.radius) / 2, 0);
    rotateX(PI);
    specularMaterial(60);
    cone(3 + lineSize, 16 + lineSize);
  
    fill(0, 0, 0)
    let fontSizeLarge = (12 + textSizeCoeff)
    textSize(fontSizeLarge);
    let testoLabel = froml + ' > ' + tol;
    // text(testoLabel, -30, -30);
  
    pop();
    pop();
    return true;
  }