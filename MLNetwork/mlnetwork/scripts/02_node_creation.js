
 /* src/02_node_creation.js - start*/ 

 function sfere(sfereData) {

  let index = 0;
  for (let obj of sfereData) {

    let obj_id = obj.name;
    let obj_displayName = obj.displayName?obj.displayName:false;
    let obj_x = parseFloat(obj.position_z);
    let obj_y = 0 - parseFloat(obj.position_y);
    let obj_nodesize = parseFloat(obj.scale_x);
    let obj_label = obj.name;
    let obj_zlevel = parseFloat((obj.layer).replace('Layer', '')) - 1;
    let obj_hide = obj.hide?obj.hide:false;
    let obj_isBold = obj.isBold;
    let obj_labelSizeCoef = obj.labelSize;

    let obj_labelAdjustX = obj.labelAdjustX?parseInt(obj.labelAdjustX):0;
    let obj_labelAdjustY = obj.labelAdjustY?parseInt(obj.labelAdjustY):0;

    let obj_nodecolor = obj.color;

    font = (obj_isBold=="true" || obj_isBold==true)?font=font_bold:font=font_light;

    if (obj && obj_zlevel in zzzz && !obj_hide && !zzzz[obj_zlevel].hide) {

      let obj_z = zzzz[obj_zlevel].position_x;

      let layername = `${obj_id}+Layer${obj_zlevel + 1}`;
      nodeListPosition[layername] = { x: obj_x, y: obj_y, z: obj_z, layer: obj_zlevel };

      push();
      translate(obj_x, obj_y, obj_z * layerZcoef);
      noStroke();
      fill(obj_nodecolor);
      let sphereSizeR = 6 * obj_nodesize * sfereCoef;

      specularMaterial(40);
      sphere(sphereSizeR);

      if (label) {
        labelsize = obj_labelSizeCoef?parseInt(obj_labelSizeCoef):0;
        let labelSizeCoef_ = labelSizeCoef+labelsize;

        textFont(font, 11 + labelSizeCoef_);
        rotateY(textYrotation);
        rotateX(textXrotation);
        translate(0, 0, 10);

        let fontSizeLarge = 11 + labelSizeCoef_;
        textSize(fontSizeLarge);
        
        

        fill(obj_nodecolor);
        let testoLabel = obj_displayName?obj_displayName:obj_label;
        let bounding_box2 = font.textBounds(testoLabel, 20, 60, fontSizeLarge+(labelSizeCoef_/2));

        noStroke();
        let verticalbox = -(sphereSizeR + 20);
        let vertical = -(sphereSizeR + 10);

        if (obj.labelPosition === 'top') {
          vertical = -(sphereSizeR + 10);
          verticalbox = -(sphereSizeR + 20);
        }

        if (obj.labelPosition === 'bottom') {
          vertical = sphereSizeR + 20;
          verticalbox = sphereSizeR + 10;
        }

        rect(obj_labelAdjustX-2-(labelSizeCoef_), verticalbox+obj_labelAdjustY-(labelSizeCoef_), bounding_box2.w + 4, bounding_box2.h + 4);
        fill(0);
        translate(0, 0, 1);
        

        text(testoLabel, obj_labelAdjustX, vertical+obj_labelAdjustY);
      }

      translate(0, 0, 1);
      pop();
    }
    index++;
  }

  let makeConnection = linee(connection);
  if (document.querySelectorAll('.menu-layer.listanodi button').length === 0) nodeList();
}
