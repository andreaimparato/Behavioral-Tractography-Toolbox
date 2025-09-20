
 /* src/00_frontend_function.js - start*/ 
 function layerList() {
    
  let template = '';
  const layerListEl = document.querySelector('.drawer_layers .dynamic-content');
  for (const [index, obj] of zzzz.entries()) {
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
  const buttonEl = document.createElement('button');
  buttonEl.setAttribute('xx-data-id', index);
  buttonEl.setAttribute('xx-data-level', obj_level);

  buttonEl.textContent = obj_level;
  //template += buttonEl.outerHTML;

  template += `<div class="details-item layers-details-item">
                  <div class="title">${obj_level}</div>
                  <div class="info">
                      <div class="manage">
                          <button class="manage-show" xx-data-id="${index}" xx-data-level="${obj_level}">show/hide</button>
                      </div>
                  </div>
              </div>`;
  }


  var hasDataBuilt = $('.drawer.drawer_layers .wrapper .dynamic-content').is('[data-built]');

  if (!hasDataBuilt && template) {
      document.querySelector('.drawer.drawer_layers .wrapper .dynamic-content').innerHTML = template;
      $('.drawer.drawer_layers .wrapper .dynamic-content').attr('data-built', true)
  }

  //layerListEl.innerHTML = template;
  layerListEl.setAttribute('fill', 'true');
}






function openLabelPanel(index) {
  $('.textPopup').prependTo('.node-details-item:eq(' + index + ') .panel');


  let nodeColor = sfereJson[index].color.replace('#', '');
  console.log('>>>>>>>>>>>>>>> colore: '+nodeColor)
  $('.position_item_controls input').val('#'+nodeColor)

  var isBold = sfereJson[index].isBold ? 'bold' : 'regular';

  var displayName = sfereJson[index].displayName ? sfereJson[index].displayName : sfereJson[index].name;
  $('.position_item .item_name input').val(displayName);

  $('input#labelAdjustX').val(sfereJson[index].labelAdjustX);
  $('input#labelAdjustY').val(sfereJson[index].labelAdjustY);


  const label = sfereJson[index].displayName?sfereJson[index].displayName:sfereJson[index].name;
  const labelPosition = sfereJson[index].labelPosition || false;
  const textPopupEl = document.querySelector('.textPopup');
  textPopupEl.setAttribute('xxx-label-style', isBold);
  textPopupEl.setAttribute('xxx-node-color', '#' + nodeColor);
  textPopupEl.setAttribute('xxx-label-position', labelPosition);
  textPopupEl.setAttribute('xxx-node-index', index);
  textPopupEl.setAttribute('xxx-display-name', displayName);
  textPopupEl.setAttribute('xxx-axis-x', sfereJson[index].labelAdjustX);
  textPopupEl.setAttribute('xxx-axis-y', sfereJson[index].labelAdjustY);
  textPopupEl.style.display = 'block';
  $('.control_button').attr('xxx-node-index', index);
  //const labelAdjustX = document.querySelector('#labelAdjustX');
  //const labelAdjustY = document.querySelector('#labelAdjustY');

  //const labelEl = textPopupEl.querySelector('span.label');
  //labelEl.textContent = label;
  //labelEl.style.backgroundColor = '#' + sfereJson[index].color.replace('#', '');
  //textField = document.getElementById("nome_nodo");
  //textField.value =label;

  console.log(sfereJson[index]);

}

function labelSetSize(mode){
      switch (mode) {
          case '+':
              labelSizeCoef=labelSizeCoef<1000?labelSizeCoef+=1:false
              console.log('labelSizeCoef: ',labelSizeCoef)
            break;
          case '-':
              labelSizeCoef=labelSizeCoef>0?labelSizeCoef-=1:false
              console.log('labelSizeCoef: ',labelSizeCoef)
          default:
      }
}
function labelSetSizeSingular(mode, id){
  let coefficient;
  switch (mode) {
      case '+':
          coefficient = sfereJson[id].labelSize?sfereJson[id].labelSize:0;
          sfereJson[id].labelSize=coefficient+1
        break;
      case '-':
          coefficient = sfereJson[id].labelSize?sfereJson[id].labelSize:0;
          sfereJson[id].labelSize=(coefficient>0)?coefficient-1:0;
      default:
  }
}


function nodeColorChange(element) {
  var newColorValue = element.value;
  $('.textPopup').attr('xxx-node-color', newColorValue)
  //console.log('Ultimo valore selezionato:', newColorValue);
}

function savec () {
  let fileName = prompt("Please enter your image name");
  saveCanvas(c, fileName, "png");
};


function closeLabelPopup() {
  document.querySelector('.textPopup').style.display = 'none';
}

function saveLabelPopup(){
  //console.log()
  let labelPosition = $('.textPopup').attr('xxx-label-position');
  let textStyle = $('.textPopup').attr('xxx-label-style');
  let index =  parseInt($('.textPopup').attr('xxx-node-index'));
  let nodeColor =  $('.textPopup').attr('xxx-node-color');
  let displayName =  $('.textPopup').attr('xxx-display-name');
  //let displayName = $('input#nome_nodo').val();

  //const labelAdjustX = document.querySelector('#labelAdjustX');
  //const labelAdjustY = document.querySelector('#labelAdjustY');

  const labelAdjustX = $('.textPopup').attr('xxx-axis-x');
  const labelAdjustY = $('.textPopup').attr('xxx-axis-y');

  /*sfereJson[index].labelAdjustX = labelAdjustX;
  sfereJson[index].labelAdjustY = labelAdjustY;*/
  sfereJson[index].labelAdjustX = labelAdjustX;
  sfereJson[index].labelAdjustY = labelAdjustY;
  sfereJson[index].labelPosition = labelPosition;
  sfereJson[index].displayName = displayName;

  sfereJson[index].isBold = textStyle == "bold" ? true : false;
  //sfereJson[index].name = name;

  console.log('>>>>>>>>>>>>> colore 2: ', nodeColor)
  sfereJson[index].color = nodeColor;


  document.querySelector('.textPopup').style.display = 'none';

  $('.position_contaniner').attr('data-changed', false);
}


function showHideSphere(index) {
try{
  //console.log('showHideSphere ' + index);
  let hidden = sfereJson[index].hide ? false : true;
  sfereJson[index].hide = hidden;
  $('.manage button[xx-data-id="' + index + '"]').attr('cc-item-hide', hidden);
  let nameItem = $('.manage button[xx-data-id="' + index + '"]').attr('xx-data-node') + '+' + $('.manage button[xx-data-id="' + index + '"]').attr('xx-data-layer')
  //console.log('nameItem:', nameItem)
  nodeListPosition[nameItem].hide = hidden;
}
catch(e){
  console.log(e)
}

}

function keyPressedShowAllNode(){
  let buttons = document.querySelectorAll('.listanodi button[cc-item-hide]');
  for (let i = 0; i < buttons.length; i++) {
    buttons[i].setAttribute('cc-item-hide', 'false');
  }
  Object.keys(nodeListPosition).map(function (objectKey, index) {
    showHideSphere(index)
  });
}







function showHideLayer(index) {
  /*
  let hidden = zzzz[index].hide ? false : true;
  zzzz[index].hide = hidden;
  const layerButtons = document.querySelectorAll('.listalayer button');
  layerButtons.forEach((button) => {
      if (button.getAttribute('xx-data-id') == index) {
          button.setAttribute('cc-item-hide', hidden);
      }
  });
  */


let hidden = zzzz[index].hide ? false : true;
zzzz[index].hide = hidden;
$('.drawer_layers button[xx-data-id="' + index + '"]').attr('cc-item-hide', hidden);
}

function isNumeric(str) {
  // Verifica se la stringa è vuota o contiene solo spazi
  if (str.trim() === "") {
    return false;
  }
  
  // Controlla se la stringa contiene solo numeri usando una regex
  return /^\d+$/.test(str);
}


function nodeList() {
  //console.log('*** nodeList');
  let template = '';
  let template_ = '';
  const nodeButtons = document.querySelectorAll('.listanodi button[cc-item-hide]');

  //template += '<div class="drawer-title">Selezionare nodi</div>';

  for (let i = 0; i < sfereJson.length; i++) {
    const obj = sfereJson[i];
    if(! isNumeric(obj?.name)){
      let obj_id = obj.name;
      let obj_label = obj.name;
      let obj_zlayer = parseFloat((obj.layer).replace('Layer', '')) - 1;
      let obj_zlevel = obj.layer;
      let obj_color = obj.color;
      template_ += `
              <button xx-data-id="${i}" xx-data-zlevel="${obj_zlevel}" xx-data-layer="${obj_zlevel}" xx-data-id="${obj_id}" xx-data-node="${obj_id}">
                  <span style="background-color:#${obj_color.replace('#', '')}" class="node_color"></span>
                  <span class="node_label">${obj_label}</span>
                  <span class="node_layer">L${obj_zlayer}</span>
                  <span class="node_sh" onclick="showHideSphere(${i})">O</span>
                  <span class="node_t" onclick="openLabelPanel(${i})">T</span>
              </button>`;
      template += `<div class="details-item node-details-item">
                      <div class="title">${obj_label}</div>
                      <div class="info">
                          <div class="layer">Layer ${obj_zlayer}</div>
                          <div class="manage">
                              <button class="manage-show" xx-data-id="${i}" xx-data-zlevel="${obj_zlevel}" xx-data-layer="${obj_zlevel}" xx-data-id="${obj_id}" xx-data-node="${obj_id}" onclick="showHideSphere(${i})">show/hide</button>
                              <button class="manage-config" onclick="openLabelPanel(${i})">config</button>
                          </div>
                      </div>
                      <div class="panel"></div>
                  </div>`;
 
                  
      if (nodeButtons[i]) {
          nodeButtons[i].setAttribute('cc-item-hide', 'false');
      }
    }
  }

  

  //document.querySelector('.menu-layer.listanodi .nodi_wrapper').innerHTML = template_;

  var hasDataBuilt = $('.drawer.drawer_nodi .wrapper .dynamic-content').is('[data-built]');

  if (!hasDataBuilt && template) {
      document.querySelector('.drawer.drawer_nodi .wrapper .dynamic-content').innerHTML = template;
      $('.drawer.drawer_nodi .wrapper .dynamic-content').attr('data-built', true)
  }

  document.querySelector('.menu-layer.listanodi').setAttribute('fill', 'true');
}

