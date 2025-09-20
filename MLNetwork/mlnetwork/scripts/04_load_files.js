


 /* src/04_load_file.js - start*/ 
 function handleClick(elemento) {
    elemento.checked ? label = true : label = false;
  }
  
  function gotData(data) {
    console.log('data: ' + data);
  }
  
  function handleFile(file) {
    console.log(file.data[0]),
    console.log(file.data[1]),
    console.log(file.data[2]),
    objectLoaded = file.data,
    parseobjectLoaded();
  }
  
  function handleColor(data) {
    console.table(data);
  }
  
  function newFile() {
    const obj = [[{ id: 'nodo1.1', x: 0, y: 0, z: 0, nodesize: 10, label: 'nodo1-layer1', zlevel: 0, nodecolor: '#cccccc' }], [{ level: 0, z: 1 }], []];
    objectLoaded = obj, parseobjectLoaded();
  }
  
  
  
  
  function exportFile(){
  let fileName = prompt("Please enter your file name");
  console.log(' nome per il file', fileName)
  
  if (fileName != null) {
  
  
    let objSchema = {
        "scene_pan":{
            "position_x":"0", 
            "position_y":"0", 
            "scale_x":"0.336569927599661", 
            "color":"#000000"
        }, 
        "scene_sphere":{
            "rotation_x":"-0.386047223242714", 
            "rotation_y":"0.378081984494371", 
            "rotation_z":"-0.0257832291643735"
        }, 
        "layers":[], 
        "nodes":[], 
        "edges":[]
    };
  
    objSchema.layers = zzzz;
    objSchema.nodes = sfereJson;
    objSchema.edges = connection;
  
    const content = JSON.stringify(objSchema);
  
    // Crea un oggetto Blob dal contenuto
    const blob = new Blob([content], {type: "application/json"});
  
    // Crea un URL temporaneo per il blob
    const url = URL.createObjectURL(blob);
  
    // Crea un elemento 'a' HTML per il link di download
    const link = document.createElement('a');
    link.href = url;
    link.download = fileName+'.json';
    link.click();
  
    return objSchema;
  
  }else{
    console.log('nessun nome per il file')
  }
  
  }
  
  
  
  function trasforma(mathLabObj) {
  
    let mioOggetto = [];
    let newNode = [];
    let newLayers = [];
    let newConnect = [];
  
    mioOggetto.push(mathLabObj.nodes)
    mioOggetto.push(mathLabObj.layers)
    mioOggetto.push(mathLabObj.edges)
  
    objectLoaded = mioOggetto,
    parseobjectLoaded();
    document.body.setAttribute("file-load", "");
    console.log(mioOggetto);
  }
  
  function loadFile() {
    if(!window.jsonObj){
      var input, file, fr;
  
      if (typeof window.FileReader !== 'function') {
        alert("The file API isn't supported on this browser yet.");
        return;
      }
      input = document.getElementById('fileinput');
      if (!input) {
        alert("Um, couldn't find the fileinput element.");
      }
      else if (!input.files) {
        alert("This browser doesn't seem to support the `files` property of file inputs.");
      }
      else if (!input.files[0]) {
        alert("Please select a file before clicking 'Load'");
      }
      else {
        file = input.files[0];
        fr = new FileReader();
        fr.onload = receivedText;
        fr.readAsText(file);
      }
    }else{
      trasforma(jsonObj)
    }
  
    function receivedText(e) {
      let lines = e.target.result;
      var newArr = JSON.parse(lines);
      console.log(newArr)
      trasforma(newArr)
    }
  
  }
  
  
  document.addEventListener("DOMContentLoaded", function() {
  
    window.jsonObj=(window.jsonObj?window.jsonObj:false);
    if(window.jsonObj){
      trasforma(jsonObj);
    }
  
    const fileInput = document.getElementById('fileinput');
    fileInput.addEventListener('change', function() {
      // callback function
      loadFile();
      // esegui altre operazioni qui
    });
  });
  
  