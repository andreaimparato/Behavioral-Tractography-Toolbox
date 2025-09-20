
 /* src/06_p5_scheme.js - start*/ 
 var fr = 28;
 var lineeCoef = 2;
 var sfereCoef = 1;
 var labelSizeCoef=0;
 var layerZcoef = 1;
 var textYrotation = 0;
 var textXrotation = 0;
 var textSizeCoeff = 0;
 var layerMenuBuild;
 var offset = 16;
 var oggettoSfere;
 var texts = [];
 var connectionCount = 0;
 var sfereLabelConf = {};
 var setFirstViewStatus = true;
 
 
 var sfereJson = [],
   connection = [],
   font,
   font_bold,
   font_light,
   inputx,
   img,
   zzzz = [],
   bgcolorval = '#ffffff',
   nodeListPosition = {},
   label = !![],
   objectLoaded,
   c,
   xxx,
   yyy,
   zzz;
 
 
 
 function setup() {
   frameRate(fr);
   pixelDensity(4);
   c = createCanvas(windowWidth, windowHeight, WEBGL);
   setAttributes('antialias', !![]);
   easycam = new Dw['EasyCam'](this['_renderer'], { 
     distance: 5000 
   });
   textFont(font_light);
   background(0, 0, 0);
   registra.start();
   //noLoop();
 }
 
 
 
 function preload() {
   /*
   font = loadFont(
     "https://cdnjs.cloudflare.com/ajax/libs/topcoat/0.8.0/font/SourceCodePro-Bold.otf"
   );
   */
 
   font_light = loadFont(
     "/font/Roboto-Regular.ttf"
   );
 
   font_bold = loadFont(
     "/font/Roboto-Bold.ttf"
   );
 }
 
 
 
 function translaCanvas() {
   //const zavannah = merilou;
   if (xxx) easycam.panX(xxx);
   if (yyy) easycam.panY(yyy);
   if (zzz) easycam.zoom(zzz);
   xxx = ![], yyy = ![], zzz = ![];
 }
 
 
 
 function draw() {
   ambientLight(60, 60, 60);
   pointLight(255, 255, 255, 300, 300, 3550);
   background(bgcolorval);
   sfere(sfereJson);
   layer();
   if (xxx || yyy || zzz) translaCanvas();
   setFirstView();
 }
 
 
 function setFirstView(){
   if(setFirstViewStatus){
     setFirstViewStatus=false;
     console.log('***** setFirstView')
     easycam.rotateY(-.80);
     easycam.zoom(-1300);
   }
 }
 
 
 
 
 function keyPressed() {
   //console.log('keyPressed');
   var elementoCliccato = event.target;
   if ( $(elementoCliccato).closest('.position_item').length ) {
     console.log('dentro if');
     return;
   }
   if (keyCode === LEFT_ARROW) xxx = 50; else {
     if (keyCode === RIGHT_ARROW) xxx = -50; else {
       if (keyCode === UP_ARROW) yyy = 50; else keyCode === DOWN_ARROW && (yyy = -50);
     }
   }
 
   if (key === 'c' || key === 'C') {
     let buttons = document.querySelectorAll('.listanodi button[cc-item-hide]');
     buttons.forEach(button => {
       button.setAttribute('cc-item-hide', 'true');
     });
     Object.keys(nodeListPosition).map(function (objectKey, index) {
       showHideSphere(index)
     });
   }
 
   if (key === 's' || key === 'S') {
     keyPressedShowAllNode();
   }
 
 }
 
 
 function windowResized() {
   resizeCanvas(windowWidth, windowHeight), easycam.setViewport([0, 0, windowWidth, windowHeight]);
 }
 
 /*
 function mousePressed() {
   loop();
 }
 function touchStarted() {
   loop();
 }
 
 function mouseReleased() {
   noLoop();
 }
 function touchEnded() {
   noLoop();
 }
 
 function mouseWheel(event) {
   loop();
 }
 */
 