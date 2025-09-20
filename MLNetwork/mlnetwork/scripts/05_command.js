

 /* src/05_command.js - start*/ 
 var comands = {};

 comands.handlers = function () {
   $('.comand-item.comand-item-button.layer').on('click', function (e) {
     e.preventDefault();
     $('.menu-layer.listalayer').toggleClass('selected');
     $('.drawer.drawer_layers').addClass('open');
   });
 
   $('.comand-item.comand-item-button.nodi').on('click', function (e) {
     e.preventDefault();
     $('.menu-layer.listanodi').toggleClass('selected');
     $('.drawer.drawer_nodi').addClass('open');
   });
 
   $('.drawer .close a').on('click', function (e) {
     e.preventDefault();
     $(this).parents('.drawer').removeClass('open');
   });
 
   $('.menu-collapse button').on('click', function (e) {
     e.preventDefault();
     $('body').toggleClass('menu-collapsed');
   });
   $('.comandi-collapse button').on('click', function (e) {
     e.preventDefault();
     $('.comandi-bottom').toggleClass('collapsed');
   });
 
 
   $('body').on('click', '.input-minus-plus.layers button.meno', function () {
     console.log('layer meno')
     layerZcoef -= .1;
   })
 
 
   $('body').on('click', '.input-minus-plus.layers button.piu', function () {
     console.log('layer piu')
     layerZcoef += .1;
   })
 
 
   $('body').on('click', '.input-minus-plus.nodi button.meno', function () {
     console.log('nodi meno')
     sfereCoef -= .1;
   })
 
   $('body').on('click', '.input-minus-plus.nodi button.piu', function () {
     console.log('nodi piu')
     sfereCoef += .1;
   })
 
 
   $('body').on('click', '.input-minus-plus.edge button.meno', function () {
     console.log('edge meno')
     lineeCoef -= .5;
   })
 
   $('body').on('click', '.position_button', function(){
     $(this).parents().find('.textPopup').attr('xxx-label-position',  $(this).attr('xxx-data-position'));
 })
 
   $('body').on('click', '.style_button', function(){
     $(this).parents().find('.textPopup').attr('xxx-label-style',  $(this).attr('xxx-data-style'))
   });
 
   $('body').on('click', '.textPopup input, .textPopup .control_button', function(){
     $('.position_contaniner').attr('data-changed', true);
   })
 
   $('body').on('click', '.input-minus-plus.edge button.piu', function () {
     console.log('edge piu')
     lineeCoef += .5;
   })
 
 
   $('body').on('click', '.layers-details-item .manage-show', function(){
     let index = $(this).attr('xx-data-id');
     showHideLayer(index)
   });
 
   $(".item_name input").on('change keydown paste input', function(){
     var value = $(this).val();
     $('.textPopup').attr('xxx-display-name', value);
   });
 
   $("input#labelAdjustY").on('change keydown paste input', function(){
     var value = $(this).val();
     $('.textPopup').attr('xxx-axis-y', value);
   });
 
   $("input#labelAdjustX").on('change keydown paste input', function(){
     var value = $(this).val();
     $('.textPopup').attr('xxx-axis-x', value);
   });
 
 
   $('body').on('click', 'button.start-recording', function () {
     $('.screen-recording').attr('xxx-video-recorder', 'true')
   })
 
   $('body').on('click', 'button.stop-recording', function () {
     $('.screen-recording').attr('xxx-video-recorder', 'false')
   })
 
   $('body').on('click', '.closeVideoBox', function () {
     $('#videoContainer').remove();
   })
 
 
 
   colorPicker = document.querySelector('#canvasColor');
   colorPicker.addEventListener("change", watchColorPicker, false);
   function watchColorPicker(event) {
     bgcolorval = event.target.value;
     console.log('color: ' + bgcolorval)
   }
 
 
 
 
 
 }
 
 comands.initial = function () {
   comands.handlers();
 }
 
 $(function () {
   comands.initial();
 })